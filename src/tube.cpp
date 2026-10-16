#include "se3pf/tube.hpp"
#include "se3pf/se3.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace Eigen;

namespace se3pf {

    namespace {
        using RowVec5 = Eigen::Matrix<double, 1, 5>;

        // Right-hand side of eq:transport, N0' = -zeta_gamma zeta_gamma'^T W N0. N0 moves only along zeta_gamma,
        // just enough to stay W-orthogonal to it: with |zeta_gamma|_W = 1, (N0^T W zeta_gamma)' = 0 and
        // (N0^T W N0)' = 0, so N0 stays a W-orthonormal frame of the normal space.
        Mat65 transportRate(const Path& path, const Metric& W, double s, const Mat65& N)
        {
            const Vec6 z = path.velocity(s);
            const Vec6 dz = path.acceleration(s);
            const RowVec5 row = W.flat(dz).transpose() * N;  // zeta_gamma'^T W N
            return -z * row;
        }

        // Step 1 of Lemma 1: Gram-Schmidt in <.,.>_W on zeta_gamma(0), e_1, ..., e_6, in that order, skipping
        // every vector whose residual has W-norm below 1e-8. The five vectors after zeta_gamma are N0(0). The
        // order fixes the frame, and with it the meaning of xi (the start poses Psi(0.1, xi0) assume this order).
        Mat65 initialFrame(const Vec6& z0, const Metric& W)
        {
            std::vector<Vec6> basis;
            auto add = [&](Vec6 x) {
                for (const Vec6& b : basis)
                    x -= W.inner(b, x) * b;
                const double nx = W.norm(x);
                if (nx > 1e-8)
                    basis.push_back(x / nx);
            };
            add(z0);
            for (int i = 0; i < 6; ++i)
                add(Vec6::Unit(i));
            // seven vectors in R^6: exactly one e_i drops out (e_5 on the rigid-body path, e_6 on the Panda path)
            if (basis.size() != 6)
                throw std::runtime_error("Tube: Gram-Schmidt did not give six W-orthonormal vectors");

            Mat65 N0 = Mat65::Zero();
            for (int j = 0; j < 5; ++j)
                N0.col(j) = basis[j + 1];
            return N0;
        }
    }  // namespace

    // ---------------------------------------------------------------------------------------------------------
    // Rotations in SO(5) by invariant planes

    Mat5 RotationPlanes::log() const
    {
        Mat5 Phi = Mat5::Zero();
        for (int i = 0; i < count; ++i)
            Phi += theta[i] * (v[i] * u[i].transpose() - u[i] * v[i].transpose());
        return Phi;
    }

    Mat5 RotationPlanes::exp(double t) const
    {
        Mat5 E = Mat5::Identity();
        for (int i = 0; i < count; ++i) {
            const double c = std::cos(t * theta[i]), s = std::sin(t * theta[i]);
            E += (c - 1.0) * (u[i] * u[i].transpose() + v[i] * v[i].transpose())
               + s * (v[i] * u[i].transpose() - u[i] * v[i].transpose());
        }
        return E;
    }

    RotationPlanes rotationPlanes(const Mat5& H)
    {
        if (!(H.determinant() > 0.0))
            throw std::invalid_argument("rotationPlanes: H is not a rotation (det H <= 0)");

        // (H + H^T)/2 = sum_i cos(theta_i) (u_i u_i^T + v_i v_i^T) + (projection onto the fixed directions)
        const Mat5 S = 0.5 * (H + H.transpose());
        const SelfAdjointEigenSolver<Mat5> eig(S);  // eigenvalues in increasing order: largest angles first
        const Vec5& lambda = eig.eigenvalues();
        const Mat5& V = eig.eigenvectors();

        RotationPlanes P;
        std::array<Vec5, 4> used = {Vec5::Zero(), Vec5::Zero(), Vec5::Zero(), Vec5::Zero()};
        int nUsed = 0;
        // x minus its components along the plane vectors found so far
        auto orthogonalize = [&](Vec5 x) {
            for (int j = 0; j < nUsed; ++j)
                x -= used[j].dot(x) * used[j];
            return x;
        };

        for (int k = 0; k < 5; ++k) {
            // A fixed direction only if 1 - cos(theta) < 1e-12: the Panda path's second plane has theta = 8.0e-5,
            // 1 - cos(theta) = 3.2e-9, and dropping it would leave the frame open by about 1e-4.
            if (1.0 - lambda(k) < 1e-12)
                continue;
            Vec5 u = orthogonalize(V.col(k));
            if (u.norm() < 0.5)  // the second eigenvector of a plane found already
                continue;
            u.normalize();

            const Vec5 Hu = H * u;
            const double c = u.dot(Hu);
            const Vec5 r = Hu - c * u;  // sin(theta) v
            const double s = r.norm();
            Vec5 v = Vec5::Zero();
            if (s > 1e-9) {
                v = r / s;
            } else {
                // A rotation by pi (to within sin(theta) < 1e-9): H u = -u, and r has no reliable direction. Every
                // unit vector of the (-1)-eigenspace orthogonal to u (and to the planes found so far) completes the
                // plane, and each choice gives a logarithm of H (to within 2 sin(theta)). The frame does not need a
                // particular one: N0(s) exp(-(s/L) Phi) closes for any Phi with exp(Phi) = H (eq:holonomy), so take
                // the first such eigenvector.
                bool found = false;
                for (int j = 0; j < 5 && !found; ++j) {
                    if (j == k || lambda(j) > -1.0 + 1e-6)
                        continue;
                    Vec5 x = orthogonalize(V.col(j));
                    x -= u.dot(x) * u;
                    if (x.norm() > 0.5) {
                        v = x.normalized();
                        found = true;
                    }
                }
                if (!found)
                    throw std::runtime_error("rotationPlanes: no partner for a rotation by pi (H is not orthogonal)");
            }
            if (P.count == 2)
                throw std::runtime_error("rotationPlanes: more than two planes found (H is not orthogonal)");

            P.u[P.count] = u;
            P.v[P.count] = v;
            P.theta[P.count] = std::atan2(s, c);
            used[nUsed++] = P.u[P.count];
            used[nUsed++] = P.v[P.count];
            ++P.count;
        }
        return P;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Offline: the normal frame (Lemma 1)

    Tube::Tube(const Path& path, const Metric& metric, int gridIntervals)
        : path_(path), metric_(metric), L_(path.length()), n_(gridIntervals)
    {
        if (n_ < 4)
            throw std::invalid_argument("Tube: need at least 4 grid intervals");
        const double h = L_ / n_;

        // Step 1: initial frame N0(0) by Gram-Schmidt.
        std::vector<Mat65> N0(n_ + 1, Mat65::Zero());
        N0[0] = initialFrame(path_.velocity(0.0), metric_);

        // Step 2: transport (eq:transport) by RK4, one step per grid interval.
        for (int k = 0; k < n_; ++k) {
            const double s = k * h;
            const Mat65& N = N0[k];
            const Mat65 k1 = transportRate(path_, metric_, s, N);
            const Mat65 k2 = transportRate(path_, metric_, s + 0.5 * h, N + (0.5 * h) * k1);
            const Mat65 k3 = transportRate(path_, metric_, s + 0.5 * h, N + (0.5 * h) * k2);
            const Mat65 k4 = transportRate(path_, metric_, s + h, N + h * k3);
            N0[k + 1] = N + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }

        // Step 3: holonomy H = N0(0)^T W N0(L), so that N0(L) = N0(0) H. It is a rotation of the normal space.
        const Mat5 H = N0[0].transpose() * metric_.matrix() * N0[n_];
        diag_.detH = H.determinant();
        diag_.orthogonalityError = (H.transpose() * H - Mat5::Identity()).cwiseAbs().maxCoeff();
        // A loose guard (the default grid gives about 1e-11, which the tests check to 1e-10): only a grossly
        // inaccurate transport, e.g. far too few grid intervals, stops the construction here.
        if (std::abs(diag_.detH - 1.0) > 1e-6 || diag_.orthogonalityError > 1e-6)
            throw std::runtime_error("Tube: the holonomy is not a rotation; the transport is too inaccurate");

        // Step 4: Phi = log H. RK4 leaves H orthogonal only to about 1e-11, so first replace it by the nearest
        // rotation H_r = U V^T (from the SVD H = U Sigma V^T), then take the logarithm from its invariant planes.
        const JacobiSVD<Mat5> svd(H, ComputeFullU | ComputeFullV);
        const Mat5 Hr = svd.matrixU() * svd.matrixV().transpose();
        const RotationPlanes planes = rotationPlanes(Hr);
        const Mat5 Phi = planes.log();

        diag_.planes = planes.count;
        for (int i = 0; i < planes.count; ++i)
            diag_.angles[i] = planes.theta[i];
        std::sort(diag_.angles.begin(), diag_.angles.end(), std::greater<double>());
        diag_.logError = (planes.exp(1.0) - Hr).cwiseAbs().maxCoeff();
        diag_.skewError = (Phi + Phi.transpose()).cwiseAbs().maxCoeff();

        // Step 5: the closed frame N(s) = N0(s) exp(-(s/L) Phi) (eq:holonomy), and its exact derivative from
        // eq:frame-ode, N' = -zeta_gamma zeta_gamma'^T W N - (1/L) N Phi, at the nodes.
        std::vector<Mat65> Y(n_ + 1, Mat65::Zero()), D(n_ + 1, Mat65::Zero());
        for (int k = 0; k <= n_; ++k) {
            const double s = k * h;
            Y[k] = N0[k] * planes.exp(-s / L_);
            D[k] = transportRate(path_, metric_, s, Y[k]) - (1.0 / L_) * (Y[k] * Phi);
        }
        // N(L) must equal N(0). Measured here, on the node values: the spline replaces node n by node 0.
        diag_.closureGap = (Y[n_] - Y[0]).cwiseAbs().maxCoeff();

        N_.build(L_, std::move(Y), std::move(D));
    }

    double Tube::wrap(double s) const
    {
        double x = std::fmod(s, L_);  // in (-L, L)
        if (x < 0.0)
            x += L_;
        if (x >= L_)  // -1e-17 + L rounds to L
            x -= L_;
        return x;
    }

    void Tube::frame(double s, Mat65& N, Mat65* dN) const
    {
        N_.eval(s, N, dN);
    }

    // ---------------------------------------------------------------------------------------------------------
    // Online: the tubular map, the projection and J_Psi

    Pose Tube::Psi(double s, const Vec5& xi) const
    {
        Mat65 N;
        frame(s, N);
        const Vec6 e = N * xi;
        return path_.pose(s) * metric_.Exp(e);
    }

    bool Tube::projectionCondition(const Pose& g, double s, double& phi, double& dphi, Vec6& e) const
    {
        const Pose gErr = inverse(path_.pose(s)) * g;  // g_err(s) = gamma(s)^-1 g
        // Check if log is defined
        if (!metric_.Log(gErr, e))
            return false;

        const Vec6 z = path_.velocity(s);
        const Vec6 dz = path_.acceleration(s);

        // e_g'(s) = -Gamma(e_g)^-1 Ad_{g_err^-1} zeta_gamma, so phi_g' = <e_g', zeta_gamma>_W + <e_g, zeta_gamma'>_W
        const Vec6 de = -(metric_.GammaInv(e) * (Ad(inverse(gErr)) * z));
        phi = metric_.inner(e, z);
        dphi = metric_.inner(de, z) + metric_.inner(e, dz);

        return std::isfinite(phi) && std::isfinite(dphi);
    }

    bool Tube::newton(const Pose& g, double guess, double& eta, Vec6& e, int* iterations) const
    {
        double s = wrap(guess);
        double phi = 0.0, dphi = 0.0;

        if (!projectionCondition(g, s, phi, dphi, e))
            return false;
        
        // Newton iteration of at most 20 steps
        for (int k = 1; k <= 20; ++k) {
            // stop of phi' gets too small
            if (std::abs(dphi) < 1e-12)
                return false;
            s = wrap(s - phi / dphi); // Newton step
            // check the new s satisfies the projection condition still
            if (!projectionCondition(g, s, phi, dphi, e))
                return false;
            // check early convergence
            if (std::abs(phi) < 1e-10) {
                eta = s;
                if (iterations != nullptr)
                    *iterations = k;
                return true;
            }
        }
        return false;
    }

    bool Tube::outputs(const Pose& g, double guess, double& eta, Vec5& xi, int* iterations) const
    {
        Vec6 e = Vec6::Zero();
        double s = 0.0;
        
        // check if newton iteration failed
        if (!newton(g, guess, s, e, iterations))
            return false;

        Mat65 N;
        frame(s, N);
        
        const Vec5 x = N.transpose() * metric_.flat(e);  // xi = N(eta)^T W e_g(eta) (Proposition 1(a))

        // Check that xi is finite in all entries
        if (!x.allFinite())
            return false;

        eta = s;
        xi = x;
        return true;
    }

    bool Tube::globalProject(const Pose& g, double& eta, int samples) const
    {
        // Newton's method alone may converge to another root of phi_g, which is not the closest point.
        // To prevent this, search over samples to find minimizing s_k of dist(gamma(s_k), g).
        double best = std::numeric_limits<double>::infinity(), sBest = 0.0;
        for (int k = 0; k < samples; ++k) {
            const double s = L_ * k / samples;
            const double d = metric_.distance(path_.pose(s), g);
            if (d < best) {
                best = d;
                sBest = s;
            }
        }

        if (!std::isfinite(best))
            return false;

        Vec6 e = Vec6::Zero();
        return newton(g, sBest, eta, e, nullptr);
    }

    Mat6 Tube::JPsi(double eta, const Vec5& xi) const
    {
        Mat65 N, dN;
        frame(eta, N, &dN);
        const Vec6 e = N * xi;
        const Mat6 G = metric_.Gamma(e);

        //   J_Psi = [ Ad_{Exp_I(N xi)^-1} zeta_gamma + Gamma(N xi) N' xi  |  Gamma(N xi) N ]
        Mat6 J = Mat6::Zero();
        J.col(0) = Ad(inverse(metric_.Exp(e))) * path_.velocity(eta) + G * (dN * xi);
        J.rightCols<5>() = G * N;
        return J;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Offline: tube-radius estimate

    namespace {
        double identityDistance(const Pose& g1, const Pose& g2)
        {
            const Pose g = inverse(g1) * g2;
            const Vector3d w = logSO3(Matrix3d(g.block<3,3>(0,0)));
            const Vector3d p = g.block<3,1>(0,3);
            return std::sqrt(w.squaredNorm() + p.squaredNorm());
        }

        // The test for one epsilon: false as soon as a sampled path point is closer than epsilon to some pose
        // Psi(s, epsilon d). That happens beyond a focal point of the normal geodesics (first seen right beside
        // the foot point s, hence the extra samples there), beyond the cut locus, or where distant parts of the
        // path come close.
        bool radiusPasses(const Tube& tube, double eps, const RadiusOptions& opt, const std::vector<Vec5>& dirs,
                          const std::vector<Pose>& samples, double sqrtLambdaMin)
        {
            const double L = tube.length();
            const double h = L / opt.checkPoints;  // spacing of the path samples
            const double limit = eps * (1.0 - opt.relTol);
            const Metric& W = tube.metric();

            // Most path points are far from g, and a lower bound rules them out: W >= lambda_min(W) I, so
            // dist_W >= sqrt(lambda_min) dist_I. Without it, a far point at a relative rotation angle of pi, where
            // Log_I fails and Metric::distance returns NaN, would count as closer. For the points the bound cannot
            // rule out, the distance itself decides, and a failed logarithm there counts as closer.
            auto closer = [&](const Pose& gp, const Pose& g) {
                if (sqrtLambdaMin * identityDistance(gp, g) >= limit)
                    return false;
                return !(W.distance(gp, g) >= limit);
            };

            for (int i = 0; i < opt.footPoints; ++i) {
                const double s = L * (i + 0.5) / opt.footPoints;
                for (const Vec5& d : dirs) {
                    const Pose g = tube.Psi(s, eps * d);
                    for (const Pose& gp : samples)
                        if (closer(gp, g))
                            return false;
                    for (double off : {0.5, 0.25, 0.125})  // beside the foot point
                        if (closer(tube.path().pose(s + off * h), g) || closer(tube.path().pose(s - off * h), g))
                            return false;
                }
            }
            return true;
        }
    }  // namespace

    double estimateTubeRadius(const Tube& tube, const RadiusOptions& opt)
    {
        if (opt.footPoints < 1 || opt.checkPoints < 1 || opt.randomDirections < 0 || opt.bisections < 0
            || !(opt.startFraction > 0.0))
            throw std::invalid_argument("estimateTubeRadius: invalid options");

        // Directions: +-e_1, ..., +-e_5, then random unit vectors (normal entries, normalized; fixed seed).
        std::vector<Vec5> dirs;
        for (int i = 0; i < 5; ++i) {
            dirs.push_back(Vec5::Unit(i));
            dirs.push_back(-Vec5::Unit(i));
        }
        std::mt19937 rng(opt.seed);
        std::normal_distribution<double> normal(0.0, 1.0);
        for (int k = 0; k < opt.randomDirections; ++k) {
            Vec5 d = Vec5::Zero();
            for (int i = 0; i < 5; ++i)
                d[i] = normal(rng);
            dirs.push_back(d.normalized());
        }

        // The path samples do not depend on epsilon: compute them once.
        std::vector<Pose> samples(opt.checkPoints);
        for (int k = 0; k < opt.checkPoints; ++k)
            samples[k] = tube.path().pose(tube.length() * k / opt.checkPoints);

        const SelfAdjointEigenSolver<Mat6> eigW(tube.metric().matrix(), EigenvaluesOnly);
        const double sqrtLambdaMin = std::sqrt(eigW.eigenvalues().minCoeff());

        // Double epsilon from startFraction * L until the test fails, then bisect.
        double lo = 0.0, hi = opt.startFraction * tube.length();
        while (radiusPasses(tube, hi, opt, dirs, samples, sqrtLambdaMin)) {
            lo = hi;
            hi *= 2.0;
            if (hi > 1e3 * tube.length())  // never fails: not a closed path the test can resolve
                throw std::runtime_error("estimateTubeRadius: the test passes for every epsilon");
        }
        for (int k = 0; k < opt.bisections; ++k) {
            const double mid = 0.5 * (lo + hi);
            if (radiusPasses(tube, mid, opt, dirs, samples, sqrtLambdaMin))
                lo = mid;
            else
                hi = mid;
        }
        return lo;
    }

}  // namespace se3pf
