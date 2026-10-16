// Small helpers shared by the tests (namespace testing). Fixtures that build something expensive live in
// panda_fixture.hpp (the paper's Panda example) and tube_fixture.hpp (the three example tubes).
#pragma once
#include "se3pf/panda.hpp"
#include "se3pf/se3.hpp"
#include "se3pf/types.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>

// A closed loop of 10 000-20 000 control steps takes a minute in a Debug build (unoptimized Eigen is 30-100 times
// slower) and a fraction of a second in Release. Tests tagged SLOW_IN_DEBUG are hidden in Debug ("[.slow]"; run
// them with ./build-dbg/tests "[slow]") and run by default in Release, so tools/ci.sh still checks every number of
// the design doc's Closed loops table. Each hidden closed loop has a short stand-in that runs in Debug too.
#ifdef NDEBUG
#define SLOW_IN_DEBUG "[slow]"
#else
#define SLOW_IN_DEBUG "[.slow]"
#endif

namespace testing {

    using se3pf::Matrix3d;
    using se3pf::Pose;
    using se3pf::Vec6;
    using se3pf::Vec7;
    using se3pf::Vector3d;

    // Largest absolute entry of a - b: an absolute tolerance check that also works against zero
    // (Eigen's isApprox is relative and never passes against an exact zero).
    template <class A, class B>
    double maxDiff(const Eigen::MatrixBase<A>& a, const Eigen::MatrixBase<B>& b)
    {
        return (a - b).cwiseAbs().maxCoeff();
    }

    // Reproducible random numbers. tests/rng_listener.cpp reseeds the generator with kRngSeed before every run of a
    // test case (each leaf section), so a test sees the same numbers whether it runs in the full suite, alone, with
    // a name filter or in random order (--order rand): a failure always reproduces.
    constexpr unsigned kRngSeed = 12345;

    inline std::mt19937& rng()
    {
        static std::mt19937 gen(kRngSeed);
        return gen;
    }

    // Random matrix or vector with standard normal entries times scale.
    template <class M>
    M randn(double scale = 1.0)
    {
        std::normal_distribution<double> nd(0.0, scale);
        M m;
        for (Eigen::Index i = 0; i < m.size(); ++i)
            m.data()[i] = nd(rng());
        return m;
    }

    inline double uniform(double a, double b)
    {
        std::uniform_real_distribution<double> u(a, b);
        return u(rng());
    }

    // Catch2 prints a double in fixed notation with 10 decimals, so an error of 3e-12 would show as 0.0.
    // CAPTURE(sci(err)) prints it in scientific notation instead.
    inline std::string sci(double x)
    {
        std::ostringstream os;
        os << std::scientific << std::setprecision(3) << x;
        return os.str();
    }

    // ---------------------------------------------------------------------------------------------------------
    // Poses and twists

    // Random pose: rotation exp(w) and translation with standard normal entries
    inline Pose randomPose()
    {
        return se3pf::makePose(se3pf::expSO3(randn<Vector3d>()), randn<Vector3d>());
    }

    // The twist (omega, v)
    inline Vec6 twist(const Vector3d& w, const Vector3d& v)
    {
        Vec6 z;
        z << w, v;
        return z;
    }

    // Flange attitude of the Panda in its rest posture (approach axis pointing down): R0 = FK(q_r), typed here so
    // that the path and tube tests do not depend on panda.cpp (test_panda.cpp checks FK(q_r) against it).
    inline Matrix3d restAttitude()
    {
        const double c = std::sqrt(0.5);
        Matrix3d R0;
        R0 << c, -c, 0,
             -c, -c, 0,
              0,  0, -1;
        return R0;
    }

    // Angle (degrees) and distance (m) between the attitudes and positions of two poses
    inline void poseError(const Pose& a, const Pose& b, double& degrees, double& metres)
    {
        const Matrix3d R = se3pf::rotation(a).transpose() * se3pf::rotation(b);
        degrees = se3pf::logSO3(R).norm() * 180.0 / EIGEN_PI;
        metres = (se3pf::translation(b) - se3pf::translation(a)).norm();
    }

    // a - b on the circle R/LZ, wrapped into (-L/2, L/2]. Written here again rather than taken from
    // se3pf::wrapSymmetric, so that the tests do not share code with what they check. It accepts eta = s + L as
    // equal to s, so every test that uses it on an eta also checks that eta lies in [0, L).
    inline double arcDiff(double a, double b, double L)
    {
        double d = a - b;
        while (d > 0.5 * L)
            d -= L;
        while (d <= -0.5 * L)
            d += L;
        return d;
    }

    // ---------------------------------------------------------------------------------------------------------
    // The Panda

    // Random configuration inside the joint limits, at least `margin` rad away from them
    inline Vec7 randomConfig(double margin = 0.0)
    {
        const Vec7 lo = se3pf::PandaArm::qMin(), hi = se3pf::PandaArm::qMax();
        Vec7 q = Vec7::Zero();
        for (int i = 0; i < 7; ++i)
            q[i] = uniform(lo[i] + margin, hi[i] - margin);
        return q;
    }

}  // namespace testing
