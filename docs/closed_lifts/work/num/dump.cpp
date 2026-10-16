// Dump what the harmonic-lift experiments need from the project library (se3pf) to CSV:
//   data/rest_N.csv     for N in {250, 500, 1000, 2000}: s_k, gamma(s_k) (R row-major, p), zeta, zeta',
//                       the rest-posture lift node q*_rest(s_k) (Lift(path, q_r, N)), x_k = kernelVector(J_b(q*)),
//                       sigma_min(J_b(q*)), joint margin
//   data/frame_N.csv    N(s_k), N'(s_k) from Tube::frame (column-major 6x5 each) at s_k = k L / N
//   data/fk_check.csv   random q, FK(q) (R, p), J_b(q) (row-major) to validate the Python kinematics
//   data/meta.csv       L, sigma, holonomy angles, tube closure, lift diagnostics for N = 2000
#include "se3pf/panda_example.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/se3.hpp"

#include <cstdio>
#include <random>
#include <string>

using namespace se3pf;

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : "data";
    const PandaExample ex;  // l = 0.1 m: W = diag(0.01 I, I)
    const Path& path = ex.path();
    const double L = path.length();

    {
        FILE* f = std::fopen((dir + "/meta.csv").c_str(), "w");
        const auto& td = ex.tube().diagnostics();
        const auto& ld = ex.lift().diagnostics();
        std::fprintf(f, "key,value\n");
        std::fprintf(f, "L,%.17g\n", L);
        std::fprintf(f, "planes,%d\n", td.planes);
        std::fprintf(f, "theta1,%.17g\ntheta2,%.17g\n", td.angles[0], td.angles[1]);
        std::fprintf(f, "frame_closure,%.3e\n", td.closureGap);
        std::fprintf(f, "lift2000_closure,%.3e\nlift2000_minMargin,%.17g\nlift2000_minMarginJoint,%d\n",
                     ld.closureGap, ld.minJointMargin, ld.minMarginJoint);
        std::fprintf(f, "lift2000_minSigma,%.17g\nlift2000_maxDerivative,%.17g\n", ld.minSigma, ld.maxDerivative);
        const Matrix3d R0 = ex.R0();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                std::fprintf(f, "R0_%d%d,%.17g\n", i, j, R0(i, j));
        std::fclose(f);
    }

    for (int N : {250, 500, 1000, 2000}) {
        const Lift lift(path, PandaArm::restPosture(), N);
        const auto& nodes = lift.nodes();
        FILE* f = std::fopen((dir + "/rest_" + std::to_string(N) + ".csv").c_str(), "w");
        std::fprintf(f, "k,s");
        for (int i = 0; i < 9; ++i) std::fprintf(f, ",R%d", i);
        for (int i = 0; i < 3; ++i) std::fprintf(f, ",p%d", i);
        for (int i = 0; i < 6; ++i) std::fprintf(f, ",z%d", i);
        for (int i = 0; i < 6; ++i) std::fprintf(f, ",dz%d", i);
        for (int i = 0; i < 7; ++i) std::fprintf(f, ",q%d", i);
        for (int i = 0; i < 7; ++i) std::fprintf(f, ",x%d", i);
        std::fprintf(f, ",sigmin,margin,closure\n");
        for (int k = 0; k < N; ++k) {
            const double s = k * L / N;
            const Pose g = path.pose(s);
            const Vec6 z = path.velocity(s), dz = path.acceleration(s);
            Mat67 J;
            PandaArm::output(nodes[k], &J);
            const Vec7 x = kernelVector(J);
            Eigen::JacobiSVD<Mat67> svd(J);
            std::fprintf(f, "%d,%.17g", k, s);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) std::fprintf(f, ",%.17g", g(i, j));
            for (int i = 0; i < 3; ++i) std::fprintf(f, ",%.17g", g(i, 3));
            for (int i = 0; i < 6; ++i) std::fprintf(f, ",%.17g", z(i));
            for (int i = 0; i < 6; ++i) std::fprintf(f, ",%.17g", dz(i));
            for (int i = 0; i < 7; ++i) std::fprintf(f, ",%.17g", nodes[k](i));
            for (int i = 0; i < 7; ++i) std::fprintf(f, ",%.17g", x(i));
            std::fprintf(f, ",%.17g,%.17g,%.3e\n", svd.singularValues().minCoeff(), PandaArm::jointMargin(nodes[k]),
                         lift.diagnostics().closureGap);
        }
        std::fclose(f);

        FILE* g = std::fopen((dir + "/frame_" + std::to_string(N) + ".csv").c_str(), "w");
        std::fprintf(g, "k,s");
        for (int i = 0; i < 30; ++i) std::fprintf(g, ",N%d", i);
        for (int i = 0; i < 30; ++i) std::fprintf(g, ",dN%d", i);
        std::fprintf(g, "\n");
        for (int k = 0; k < N; ++k) {
            const double s = k * L / N;
            Mat65 Nf, dN;
            ex.tube().frame(s, Nf, &dN);
            std::fprintf(g, "%d,%.17g", k, s);
            for (int i = 0; i < 30; ++i) std::fprintf(g, ",%.17g", Nf.data()[i]);
            for (int i = 0; i < 30; ++i) std::fprintf(g, ",%.17g", dN.data()[i]);
            std::fprintf(g, "\n");
        }
        std::fclose(g);
    }

    {
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        FILE* f = std::fopen((dir + "/fk_check.csv").c_str(), "w");
        for (int t = 0; t < 20; ++t) {
            Vec7 q = PandaArm::restPosture();
            for (int i = 0; i < 7; ++i) q(i) += 0.8 * U(rng);
            Mat67 J;
            const Pose g = PandaArm::output(q, &J);
            for (int i = 0; i < 7; ++i) std::fprintf(f, "%.17g,", q(i));
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 4; ++j) std::fprintf(f, "%.17g,", g(i, j));
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 7; ++j) std::fprintf(f, "%.17g%s", J(i, j), (i == 5 && j == 6) ? "\n" : ",");
        }
        std::fclose(f);
    }
    std::printf("L = %.10f, holonomy angles %.6g %.6g, frame closure %.2e, lift(2000) closure %.2e\n", L,
                ex.tube().diagnostics().angles[0], ex.tube().diagnostics().angles[1],
                ex.tube().diagnostics().closureGap, ex.lift().diagnostics().closureGap);
}
