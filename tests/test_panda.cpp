#include <catch2/catch_test_macros.hpp>
#include "se3pf/panda.hpp"
#include "se3pf/se3.hpp"
#include "test_helpers.hpp"

#include <chrono>
#include <cmath>
#include <limits>

using namespace Eigen;
using namespace se3pf;
using testing::maxDiff;
using testing::sci;

namespace {
    // Random configuration inside the joint limits, 0.1 rad away from them.
    Vec7 randomConfig() { return testing::randomConfig(0.1); }

    // Reference values from the earlier implementation (~/Desktop/panda_path_following, panda::bodyJacobian),
    // printed with 15 decimals, at three configurations inside the limits.
    struct ReferenceValue {
        Vec7 q;
        Matrix<double, 3, 4> g;  // top three rows of the flange pose
        Mat67 Jb;
    };

    // clang-format off
    ReferenceValue reference(int k)
    {
        ReferenceValue r;
        if (k == 0) {
            r.q << 0.3, -0.5, 0.7, -1.9, 0.4, 1.8, -0.6;
            r.g << 0.006495809377131, 0.995159515878550, 0.098057852397989, 0.149833601876031,
                   0.958396822475364, -0.034179175684509, 0.283385452376535, 0.444551035812734,
                   0.285365266158413, 0.092137516278101, -0.953979739284038, 0.666388305172192;
            r.Jb << 0.285365266158413, 0.913671812643188, 0.111670819296981, -0.447047216669530, 0.803751133259188, -0.564642473395035, 0.0,
                    0.092137516278101, -0.326742359493074, -0.370094926647863, 0.810142123736426, 0.549875735008960, 0.825335614909678, 0.0,
                    -0.953979739284038, 0.241750386337250, -0.922257758648337, -0.379234130077793, 0.227202094693087, 0.0, 1.0,
                    0.140712329150981, 0.018155876477127, 0.217477922940341, 0.189976819374921, 0.047547363807651, 0.088310910795335, 0.0,
                    -0.447520382584671, 0.288295734262497, -0.461765094246760, -0.103118376483897, -0.069499788971892, 0.060416744653269, 0.0,
                    -0.001131056804254, 0.321033265135990, 0.211635906217381, -0.444235197975579, 0.0, -0.088, 0.0;
        } else if (k == 1) {
            r.q << -1.2, 0.9, -2.1, -0.4, 2.0, 0.3, 1.5;
            r.g << -0.974082878427516, -0.181287591236299, 0.135267716836299, 0.080236165827568,
                   -0.106284445876585, 0.894735076580754, 0.433765788531945, -0.548344612238793,
                   -0.199665125946100, 0.408146973527642, -0.890814170004695, 0.698578990636545;
            r.Jb << -0.199665125946100, -0.946396308809544, -0.323004629184119, -0.476552729977983, 0.020904272455485, 0.997494986604054, 0.0,
                    0.408146973527642, 0.155247072368236, -0.450987565337682, 0.837071881024182, -0.294779924584880, 0.070737201667703, 0.0,
                    -0.890814170004695, 0.283253196262217, -0.832032586763535, -0.268715763492150, -0.955336489125606, 0.0, 1.0,
                    -0.542660954684029, 0.015027261952694, -0.046925492814079, -0.273242234968413, -0.115400467471781, 0.007568880578444, 0.0,
                    -0.027617961943898, -0.549343549860207, -0.101242125729398, -0.103342275293434, -0.008183606183215, -0.106731963566634, 0.0,
                    0.108977027476002, 0.351296029124413, 0.073093400625774, 0.162660424837955, 0.0, -0.088, 0.0;
        } else {
            r.q << 2.5, -1.5, 1.1, -2.9, -2.4, 3.3, -2.7;
            r.g << -0.056832067446435, 0.469810021658104, -0.880936240405272, -0.243994648288803,
                   0.125110075546662, 0.878758262316029, 0.460577230666065, -0.082273914264048,
                   0.990513798544206, -0.084038443351390, -0.108719616135281, 0.351917753174064;
            r.Jb << 0.990513798544206, -0.066218729006274, -0.050037751333651, 0.918168985413495, 0.142613487598055, -0.427379880233830, 0.0,
                    -0.084038443351390, -0.985179782446856, -0.155097119097840, 0.381592074888088, -0.067417335870344, -0.904072142017061, 0.0,
                    -0.108719616135281, 0.158227292799817, -0.986631191017709, -0.106551408284247, 0.987479769908865, 0.0, 1.0,
                    -0.035201985524930, -0.142571081296099, -0.213140534435541, 0.175010397176213, -0.044352205688261, -0.096735719195826, 0.0,
                    -0.175759203702390, 0.015118170621797, -0.013257157586698, -0.438202780253393, -0.093821902841634, 0.045729647185020, 0.0,
                    -0.184856452121393, 0.034464599319146, 0.012893591978587, -0.061243576519582, 0.0, -0.088, 0.0;
        }
        return r;
    }
    // clang-format on
}

TEST_CASE("Panda FK at q = 0: arm straight up, flange at (0.088, 0, 0.926), R = diag(1, -1, -1)")
{
    const Pose g = PandaArm::output(Vec7::Zero(), nullptr);
    const double errR = maxDiff(rotation(g), Matrix3d(Vector3d(1.0, -1.0, -1.0).asDiagonal()));
    const double errP = maxDiff(translation(g), Vector3d(0.088, 0.0, 0.926));
    CAPTURE(sci(errR), sci(errP));
    CHECK(errR < 1e-12);
    CHECK(errP < 1e-12);
    CHECK(maxDiff(Vector4d(g.row(3).transpose()), Vector4d(0, 0, 0, 1)) == 0.0);
}

TEST_CASE("Panda FK at the rest posture: flange at (0.306891, 0, 0.590282), approach axis down")
{
    const Pose g = PandaArm::output(PandaArm::restPosture(), nullptr);
    const double c = std::sqrt(0.5);
    Matrix3d R;
    R <<  c, -c,  0,
         -c, -c,  0,
          0,  0, -1;
    const double errR = maxDiff(rotation(g), R);
    const double errP = maxDiff(translation(g), Vector3d(0.306891, 0.0, 0.590282));
    CAPTURE(sci(errR), sci(errP));
    // and it is the attitude the path tests use for R0
    CHECK(maxDiff(rotation(g), testing::restAttitude()) < 1e-12);
    CHECK(errR < 1e-6);
    CHECK(errP < 1e-6);
}

TEST_CASE("Panda body Jacobian: column j against a central difference of FK along joint j")
{
    const double h = 1e-6;
    for (int k = 0; k < 20; ++k) {
        const Vec7 q = randomConfig();
        Mat67 Jb = Mat67::Zero();
        const Pose g = PandaArm::output(q, &Jb);

        // the pose does not depend on whether the Jacobian is asked for; R is a rotation
        CHECK(maxDiff(PandaArm::output(q, nullptr), g) == 0.0);
        CHECK(maxDiff(Matrix3d(rotation(g).transpose() * rotation(g)), Matrix3d::Identity()) < 1e-14);

        for (int j = 0; j < 7; ++j) {
            const Vec7 e = h * Vec7::Unit(j);
            const Matrix4d dg =
                (PandaArm::output(Vec7(q + e), nullptr) - PandaArm::output(Vec7(q - e), nullptr)) / (2 * h);
            const Vec6 fd = vee(Matrix4d(inverse(g) * dg));  // (g^-1 dg/dq_j)^vee, ordered (omega, v)
            const double err = maxDiff(fd, Vec6(Jb.col(j)));
            CAPTURE(k, j, sci(err), q.transpose());
            CHECK(err < 1e-8);
        }
    }
}

TEST_CASE("Panda body Jacobian: rows (omega, v) in the flange frame, two columns by hand")
{
    // Joint 7 spins the flange about its own z-axis, through its origin: (omega, v) = (e3, 0) at every q.
    for (int k = 0; k < 5; ++k) {
        Mat67 Jb = Mat67::Zero();
        PandaArm::output(randomConfig(), &Jb);
        Vec6 spin;
        spin << 0, 0, 1, 0, 0, 0;
        CHECK(maxDiff(Vec6(Jb.col(6)), spin) < 1e-14);
    }

    // At q = 0, joint 1 turns about the base z-axis. The flange (R = diag(1, -1, -1), p = (0.088, 0, 0.926))
    // turns about -e3 in its own frame, and its origin moves with e3 x p = (0, 0.088, 0) in the base frame,
    // which is (0, -0.088, 0) in the flange frame.
    Mat67 Jb = Mat67::Zero();
    PandaArm::output(Vec7::Zero(), &Jb);
    Vec6 col0;
    col0 << 0, 0, -1, 0, -0.088, 0;
    const double err = maxDiff(Vec6(Jb.col(0)), col0);
    CAPTURE(sci(err));
    CHECK(err < 1e-15);
}

TEST_CASE("Panda FK and body Jacobian agree with the reference implementation")
{
    for (int k = 0; k < 3; ++k) {
        const ReferenceValue r = reference(k);
        Mat67 Jb = Mat67::Zero();
        const Pose g = PandaArm::output(r.q, &Jb);
        const double errG = maxDiff(Matrix<double, 3, 4>(g.topRows<3>()), r.g);
        const double errJ = maxDiff(Jb, r.Jb);
        CAPTURE(k, sci(errG), sci(errJ));
        CHECK(errG < 1e-12);
        CHECK(errJ < 1e-12);
    }
}

TEST_CASE("Panda jointMargin: distance to the nearest position limit, and which joint")
{
    const double pi = EIGEN_PI;
    int joint = -1;

    // At the rest posture the nearest limit is joint 4's lower one: q_4 = -3pi/4, q_min = -3.0718.
    const Vec7 qr = PandaArm::restPosture();
    CHECK(std::abs(PandaArm::jointMargin(qr, &joint) - (3.0718 - 3 * pi / 4)) < 1e-15);
    CHECK(joint == 3);
    CHECK(PandaArm::jointMargin(qr) == PandaArm::jointMargin(qr, &joint));  // joint may be null

    // 0.01 rad below joint 6's upper limit
    Vec7 q = qr;
    q[5] = PandaArm::qMax()[5] - 0.01;
    CHECK(std::abs(PandaArm::jointMargin(q, &joint) - 0.01) < 1e-14);
    CHECK(joint == 5);

    // 0.02 rad past joint 2's lower limit: negative
    q = qr;
    q[1] = PandaArm::qMin()[1] - 0.02;
    CHECK(std::abs(PandaArm::jointMargin(q, &joint) + 0.02) < 1e-14);
    CHECK(joint == 1);
}

TEST_CASE("Panda jointMargin fails closed: NaN if any joint is NaN, -inf if a joint is infinite")
{
    const Vec7 qr = PandaArm::restPosture();
    for (int i = 0; i < 7; ++i) {
        Vec7 q = qr;
        q[i] = std::numeric_limits<double>::quiet_NaN();
        int joint = -1;
        const double margin = PandaArm::jointMargin(q, &joint);
        CAPTURE(i, margin, joint);
        CHECK(std::isnan(margin));
        CHECK_FALSE(margin >= 0.10);  // what a safety check would test
        CHECK(joint == i);

        q[i] = std::numeric_limits<double>::infinity();
        CHECK(PandaArm::jointMargin(q, &joint) == -std::numeric_limits<double>::infinity());
        CHECK(joint == i);
    }
}

TEST_CASE("Panda retract is q + h dq")
{
    const Vec7 q = randomConfig(), dq = testing::randn<Vec7>();
    CHECK(maxDiff(PandaArm::retract(q, dq, 0.0), q) == 0.0);
    CHECK(maxDiff(PandaArm::retract(q, dq, 1e-3), Vec7(q + 1e-3 * dq)) == 0.0);
    CHECK(maxDiff(PandaArm::retract(q, dq, -2.0), Vec7(q - 2.0 * dq)) == 0.0);
}

TEST_CASE("Panda limits: the table of the design doc, entry by entry; the rest posture is well inside")
{
    // The appendix table "Panda limits", typed once more here so that a typo in either copy fails the test.
    // Literals compare exactly, so ==.
    //                              joint 1  joint 2  joint 3  joint 4  joint 5  joint 6  joint 7
    const Vec7 qMin    = (Vec7() << -2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973).finished();
    const Vec7 qMax    = (Vec7() <<  2.8973,  1.7628,  2.8973, -0.0698,  2.8973,  3.7525,  2.8973).finished();
    const Vec7 dqMax   = (Vec7() <<  2.1750,  2.1750,  2.1750,  2.1750,  2.6100,  2.6100,  2.6100).finished();
    const Vec7 ddqMax  = (Vec7() <<  15.0,    7.5,     10.0,    12.5,    15.0,    20.0,    20.0).finished();
    const Vec7 jerkMax = (Vec7() <<  7500.0,  3750.0,  5000.0,  6250.0,  7500.0,  10000.0, 10000.0).finished();
    const Vec7 tauMax  = (Vec7() <<  87.0,    87.0,    87.0,    87.0,    12.0,    12.0,    12.0).finished();

    for (int i = 0; i < 7; ++i) {
        CAPTURE(i);
        CHECK(PandaArm::qMin()[i] == qMin[i]);
        CHECK(PandaArm::qMax()[i] == qMax[i]);
        CHECK(PandaArm::dqMax()[i] == dqMax[i]);
        CHECK(PandaArm::ddqMax()[i] == ddqMax[i]);
        CHECK(PandaArm::jerkMax()[i] == jerkMax[i]);
        CHECK(PandaArm::tauMax()[i] == tauMax[i]);

        // Two facts about Franka's table, as documentation: every range is non-empty, and jerk = 500 x acceleration.
        CHECK(PandaArm::qMin()[i] < PandaArm::qMax()[i]);
        CHECK(PandaArm::jerkMax()[i] == 500.0 * PandaArm::ddqMax()[i]);
    }
    CHECK(PandaArm::tauRateMax == 1000.0);
    CHECK(PandaArm::flangeSpeedMax == 1.7);
    CHECK(PandaArm::kFlangeOffset == 0.107);
    CHECK(PandaArm::jointMargin(PandaArm::restPosture()) > 0.7);

    // The hardware margins of the design doc (Hardware), typed once more
    STATIC_REQUIRE(HardwareBounds::jointMargin == 0.10);
    STATIC_REQUIRE(HardwareBounds::dqFrac == 0.5);
    STATIC_REQUIRE(HardwareBounds::ddqFrac == 0.3);
    STATIC_REQUIRE(HardwareBounds::tauFrac == 0.6);
    STATIC_REQUIRE(HardwareBounds::tauRateFrac * PandaArm::tauRateMax == 600.0);
    STATIC_REQUIRE(HardwareBounds::flangeSpeed == 0.5);
    STATIC_REQUIRE(HardwareBounds::xi == 0.04);
    STATIC_REQUIRE(HardwareBounds::beta == 0.05);
    STATIC_REQUIRE(HardwareBounds::sigmaMin == 1e-3);
}

TEST_CASE("Timing: PandaArm::output(q, &Jb)", "[.timing]")
{
    // Release build only: ./build-rel/tests "[.timing]". The controller calls output three times per cycle.
    // The name does not start with "Panda", so the filter "Panda*" selects only the correctness tests.
#ifndef NDEBUG
    SKIP("timing is meaningful in a Release build only");
#endif
    const int calls = 1000000;
    Vec7 q = PandaArm::restPosture();
    Mat67 Jb = Mat67::Zero();
    double sink = 0.0;  // uses every result, so the optimizer cannot drop the calls

    const auto start = std::chrono::steady_clock::now();
    for (int k = 0; k < calls; ++k) {
        q[k % 7] += 1e-7;
        const Pose g = PandaArm::output(q, &Jb);
        sink += g(0, 3) + Jb(5, 3);
    }
    const auto stop = std::chrono::steady_clock::now();

    const double us = std::chrono::duration<double, std::micro>(stop - start).count() / calls;
    WARN("PandaArm::output(q, &Jb): " << us << " us per call (checksum " << sink << ")");
    CHECK(std::isfinite(sink));
    CHECK(us < 2.0);
}
