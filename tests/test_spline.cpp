#include <catch2/catch_test_macros.hpp>
#include "se3pf/spline.hpp"
#include "se3pf/types.hpp"
#include "test_helpers.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace se3pf;
using testing::maxDiff;

namespace {
    const double L = 0.77;
    const double w = 2.0 * EIGEN_PI / L;  // one period over [0, L)
    double f(double s) { return std::sin(w * s) + 0.3 * std::cos(2 * w * s); }
    double df(double s) { return w * std::cos(w * s) - 0.6 * w * std::sin(2 * w * s); }

    // Largest value and derivative errors of the spline of f on n intervals, sampled over several periods.
    void splineErrors(int n, double& ev, double& ed)
    {
        std::vector<double> y(n + 1), d(n + 1);
        for (int k = 0; k <= n; ++k) {
            y[k] = f(k * L / n);
            d[k] = df(k * L / n);
        }
        PeriodicHermite<double> sp;
        sp.build(L, y, d);
        ev = ed = 0.0;
        for (int i = 0; i < 3000; ++i) {
            const double s = -2.0 * L + 5.0 * L * i / 2999.0;  // also negative s and s beyond L
            double v, dv;
            sp.eval(s, v, &dv);
            ev = std::max(ev, std::abs(v - f(s)));
            ed = std::max(ed, std::abs(dv - df(s)));
        }
    }
}

TEST_CASE("PeriodicHermite: fourth-order values, third-order derivatives")
{
    double ev1, ed1, ev2, ed2;
    splineErrors(100, ev1, ed1);
    splineErrors(200, ev2, ed2);
    CHECK(ev1 < 1e-6);
    CHECK(ev1 / ev2 > 14.0);  // about 16 when n doubles
    CHECK(ed1 / ed2 > 7.0);   // about 8
}

TEST_CASE("PeriodicHermite: closes exactly, wraps s, works for Eigen types")
{
    const int n = 100;
    std::vector<Vec7> Y(n + 1), D(n + 1);
    for (int k = 0; k <= n; ++k)
        for (int i = 0; i < 7; ++i) {
            Y[k][i] = std::sin(w * k * L / n + i);
            D[k][i] = w * std::cos(w * k * L / n + i);
        }
    PeriodicHermite<Vec7> sp;
    sp.build(L, Y, D);

    Vec7 a, b, da;
    sp.eval(L, a);
    CHECK(maxDiff(a, Y[0]) == 0.0);  // s = L is node 0, exactly
    sp.eval(0.123, a, &da);
    sp.eval(0.123 + 3 * L, b);
    CHECK(maxDiff(a, b) < 1e-12);  // periodic
    sp.eval(-L + 0.123, b);
    CHECK(maxDiff(a, b) < 1e-12);
    CHECK(std::abs(a[3] - std::sin(w * 0.123 + 3)) < 1e-6);
    CHECK(std::abs(da[3] - w * std::cos(w * 0.123 + 3)) < 1e-3);

    // a matrix type, as for the normal frame N
    std::vector<Mat65> M(n + 1), DM(n + 1);
    for (int k = 0; k <= n; ++k) {
        M[k] = Mat65::Constant(std::sin(w * k * L / n));
        DM[k] = Mat65::Constant(w * std::cos(w * k * L / n));
    }
    PeriodicHermite<Mat65> sm;
    sm.build(L, M, DM);
    Mat65 m;
    sm.eval(0.2, m);
    CHECK(std::abs(m(4, 2) - std::sin(w * 0.2)) < 1e-6);
}

TEST_CASE("PeriodicHermite::centralDifferences is fourth order")
{
    double err[2];
    for (int j = 0; j < 2; ++j) {
        const int n = j == 0 ? 200 : 400;
        std::vector<double> y(n + 1);
        for (int k = 0; k <= n; ++k)
            y[k] = f(k * L / n);
        const std::vector<double> d = PeriodicHermite<double>::centralDifferences(y, L);
        err[j] = 0.0;
        for (int k = 0; k <= n; ++k)
            err[j] = std::max(err[j], std::abs(d[k] - df(k * L / n)));
    }
    CHECK(err[0] / err[1] > 14.0);
}

TEST_CASE("PeriodicHermite::build rejects bad input")
{
    PeriodicHermite<double> sp;
    CHECK_THROWS_AS(sp.build(0.0, {1.0, 2.0}, {0.0, 0.0}), std::invalid_argument);
    CHECK_THROWS_AS(sp.build(1.0, {1.0, 2.0, 3.0}, {0.0, 0.0}), std::invalid_argument);
}

TEST_CASE("PeriodicHermite::eval: a non-finite s gives NaN, never an out-of-range node")
{
    // Online code calls eval with whatever s it has; NaN / h would be an undefined node index. eval returns NaN
    // values and derivatives instead, which the callers' allFinite checks catch (Tube::frame, Lift::eval, ...).
    std::vector<double> y(11), d(11, 0.0);
    std::vector<Vec7> Y(11, Vec7::Ones()), D(11, Vec7::Zero());
    for (int k = 0; k <= 10; ++k)
        y[k] = std::sin(w * k * L / 10);
    PeriodicHermite<double> sp;
    sp.build(L, y, d);
    PeriodicHermite<Vec7> sv;
    sv.build(L, Y, D);
    const double inf = std::numeric_limits<double>::infinity();
    for (const double s : {std::numeric_limits<double>::quiet_NaN(), inf, -inf}) {
        CAPTURE(s);
        double v = 0.0, dv = 0.0;
        sp.eval(s, v, &dv);
        CHECK(std::isnan(v));
        CHECK(std::isnan(dv));
        Vec7 a = Vec7::Zero(), da = Vec7::Zero();
        sv.eval(s, a, &da);
        CHECK(a.array().isNaN().all());  // every entry, also where the node values are constant
        CHECK(da.array().isNaN().all());
    }
}
