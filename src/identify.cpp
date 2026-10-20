#include "se3pf/identify.hpp"
#include "se3pf/experiments.hpp"  // SafetyStop::torque: the damped hold
#include "se3pf/panda.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace se3pf {

    namespace {
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
    }

    JointIdentification::JointIdentification(const IdentifySettings& settings, const std::vector<int>& joints,
                                             const Vec7& q0)
        : s_(settings), qHold_(q0)
    {
        if (joints.empty())
            throw std::invalid_argument("JointIdentification: no joints");
        for (int j : joints)
            if (j < 0 || j > 6)
                throw std::invalid_argument("JointIdentification: joint index out of range");
        if (!(s_.rampRate > 0.0) || !(s_.stepDuration > 0.0) || !(s_.travelLimit > 0.0) || !(s_.speedLimit > 0.0)
            || !(s_.motionSpeed > 0.0) || !(s_.motionHold >= 0.0) || !(s_.settle > 0.0) || !(s_.blend > 0.0)
            || !(s_.torqueCap.array() > 0.0).all() || !(s_.stepExtra1.array() > 0.0).all()
            || !(s_.stepExtra2.array() > s_.stepExtra1.array()).all()
            || !q0.allFinite())
            throw std::invalid_argument("JointIdentification: invalid settings");
        for (int j : joints)
            for (int d : {1, -1})
                for (int kind : {0, 1, 2})
                    pulses_.push_back({j, d, kind});
        for (auto& b : breakaway_)
            b[0] = b[1] = kNaN;
        lastTau_ = Vec7::Zero();
        beginSettle(0.0, q0);
    }

    double JointIdentification::duration() const noexcept
    {
        double T = 0.0;
        for (const Pulse& p : pulses_) {
            T += s_.settle + s_.blend;
            T += p.kind == 0 ? s_.torqueCap[p.joint] / s_.rampRate : s_.stepDuration;
        }
        return T + s_.settle;
    }

    Vec7 JointIdentification::holdTorque(const Vec7& q, const Vec7& dq) const noexcept
    {
        return SafetyStop::torque(q, dq, qHold_);
    }

    void JointIdentification::beginSettle(double t, const Vec7& q) noexcept
    {
        phase_ = IdentifyPhase::Hold;
        t0_ = t;
        if (joint_ >= 0)
            qHold_[joint_] = q[joint_];  // re-latch the active joint where the pulse left it
        joint_ = -1;
        dir_ = 0;
        pulseKind_ = -1;
        pulseTau_ = 0.0;
        blendFrom_ = lastTau_;
        blendT0_ = t;
    }

    void JointIdentification::beginPulse(double t, const Vec7& q) noexcept
    {
        const Pulse& p = pulses_[next_];
        joint_ = p.joint;
        dir_ = p.dir;
        pulseKind_ = p.kind;
        t0_ = t;
        q0_ = q[joint_];
        movingSince_ = -1.0;
        phase_ = p.kind == 0 ? IdentifyPhase::Ramp : IdentifyPhase::Step;
        const double b = breakaway_[joint_][dir_ > 0 ? 0 : 1];
        stepTau_ = std::isfinite(b) ? b + (p.kind == 2 ? s_.stepExtra2[joint_] : s_.stepExtra1[joint_]) : kNaN;
        blendFrom_ = lastTau_;
        blendT0_ = t;
    }

    Vec7 JointIdentification::step(double t, const Vec7& q, const Vec7& dq, bool& done) noexcept
    {
        done = done_;
        if (done_ || !q.allFinite() || !dq.allFinite() || !std::isfinite(t)) {
            // Finished, or a bad measurement: the hold alone (the driver's monitor reports a NonFinite state)
            Vec7 tau = holdTorque(q.allFinite() ? q : qHold_, dq.allFinite() ? dq : Vec7::Zero());
            lastTau_ = tau;
            return tau;
        }

        // Phase transitions
        if (phase_ == IdentifyPhase::Hold) {
            if (t - t0_ >= s_.settle) {
                if (next_ >= pulses_.size()) {
                    done_ = true;
                    done = true;
                } else {
                    beginPulse(t, q);
                    // A step without a breakaway (the ramp reached its cap) is skipped
                    if (phase_ == IdentifyPhase::Step && !std::isfinite(stepTau_)) {
                        ++next_;
                        beginSettle(t, q);
                    }
                }
            }
        }
        Vec7 tau = holdTorque(q, dq);
        if (phase_ != IdentifyPhase::Hold) {
            const int j = joint_;
            const double travel = std::abs(q[j] - q0_);
            const double speed = std::abs(dq[j]);
            const bool moving = speed >= s_.motionSpeed;
            if (!moving)
                movingSince_ = -1.0;
            else if (movingSince_ < 0.0)
                movingSince_ = t;
            const bool brokeAway = movingSince_ >= 0.0 && t - movingSince_ >= s_.motionHold;
            bool over = travel >= s_.travelLimit || speed >= s_.speedLimit;
            if (phase_ == IdentifyPhase::Ramp) {
                pulseTau_ = dir_ * s_.rampRate * (t - t0_);
                if (brokeAway) {
                    // The torque when the motion began (the ramp's value motionHold ago)
                    breakaway_[j][dir_ > 0 ? 0 : 1] = s_.rampRate * (movingSince_ - t0_);
                    over = true;
                } else if (std::abs(pulseTau_) >= s_.torqueCap[j]) {
                    pulseTau_ = dir_ * s_.torqueCap[j];
                    over = true;  // not found: stays NaN
                }
            } else {
                pulseTau_ = dir_ * stepTau_;
                if (t - t0_ >= s_.stepDuration)
                    over = true;
            }
            tau[j] = pulseTau_;
            if (over) {
                ++next_;
                lastTau_ = tau;
                beginSettle(t, q);
                tau = holdTorque(q, dq);
            }
        }
        // Blend from the previous phase's last command over s_.blend (the torque-rate margin)
        if (blendT0_ >= 0.0) {
            const double w = std::clamp((t - blendT0_) / s_.blend, 0.0, 1.0);
            tau = blendFrom_ + w * (tau - blendFrom_);
            if (w >= 1.0)
                blendT0_ = -1.0;
        }
        lastTau_ = tau;
        return tau;
    }

    // ---------------------------------------------------------------------------------------------------------
    // The fit

    std::vector<JointIdentificationResult> fitIdentification(const std::vector<IdentifyLogRow>& rows,
                                                             const IdentifySettings& s)
    {
        std::vector<JointIdentificationResult> out(7);
        // Each pulse: the rows with the same (joint, direction, pulse kind), contiguous in time
        struct Window {
            int joint = -1, dir = 0, kind = -1;
            std::size_t first = 0, last = 0;
        };
        std::vector<Window> windows;
        for (std::size_t k = 0; k < rows.size(); ++k) {
            const IdentifyLogRow& r = rows[k];
            if (r.phase != 1 && r.phase != 2)
                continue;
            if (windows.empty() || windows.back().joint != r.joint - 1 || windows.back().dir != r.direction
                || windows.back().kind != r.pulse || windows.back().last + 1 != k)
                windows.push_back({r.joint - 1, r.direction, r.pulse, k, k});
            windows.back().last = k;
        }
        auto onset = [&](const Window& w) -> std::size_t {  // the first row of sustained motion, or npos
            std::size_t since = std::string::npos;
            for (std::size_t k = w.first; k <= w.last; ++k) {
                const int j = w.joint;
                if (std::abs(rows[k].dq[j]) >= s.motionSpeed) {
                    if (since == std::string::npos)
                        since = k;
                    // The stepper ends the ramp in the cycle that completes motionHold, and that cycle is logged
                    // as the settle, so the ramp's rows cover motionHold less one period
                    if (rows[k].t - rows[since].t >= s.motionHold - 1.5e-3)
                        return since;
                } else {
                    since = std::string::npos;
                }
            }
            return std::string::npos;
        };
        struct StepFit {
            bool ok = false;
            double tau = 0.0, tauJ = 0.0, a = 0.0;
        };
        StepFit steps[7][2][2];
        for (const Window& w : windows) {
            if (w.joint < 0 || w.joint > 6)
                continue;
            const int d = w.dir > 0 ? 0 : 1;
            JointIdentificationResult& r = out[static_cast<std::size_t>(w.joint)];
            const std::size_t k0 = onset(w);
            if (w.kind == 0) {
                r.attempted[d] = true;
                if (k0 == std::string::npos)
                    continue;
                r.found[d] = true;
                r.breakawayRaw[d] = std::abs(rows[k0].pulseTorque);
                r.breakawayCmd[d] = r.breakawayRaw[d];
                // The sensor's torque beyond gravity just before the motion: the mean over the 20 ms before onset
                double sum = 0.0;
                int n = 0;
                for (std::size_t k = k0; k > w.first && rows[k0].t - rows[k].t <= 0.02; --k) {
                    sum += rows[k].tauJ[w.joint] - rows[k].G[w.joint];
                    ++n;
                }
                r.breakawaySensor[d] = n > 0 ? std::abs(sum / n) : 0.0;
            } else if (w.kind == 1 || w.kind == 2) {
                if (k0 == std::string::npos || w.last <= k0 + 10)
                    continue;
                // Mean acceleration from the onset to the end of the pulse, mean command and sensor torque over it
                const double dt = rows[w.last].t - rows[k0].t;
                if (!(dt > 0.02))
                    continue;
                StepFit& f = steps[w.joint][d][w.kind - 1];
                f.a = (rows[w.last].dq[w.joint] - rows[k0].dq[w.joint]) / dt;
                double st = 0.0, sj = 0.0;
                for (std::size_t k = k0; k <= w.last; ++k) {
                    st += rows[k].pulseTorque;
                    sj += rows[k].tauJ[w.joint] - rows[k].G[w.joint];
                }
                f.tau = st / static_cast<double>(w.last - k0 + 1);
                f.tauJ = sj / static_cast<double>(w.last - k0 + 1);
                f.ok = true;
                r.accel[d][w.kind - 1] = f.a;
            }
        }
        for (int j = 0; j < 7; ++j)
            for (int d = 0; d < 2; ++d) {
                const StepFit& a = steps[j][d][0];
                const StepFit& b = steps[j][d][1];
                if (!a.ok || !b.ok)
                    continue;
                const double da = b.a - a.a;
                if (!(std::abs(da) > 0.05))  // the two steps must accelerate differently
                    continue;
                JointIdentificationResult& r = out[static_cast<std::size_t>(j)];
                r.inertiaFound[d] = true;
                r.inertiaCmd[d] = (b.tau - a.tau) / da;
                r.inertiaSensor[d] = (b.tauJ - a.tauJ) / da;
                r.kineticFriction[d] = std::abs(a.tau - r.inertiaCmd[d] * a.a);
                if (r.found[d] && r.inertiaCmd[d] > 0.0)
                    r.breakawayCmd[d] = std::max(
                        0.0, r.breakawayRaw[d] - std::sqrt(2.0 * r.inertiaCmd[d] * s.motionSpeed * s.rampRate));
            }
        return out;
    }

    std::string identificationSummary(const std::vector<JointIdentificationResult>& r)
    {
        std::string s = "joint dir | breakaway cmd (sensor) N m | kinetic friction N m | inertia cmd (sensor) kg m^2 | "
                        "accelerations of the two steps rad/s^2   (breakaway corrected for the ramp's lag where the "
                        "inertia is known)\n";
        char buf[256];
        for (std::size_t j = 0; j < r.size(); ++j)
            for (int d = 0; d < 2; ++d) {
                const JointIdentificationResult& x = r[j];
                if (!x.attempted[d] && !x.found[d] && !x.inertiaFound[d])
                    continue;
                std::snprintf(buf, sizeof(buf), "  %zu   %c  | ", j + 1, d == 0 ? '+' : '-');
                s += buf;
                if (x.found[d])
                    std::snprintf(buf, sizeof(buf), "%.3f (%.3f)            | ", x.breakawayCmd[d], x.breakawaySensor[d]);
                else
                    std::snprintf(buf, sizeof(buf), "not found (ramp to cap)   | ");
                s += buf;
                if (x.inertiaFound[d])
                    std::snprintf(buf, sizeof(buf), "%.3f                | %.4f (%.4f)              | %.2f, %.2f\n",
                                  x.kineticFriction[d], x.inertiaCmd[d], x.inertiaSensor[d], x.accel[d][0],
                                  x.accel[d][1]);
                else
                    std::snprintf(buf, sizeof(buf), "-                    | not identified              | %.2f, %.2f\n",
                                  x.accel[d][0], x.accel[d][1]);
                s += buf;
            }
        if (s.find("\n  ") == std::string::npos)
            s += "  (nothing identified)\n";
        return s;
    }

}  // namespace se3pf
