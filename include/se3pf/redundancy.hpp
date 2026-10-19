#pragma once
#include "se3pf/types.hpp"

namespace se3pf {

    // Redundancy output beta for one redundant joint (n = 7, r = 1), eq:beta-flat and eq:dbeta. An interface, so
    // the controller can take any output: beta = X(eta)^T (q - q*(eta)) of Lemma 3 (in lift), or a joint coupling.
    // Implementations must not allocate or throw: the controller calls evaluate in every 1 kHz step.
    class RedundancyOutput {
        public:
            virtual ~RedundancyOutput() = default;

            // beta(q) and its differential dbeta (a row, stored as a column), given eta(q) and d eta(q), the first
            // row of J_Psi^-1 J_b, which the controller has already computed.
            virtual void evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const = 0;

            // Among configurations with eta(Pi(q)) = eta, the level set {beta = value} is the hyperplane p^T q = b.
            // Used to find start configurations.
            virtual void levelSet(double eta, double value, Vec7& p, double& b) const = 0;
    };

}  // namespace se3pf
