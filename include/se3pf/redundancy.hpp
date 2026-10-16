#pragma once
#include "se3pf/types.hpp"

namespace se3pf {

    // Redundancy output beta for one redundant joint (n = 7, r = 1)
    class RedundancyOutput {
        public:
            virtual ~RedundancyOutput() = default;

            // computes beta(q) and its differential dbeta
            virtual void evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const = 0;

            // Among configurations with eta(Pi(q)) = eta, the level set {beta = value} is the hyperplane p^T q = b.
            // Used to find start configurations.
            virtual void levelSet(double eta, double value, Vec7& p, double& b) const = 0;
    };

}  // namespace se3pf
