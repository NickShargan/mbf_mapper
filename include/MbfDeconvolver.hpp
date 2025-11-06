#pragma once

#include <Eigen/Dense>
#include <vector>

class MbfDeconvolver {
public:
    // ctor: give it AIF and threshold
    MbfDeconvolver(const std::vector<double>& aif,
                   double rel_thresh = 0.01);

    // take tissue curve (std::vector) → return h(t) as Eigen::VectorXd
    Eigen::VectorXd deconvolve(const std::vector<double>& tissue) const;

    // tissue as Eigen vector
    Eigen::VectorXd deconvolve(const Eigen::VectorXd& c) const;

    // convenience: get MBF directly (h(0))
    double computeMBF(const std::vector<double>& tissue) const;

private:
    int T_;
    double rel_thresh_;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_;

    static Eigen::MatrixXd buildConvMatrix(const std::vector<double>& aif);
};
