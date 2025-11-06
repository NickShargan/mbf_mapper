#include "MbfDeconvolver.hpp"

#include <Eigen/SVD>

// ctor
MbfDeconvolver::MbfDeconvolver(const std::vector<double>& aif,
                               double rel_thresh)
    : T_(static_cast<int>(aif.size())),
      rel_thresh_(rel_thresh)
{
    // build convolution matrix
    Eigen::MatrixXd A = buildConvMatrix(aif);

    // SVD once
    svd_ = Eigen::JacobiSVD<Eigen::MatrixXd>(
        A, Eigen::ComputeThinU | Eigen::ComputeThinV
    );
}

Eigen::VectorXd MbfDeconvolver::deconvolve(const std::vector<double>& tissue) const
{
    Eigen::VectorXd c(T_);
    for (int i = 0; i < T_; ++i)
        c(i) = tissue[i];
    return deconvolve(c);
}

Eigen::VectorXd MbfDeconvolver::deconvolve(const Eigen::VectorXd& c) const
{
    const auto& U = svd_.matrixU();
    const auto& V = svd_.matrixV();
    const auto& S = svd_.singularValues();

    Eigen::VectorXd Sinv = Eigen::VectorXd::Zero(T_);
    double smax = S(0);
    for (int i = 0; i < T_; ++i) {
        if (S(i) / smax > rel_thresh_) {
            Sinv(i) = 1.0 / S(i);
        }
    }

    Eigen::VectorXd Ut_c = U.transpose() * c;
    Eigen::VectorXd SigmaInv_Ut_c = Sinv.asDiagonal() * Ut_c;
    Eigen::VectorXd h = V * SigmaInv_Ut_c;

    // zero clamp
    for (int i = 0; i < T_; ++i) {
        if (h(i) < 0.0)
            h(i) = 0.0;
    }

    return h;
}

double MbfDeconvolver::computeMBF(const std::vector<double>& tissue) const
{
    Eigen::VectorXd h = deconvolve(tissue);
    return h(0);
}

Eigen::MatrixXd MbfDeconvolver::buildConvMatrix(const std::vector<double>& aif)
{
    int T = static_cast<int>(aif.size());
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(T, T);
    for (int i = 0; i < T; ++i) {
        for (int j = 0; j <= i; ++j) {
            A(i, j) = aif[i - j];
        }
    }
    return A;
}
