#include "MbfFermiFitter.hpp"
#include "MbfDeconvolver.hpp"   // for MbfDeconvolver::buildConvMatrix
#include <algorithm>
#include <cmath>

MbfFermiFitter::MbfFermiFitter(const Eigen::MatrixXd& A,
                               const Eigen::VectorXd& t,
                               double dt)
    : A_(A), t_(t), dt_(dt)
{
}

MbfFermiFitter::MbfFermiFitter(const std::vector<double>& aif,
                               const std::vector<double>& t,
                               double dt)
    : A_(MbfDeconvolver::buildConvMatrix(aif)),
      t_(Eigen::Map<const Eigen::VectorXd>(t.data(), static_cast<long>(t.size()))),
      dt_(dt)
{
    // A_ *= dt_;
}

Eigen::VectorXd MbfFermiFitter::fit(const Eigen::VectorXd& c,
                                    int max_iters) const
{
    const int T = static_cast<int>(c.size());

    // initial values
    double F   = 1.0;
    double t0  = 1.0;
    double k   = 0.5;
    double tau = 3.0;

    Eigen::VectorXd h(T), model(T), r(T);

    for (int it = 0; it < max_iters; ++it) {
        // h(t) from current params
        fermiImpulse(t_, F, t0, k, tau, h);

        model = A_ * h;

        r = model - c;

        Eigen::MatrixXd J(T, 4);
        const double eps = 1e-4;

        double params[4] = {F, t0, k, tau};
        for (int p = 0; p < 4; ++p) {
            double old = params[p];
            params[p] = old + eps;

            Eigen::VectorXd h_pert(T);
            fermiImpulse(t_, params[0], params[1], params[2], params[3], h_pert);
            Eigen::VectorXd model_pert = A_ * h_pert;

            J.col(p) = (model_pert - model) / eps;

            params[p] = old; // restore
        }

        Eigen::MatrixXd JTJ = J.transpose() * J;
        double lambda = 1e-2;
        JTJ += lambda * Eigen::MatrixXd::Identity(4, 4);

        Eigen::VectorXd JTr   = J.transpose() * r;
        Eigen::VectorXd delta = JTJ.ldlt().solve(-JTr);

        // clamp updates
        for (int i = 0; i < 4; ++i) {
            if (delta(i) > 1.0)  delta(i) = 1.0;
            if (delta(i) < -1.0) delta(i) = -1.0;
        }

        // update params
        F   += delta(0);
        t0  += delta(1);
        k   += delta(2);
        tau += delta(3);

        // keep params in ok ranges
        if (F   < 0.0)  F   = std::max(F, 1e-6);
        if (t0  < 0.0)  t0  = 0.0;
        if (k   < 1e-3) k   = 1e-3;
        if (tau < 1e-3) tau = 1e-3;

        // stopping criterion
        if (delta.norm() < 1e-5)
            break;
    }

    // final h
    fermiImpulse(t_, F, t0, k, tau, h);
    return h;
}

void MbfFermiFitter::fermiImpulse(const Eigen::VectorXd& t,
                                  double F, double t0, double k, double tau,
                                  Eigen::VectorXd& h_out)
{
    const int T = static_cast<int>(t.size());
    h_out.resize(T);
    for (int i = 0; i < T; ++i) {
        double ti = t(i);

        double x        = (ti - t0);
        double logistic = 1.0 / (1.0 + std::exp(-x / (k + 1e-8)));
        double decay    = std::exp(-x / tau);
        double val      = F * logistic * decay;
        if (val < 0.0) val = 0.0;
        h_out(i) = val;
    }
}
