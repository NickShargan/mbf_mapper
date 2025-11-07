#pragma once

#include <Eigen/Dense>
#include <vector>

class MbfFermiFitter {
public:
    // A is convolution matrix (T x T)
    // t is time vector (T)
    MbfFermiFitter(const Eigen::MatrixXd& A,
                   const Eigen::VectorXd& t,
                   double dt);

    // ctor from std::vector inputs
    MbfFermiFitter(const std::vector<double>& aif,
                   const std::vector<double>& t,
                   double dt);

    // fit tissue curve c → return h(t)
    Eigen::VectorXd fit(const Eigen::VectorXd& c,
                        int max_iters = 30) const;

private:
    Eigen::MatrixXd A_;
    Eigen::VectorXd t_;
    double dt_;

    // internal Fermi impulse builder
    static void fermiImpulse(const Eigen::VectorXd& t,
                             double F, double t0, double k, double tau,
                             Eigen::VectorXd& h_out);
};
