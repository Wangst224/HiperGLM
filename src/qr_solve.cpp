#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::VectorXd qr_solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::VectorXd x = qr.solve(b);

    return x;
}
