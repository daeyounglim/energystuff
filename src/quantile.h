#ifndef TSSTUFF_QUANTILE_CPP_H
#define TSSTUFF_QUANTILE_CPP_H
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo)]]

template <typename T>
double quantile_cpp(const T& x, const double& prob) {
    const int n = x.n_elem;
    double idx = 1.0 + static_cast<double>(n-1) * prob;
    int lo = static_cast<int>(std::floor(idx));
    T x_sorted = arma::sort(x);
    double qs = x_sorted(lo-1);
    double g = idx - static_cast<double>(lo);
    return (1.0 - g) * qs + g * x_sorted(lo);
}

arma::mat quantile_m_cpp(const arma::cube& x, const double& prob) {
    using namespace arma;
    int nr = x.n_rows;
    int nc = x.n_cols;
    mat out(nr, nc, fill::zeros);
    for (int ir = 0; ir < nr; ++ir) {
        for (int ic = 0; ic < nc; ++ic) {
            rowvec x_rc = x.tube(ir, ic);
            out(ir,ic) = quantile_cpp(x_rc, prob);
        }
    }
    return out;
}

arma::vec quantile_v_cpp(const arma::mat& x, const double& prob) {
    using namespace arma;
    int nr = x.n_rows;
    vec out(nr, fill::zeros);
    for (int ir = 0; ir < nr; ++ir) {
        rowvec x_r = x.row(ir);
        out(ir) = quantile_cpp(x_r, prob);
    }
    return out;
}

#endif