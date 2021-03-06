// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// hybridmonitor_single
Rcpp::List hybridmonitor_single(const arma::mat& y, const arma::vec& y_Kp1, const double& sig_level, const double& z_alpha, const int& B1, const int& B2, const int& B3, const bool verbose, const int& ncores);
RcppExport SEXP _energystuff_hybridmonitor_single(SEXP ySEXP, SEXP y_Kp1SEXP, SEXP sig_levelSEXP, SEXP z_alphaSEXP, SEXP B1SEXP, SEXP B2SEXP, SEXP B3SEXP, SEXP verboseSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_Kp1(y_Kp1SEXP);
    Rcpp::traits::input_parameter< const double& >::type sig_level(sig_levelSEXP);
    Rcpp::traits::input_parameter< const double& >::type z_alpha(z_alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int& >::type B2(B2SEXP);
    Rcpp::traits::input_parameter< const int& >::type B3(B3SEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(hybridmonitor_single(y, y_Kp1, sig_level, z_alpha, B1, B2, B3, verbose, ncores));
    return rcpp_result_gen;
END_RCPP
}
// hybridmonitor_multiple
Rcpp::List hybridmonitor_multiple(const arma::cube& y, const arma::mat& y_Kp1, const double& sig_level, const double& z_alpha, const int& B1, const int& B2, const int& B3, const int& T, const bool verbose, const int& ncores);
RcppExport SEXP _energystuff_hybridmonitor_multiple(SEXP ySEXP, SEXP y_Kp1SEXP, SEXP sig_levelSEXP, SEXP z_alphaSEXP, SEXP B1SEXP, SEXP B2SEXP, SEXP B3SEXP, SEXP TSEXP, SEXP verboseSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y_Kp1(y_Kp1SEXP);
    Rcpp::traits::input_parameter< const double& >::type sig_level(sig_levelSEXP);
    Rcpp::traits::input_parameter< const double& >::type z_alpha(z_alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int& >::type B2(B2SEXP);
    Rcpp::traits::input_parameter< const int& >::type B3(B3SEXP);
    Rcpp::traits::input_parameter< const int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(hybridmonitor_multiple(y, y_Kp1, sig_level, z_alpha, B1, B2, B3, T, verbose, ncores));
    return rcpp_result_gen;
END_RCPP
}
// hybridmonitor_multiple_select_months
Rcpp::List hybridmonitor_multiple_select_months(const arma::cube& y, const arma::mat& y_Kp1_, const double& sig_level, const double& z_alpha_K, const int& B1, const int& B2, const int& B3, const int& T, const arma::uvec& m_idx, const bool do_bonferroni, const bool do_holm, const bool do_hochberg, const bool do_hommel, const bool do_BH, const bool do_BY, const bool verbose, const int& ncores);
RcppExport SEXP _energystuff_hybridmonitor_multiple_select_months(SEXP ySEXP, SEXP y_Kp1_SEXP, SEXP sig_levelSEXP, SEXP z_alpha_KSEXP, SEXP B1SEXP, SEXP B2SEXP, SEXP B3SEXP, SEXP TSEXP, SEXP m_idxSEXP, SEXP do_bonferroniSEXP, SEXP do_holmSEXP, SEXP do_hochbergSEXP, SEXP do_hommelSEXP, SEXP do_BHSEXP, SEXP do_BYSEXP, SEXP verboseSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y_Kp1_(y_Kp1_SEXP);
    Rcpp::traits::input_parameter< const double& >::type sig_level(sig_levelSEXP);
    Rcpp::traits::input_parameter< const double& >::type z_alpha_K(z_alpha_KSEXP);
    Rcpp::traits::input_parameter< const int& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int& >::type B2(B2SEXP);
    Rcpp::traits::input_parameter< const int& >::type B3(B3SEXP);
    Rcpp::traits::input_parameter< const int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type m_idx(m_idxSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_bonferroni(do_bonferroniSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_holm(do_holmSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_hochberg(do_hochbergSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_hommel(do_hommelSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_BH(do_BHSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_BY(do_BYSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(hybridmonitor_multiple_select_months(y, y_Kp1_, sig_level, z_alpha_K, B1, B2, B3, T, m_idx, do_bonferroni, do_holm, do_hochberg, do_hommel, do_BH, do_BY, verbose, ncores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_energystuff_hybridmonitor_single", (DL_FUNC) &_energystuff_hybridmonitor_single, 9},
    {"_energystuff_hybridmonitor_multiple", (DL_FUNC) &_energystuff_hybridmonitor_multiple, 10},
    {"_energystuff_hybridmonitor_multiple_select_months", (DL_FUNC) &_energystuff_hybridmonitor_multiple_select_months, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_energystuff(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
