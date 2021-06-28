// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// getXPtrSFBM
SEXP getXPtrSFBM(std::string path, int n, int m, std::vector<size_t> p);
RcppExport SEXP _bigsparser_getXPtrSFBM(SEXP pathSEXP, SEXP nSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< std::vector<size_t> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getXPtrSFBM(path, n, m, p));
    return rcpp_result_gen;
END_RCPP
}
// getXPtrSFBM_compact
SEXP getXPtrSFBM_compact(std::string path, int n, int m, std::vector<size_t> p, const std::vector<int>& first_i);
RcppExport SEXP _bigsparser_getXPtrSFBM_compact(SEXP pathSEXP, SEXP nSEXP, SEXP mSEXP, SEXP pSEXP, SEXP first_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< std::vector<size_t> >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type first_i(first_iSEXP);
    rcpp_result_gen = Rcpp::wrap(getXPtrSFBM_compact(path, n, m, p, first_i));
    return rcpp_result_gen;
END_RCPP
}
// prodVec
NumericVector prodVec(Environment X, const NumericVector& y);
RcppExport SEXP _bigsparser_prodVec(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(prodVec(X, y));
    return rcpp_result_gen;
END_RCPP
}
// cprodVec
NumericVector cprodVec(Environment X, const NumericVector& y);
RcppExport SEXP _bigsparser_cprodVec(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(cprodVec(X, y));
    return rcpp_result_gen;
END_RCPP
}
// sp_solve_sym_eigen
Rcpp::NumericVector sp_solve_sym_eigen(Rcpp::Environment X, const Eigen::VectorXd& b, const Eigen::VectorXd& add_to_diag, double tol, int maxiter);
RcppExport SEXP _bigsparser_sp_solve_sym_eigen(SEXP XSEXP, SEXP bSEXP, SEXP add_to_diagSEXP, SEXP tolSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Environment >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type add_to_diag(add_to_diagSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_solve_sym_eigen(X, b, add_to_diag, tol, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// write_indval
void write_indval(std::string filename, const IntegerVector& i, const NumericVector& x, size_t offset_p, int offset_i);
RcppExport SEXP _bigsparser_write_indval(SEXP filenameSEXP, SEXP iSEXP, SEXP xSEXP, SEXP offset_pSEXP, SEXP offset_iSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type offset_p(offset_pSEXP);
    Rcpp::traits::input_parameter< int >::type offset_i(offset_iSEXP);
    write_indval(filename, i, x, offset_p, offset_i);
    return R_NilValue;
END_RCPP
}
// col_count_sym
IntegerVector col_count_sym(std::vector<size_t> p, const IntegerVector& i);
RcppExport SEXP _bigsparser_col_count_sym(SEXP pSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<size_t> >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(col_count_sym(p, i));
    return rcpp_result_gen;
END_RCPP
}
// write_indval_sym
NumericVector write_indval_sym(std::string filename, std::vector<size_t> p, const IntegerVector& i, const NumericVector& x, size_t offset_p, int offset_i);
RcppExport SEXP _bigsparser_write_indval_sym(SEXP filenameSEXP, SEXP pSEXP, SEXP iSEXP, SEXP xSEXP, SEXP offset_pSEXP, SEXP offset_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::vector<size_t> >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type offset_p(offset_pSEXP);
    Rcpp::traits::input_parameter< int >::type offset_i(offset_iSEXP);
    rcpp_result_gen = Rcpp::wrap(write_indval_sym(filename, p, i, x, offset_p, offset_i));
    return rcpp_result_gen;
END_RCPP
}
// range_col
ListOf<IntegerVector> range_col(const std::vector<size_t>& p, const IntegerVector& i);
RcppExport SEXP _bigsparser_range_col(SEXP pSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(range_col(p, i));
    return rcpp_result_gen;
END_RCPP
}
// range_col_sym
ListOf<IntegerVector> range_col_sym(const std::vector<size_t>& p, const IntegerVector& i);
RcppExport SEXP _bigsparser_range_col_sym(SEXP pSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(range_col_sym(p, i));
    return rcpp_result_gen;
END_RCPP
}
// write_val_compact
List write_val_compact(std::string filename, const std::vector<size_t>& p, const IntegerVector& i, const NumericVector& x, size_t offset_p, int offset_i, bool symmetric);
RcppExport SEXP _bigsparser_write_val_compact(SEXP filenameSEXP, SEXP pSEXP, SEXP iSEXP, SEXP xSEXP, SEXP offset_pSEXP, SEXP offset_iSEXP, SEXP symmetricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type offset_p(offset_pSEXP);
    Rcpp::traits::input_parameter< int >::type offset_i(offset_iSEXP);
    Rcpp::traits::input_parameter< bool >::type symmetric(symmetricSEXP);
    rcpp_result_gen = Rcpp::wrap(write_val_compact(filename, p, i, x, offset_p, offset_i, symmetric));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bigsparser_getXPtrSFBM", (DL_FUNC) &_bigsparser_getXPtrSFBM, 4},
    {"_bigsparser_getXPtrSFBM_compact", (DL_FUNC) &_bigsparser_getXPtrSFBM_compact, 5},
    {"_bigsparser_prodVec", (DL_FUNC) &_bigsparser_prodVec, 2},
    {"_bigsparser_cprodVec", (DL_FUNC) &_bigsparser_cprodVec, 2},
    {"_bigsparser_sp_solve_sym_eigen", (DL_FUNC) &_bigsparser_sp_solve_sym_eigen, 5},
    {"_bigsparser_write_indval", (DL_FUNC) &_bigsparser_write_indval, 5},
    {"_bigsparser_col_count_sym", (DL_FUNC) &_bigsparser_col_count_sym, 2},
    {"_bigsparser_write_indval_sym", (DL_FUNC) &_bigsparser_write_indval_sym, 6},
    {"_bigsparser_range_col", (DL_FUNC) &_bigsparser_range_col, 2},
    {"_bigsparser_range_col_sym", (DL_FUNC) &_bigsparser_range_col_sym, 2},
    {"_bigsparser_write_val_compact", (DL_FUNC) &_bigsparser_write_val_compact, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_bigsparser(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
