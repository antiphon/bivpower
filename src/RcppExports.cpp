// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// doublesum
NumericVector doublesum(NumericMatrix xy, NumericVector mark, NumericVector r);
RcppExport SEXP _bivpower_doublesum(SEXP xySEXP, SEXP markSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mark(markSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(doublesum(xy, mark, r));
    return rcpp_result_gen;
END_RCPP
}
// doublesum22
NumericVector doublesum22(NumericMatrix xy, NumericVector r);
RcppExport SEXP _bivpower_doublesum22(SEXP xySEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(doublesum22(xy, r));
    return rcpp_result_gen;
END_RCPP
}
// g_doublesum
NumericVector g_doublesum(NumericMatrix xy, NumericVector mark, NumericVector r, double bw, int kern);
RcppExport SEXP _bivpower_g_doublesum(SEXP xySEXP, SEXP markSEXP, SEXP rSEXP, SEXP bwSEXP, SEXP kernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mark(markSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type kern(kernSEXP);
    rcpp_result_gen = Rcpp::wrap(g_doublesum(xy, mark, r, bw, kern));
    return rcpp_result_gen;
END_RCPP
}
// g_doublesum_asym
NumericMatrix g_doublesum_asym(NumericMatrix xy, NumericVector mark, NumericVector r, double bw12, double bw21, int kern);
RcppExport SEXP _bivpower_g_doublesum_asym(SEXP xySEXP, SEXP markSEXP, SEXP rSEXP, SEXP bw12SEXP, SEXP bw21SEXP, SEXP kernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mark(markSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type bw12(bw12SEXP);
    Rcpp::traits::input_parameter< double >::type bw21(bw21SEXP);
    Rcpp::traits::input_parameter< int >::type kern(kernSEXP);
    rcpp_result_gen = Rcpp::wrap(g_doublesum_asym(xy, mark, r, bw12, bw21, kern));
    return rcpp_result_gen;
END_RCPP
}
// evaluate_lambda_c
NumericVector evaluate_lambda_c(List snfield, NumericMatrix x);
RcppExport SEXP _bivpower_evaluate_lambda_c(SEXP snfieldSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type snfield(snfieldSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(evaluate_lambda_c(snfield, x));
    return rcpp_result_gen;
END_RCPP
}
// rcox_MH
List rcox_MH(int n, NumericVector win, List snfield, int iter, int dbg);
RcppExport SEXP _bivpower_rcox_MH(SEXP nSEXP, SEXP winSEXP, SEXP snfieldSEXP, SEXP iterSEXP, SEXP dbgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type win(winSEXP);
    Rcpp::traits::input_parameter< List >::type snfield(snfieldSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type dbg(dbgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcox_MH(n, win, snfield, iter, dbg));
    return rcpp_result_gen;
END_RCPP
}
// rcox_thin
List rcox_thin(NumericVector win, List snfield, int dbg);
RcppExport SEXP _bivpower_rcox_thin(SEXP winSEXP, SEXP snfieldSEXP, SEXP dbgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type win(winSEXP);
    Rcpp::traits::input_parameter< List >::type snfield(snfieldSEXP);
    Rcpp::traits::input_parameter< int >::type dbg(dbgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcox_thin(win, snfield, dbg));
    return rcpp_result_gen;
END_RCPP
}
// mc_triple_integral_ball
NumericVector mc_triple_integral_ball(NumericMatrix bbox, NumericVector r, int n);
RcppExport SEXP _bivpower_mc_triple_integral_ball(SEXP bboxSEXP, SEXP rSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mc_triple_integral_ball(bbox, r, n));
    return rcpp_result_gen;
END_RCPP
}
// mc_triple_integral_ball_matrix
NumericMatrix mc_triple_integral_ball_matrix(NumericMatrix bbox, NumericVector r, int n);
RcppExport SEXP _bivpower_mc_triple_integral_ball_matrix(SEXP bboxSEXP, SEXP rSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mc_triple_integral_ball_matrix(bbox, r, n));
    return rcpp_result_gen;
END_RCPP
}
// mc_triple_integral_annuli
NumericVector mc_triple_integral_annuli(NumericMatrix bbox, NumericVector r, double h, int n);
RcppExport SEXP _bivpower_mc_triple_integral_annuli(SEXP bboxSEXP, SEXP rSEXP, SEXP hSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mc_triple_integral_annuli(bbox, r, h, n));
    return rcpp_result_gen;
END_RCPP
}
// mc_triple_integral_annuli_matrix
NumericMatrix mc_triple_integral_annuli_matrix(NumericMatrix bbox, NumericVector r, NumericVector h, int n);
RcppExport SEXP _bivpower_mc_triple_integral_annuli_matrix(SEXP bboxSEXP, SEXP rSEXP, SEXP hSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mc_triple_integral_annuli_matrix(bbox, r, h, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bivpower_doublesum", (DL_FUNC) &_bivpower_doublesum, 3},
    {"_bivpower_doublesum22", (DL_FUNC) &_bivpower_doublesum22, 2},
    {"_bivpower_g_doublesum", (DL_FUNC) &_bivpower_g_doublesum, 5},
    {"_bivpower_g_doublesum_asym", (DL_FUNC) &_bivpower_g_doublesum_asym, 6},
    {"_bivpower_evaluate_lambda_c", (DL_FUNC) &_bivpower_evaluate_lambda_c, 2},
    {"_bivpower_rcox_MH", (DL_FUNC) &_bivpower_rcox_MH, 5},
    {"_bivpower_rcox_thin", (DL_FUNC) &_bivpower_rcox_thin, 3},
    {"_bivpower_mc_triple_integral_ball", (DL_FUNC) &_bivpower_mc_triple_integral_ball, 3},
    {"_bivpower_mc_triple_integral_ball_matrix", (DL_FUNC) &_bivpower_mc_triple_integral_ball_matrix, 3},
    {"_bivpower_mc_triple_integral_annuli", (DL_FUNC) &_bivpower_mc_triple_integral_annuli, 4},
    {"_bivpower_mc_triple_integral_annuli_matrix", (DL_FUNC) &_bivpower_mc_triple_integral_annuli_matrix, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bivpower(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
