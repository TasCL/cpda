// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logLik_norm
double logLik_norm(arma::vec object, arma::vec pVec, int ns, double h, double m, double p);
RcppExport SEXP cpda_logLik_norm(SEXP objectSEXP, SEXP pVecSEXP, SEXP nsSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type object(objectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_norm(object, pVec, ns, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// logLik_norm2
Rcpp::List logLik_norm2(arma::vec object, arma::vec pVec, int ns, double h, double m, double p);
RcppExport SEXP cpda_logLik_norm2(SEXP objectSEXP, SEXP pVecSEXP, SEXP nsSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type object(objectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_norm2(object, pVec, ns, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// logLik_plba
double logLik_plba(arma::mat object, arma::vec pVec, int ns, double h, double m, double p);
RcppExport SEXP cpda_logLik_plba(SEXP objectSEXP, SEXP pVecSEXP, SEXP nsSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type object(objectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_plba(object, pVec, ns, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// logLik_plba2
Rcpp::List logLik_plba2(arma::mat object, arma::vec pVec, int ns, double h, double m, double p);
RcppExport SEXP cpda_logLik_plba2(SEXP objectSEXP, SEXP pVecSEXP, SEXP nsSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type object(objectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_plba2(object, pVec, ns, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// logLik_lba
double logLik_lba(arma::mat object, arma::vec pVec, int ns, double h, double m, double p);
RcppExport SEXP cpda_logLik_lba(SEXP objectSEXP, SEXP pVecSEXP, SEXP nsSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type object(objectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_lba(object, pVec, ns, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// logLik_lba2
Rcpp::List logLik_lba2(arma::mat object, arma::vec pVec, int ns, double h, double m, double p);
RcppExport SEXP cpda_logLik_lba2(SEXP objectSEXP, SEXP pVecSEXP, SEXP nsSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type object(objectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_lba2(object, pVec, ns, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// logLik_fft
double logLik_fft(arma::vec y, arma::vec yhat, double h, double m, double p);
RcppExport SEXP cpda_logLik_fft(SEXP ySEXP, SEXP yhatSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yhat(yhatSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_fft(y, yhat, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// logLik_fft2
Rcpp::List logLik_fft2(arma::vec y, arma::vec yhat, double h, double m, double p);
RcppExport SEXP cpda_logLik_fft2(SEXP ySEXP, SEXP yhatSEXP, SEXP hSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yhat(yhatSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_fft2(y, yhat, h, m, p));
    return rcpp_result_gen;
END_RCPP
}
// rlba
arma::mat rlba(int n, arma::vec pVec);
RcppExport SEXP cpda_rlba(SEXP nSEXP, SEXP pVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    rcpp_result_gen = Rcpp::wrap(rlba(n, pVec));
    return rcpp_result_gen;
END_RCPP
}
// rplba
arma::mat rplba(int n, arma::vec pVec);
RcppExport SEXP cpda_rplba(SEXP nSEXP, SEXP pVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    rcpp_result_gen = Rcpp::wrap(rplba(n, pVec));
    return rcpp_result_gen;
END_RCPP
}
// rplba_omp
arma::mat rplba_omp(int n, arma::vec pVec);
RcppExport SEXP cpda_rplba_omp(SEXP nSEXP, SEXP pVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    rcpp_result_gen = Rcpp::wrap(rplba_omp(n, pVec));
    return rcpp_result_gen;
END_RCPP
}
// choiceDT
arma::mat choiceDT(arma::mat data, arma::vec pVec);
RcppExport SEXP cpda_choiceDT(SEXP dataSEXP, SEXP pVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pVec(pVecSEXP);
    rcpp_result_gen = Rcpp::wrap(choiceDT(data, pVec));
    return rcpp_result_gen;
END_RCPP
}
// cquantile
double cquantile(arma::vec y, double q);
RcppExport SEXP cpda_cquantile(SEXP ySEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(cquantile(y, q));
    return rcpp_result_gen;
END_RCPP
}
// bwNRD0
double bwNRD0(arma::vec y, double m);
RcppExport SEXP cpda_bwNRD0(SEXP ySEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(bwNRD0(y, m));
    return rcpp_result_gen;
END_RCPP
}
