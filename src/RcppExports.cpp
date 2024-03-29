// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/RobustCalibration.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Chol_Eigen
MatrixXd Chol_Eigen(const Eigen::MatrixXd R);
RcppExport SEXP _RobustCalibration_Chol_Eigen(SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(Chol_Eigen(R));
    return rcpp_result_gen;
END_RCPP
}
// Mogihammer
Eigen::VectorXd Mogihammer(const MatrixXd obsCoords, const VectorXd m, int simul_type);
RcppExport SEXP _RobustCalibration_Mogihammer(SEXP obsCoordsSEXP, SEXP mSEXP, SEXP simul_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd >::type obsCoords(obsCoordsSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type simul_type(simul_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Mogihammer(obsCoords, m, simul_type));
    return rcpp_result_gen;
END_RCPP
}
// Accept_proposal
bool Accept_proposal(double r);
RcppExport SEXP _RobustCalibration_Accept_proposal(SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(Accept_proposal(r));
    return rcpp_result_gen;
END_RCPP
}
// Get_R_z_new
MatrixXd Get_R_z_new(const Eigen::VectorXd beta_delta, const double eta_delta, const double lambda_z, const List R0, const String kernel_type, const Eigen::VectorXd alpha, const Eigen::VectorXd inv_output_weights);
RcppExport SEXP _RobustCalibration_Get_R_z_new(SEXP beta_deltaSEXP, SEXP eta_deltaSEXP, SEXP lambda_zSEXP, SEXP R0SEXP, SEXP kernel_typeSEXP, SEXP alphaSEXP, SEXP inv_output_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type beta_delta(beta_deltaSEXP);
    Rcpp::traits::input_parameter< const double >::type eta_delta(eta_deltaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_z(lambda_zSEXP);
    Rcpp::traits::input_parameter< const List >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type inv_output_weights(inv_output_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_R_z_new(beta_delta, eta_delta, lambda_z, R0, kernel_type, alpha, inv_output_weights));
    return rcpp_result_gen;
END_RCPP
}
// Get_R_new
MatrixXd Get_R_new(const Eigen::VectorXd beta_delta, const double eta_delta, const List R0, const String kernel_type, const Eigen::VectorXd alpha, const Eigen::VectorXd inv_output_weights);
RcppExport SEXP _RobustCalibration_Get_R_new(SEXP beta_deltaSEXP, SEXP eta_deltaSEXP, SEXP R0SEXP, SEXP kernel_typeSEXP, SEXP alphaSEXP, SEXP inv_output_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type beta_delta(beta_deltaSEXP);
    Rcpp::traits::input_parameter< const double >::type eta_delta(eta_deltaSEXP);
    Rcpp::traits::input_parameter< const List >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type inv_output_weights(inv_output_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_R_new(beta_delta, eta_delta, R0, kernel_type, alpha, inv_output_weights));
    return rcpp_result_gen;
END_RCPP
}
// Sample_sigma_2_theta_m
Eigen::VectorXd Sample_sigma_2_theta_m(const Eigen::VectorXd param, const Eigen::MatrixXd L_cur, const Eigen::VectorXd output, const int p_theta, const int p_x, Eigen::MatrixXd X, bool have_mean, const VectorXd cm_obs, const double S_2_f, const int num_obs_all);
RcppExport SEXP _RobustCalibration_Sample_sigma_2_theta_m(SEXP paramSEXP, SEXP L_curSEXP, SEXP outputSEXP, SEXP p_thetaSEXP, SEXP p_xSEXP, SEXP XSEXP, SEXP have_meanSEXP, SEXP cm_obsSEXP, SEXP S_2_fSEXP, SEXP num_obs_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type L_cur(L_curSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const int >::type p_theta(p_thetaSEXP);
    Rcpp::traits::input_parameter< const int >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type have_mean(have_meanSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type cm_obs(cm_obsSEXP);
    Rcpp::traits::input_parameter< const double >::type S_2_f(S_2_fSEXP);
    Rcpp::traits::input_parameter< const int >::type num_obs_all(num_obs_allSEXP);
    rcpp_result_gen = Rcpp::wrap(Sample_sigma_2_theta_m(param, L_cur, output, p_theta, p_x, X, have_mean, cm_obs, S_2_f, num_obs_all));
    return rcpp_result_gen;
END_RCPP
}
// Log_marginal_post
double Log_marginal_post(const Eigen::VectorXd param, Eigen::MatrixXd L_cur, const Eigen::VectorXd output, const int p_theta, int p_x, Eigen::MatrixXd X, bool have_mean, const VectorXd CL, const double a, const double b, const VectorXd cm_obs, const double S_2_f, const int num_obs_all);
RcppExport SEXP _RobustCalibration_Log_marginal_post(SEXP paramSEXP, SEXP L_curSEXP, SEXP outputSEXP, SEXP p_thetaSEXP, SEXP p_xSEXP, SEXP XSEXP, SEXP have_meanSEXP, SEXP CLSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cm_obsSEXP, SEXP S_2_fSEXP, SEXP num_obs_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type L_cur(L_curSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const int >::type p_theta(p_thetaSEXP);
    Rcpp::traits::input_parameter< int >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type have_mean(have_meanSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type CL(CLSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type cm_obs(cm_obsSEXP);
    Rcpp::traits::input_parameter< const double >::type S_2_f(S_2_fSEXP);
    Rcpp::traits::input_parameter< const int >::type num_obs_all(num_obs_allSEXP);
    rcpp_result_gen = Rcpp::wrap(Log_marginal_post(param, L_cur, output, p_theta, p_x, X, have_mean, CL, a, b, cm_obs, S_2_f, num_obs_all));
    return rcpp_result_gen;
END_RCPP
}
// Sample_sigma_2_theta_m_no_discrepancy
Eigen::VectorXd Sample_sigma_2_theta_m_no_discrepancy(const Eigen::VectorXd param, const Eigen::VectorXd output, const int p_theta, Eigen::MatrixXd X, bool have_mean, VectorXd inv_output_weights, const VectorXd cm_obs, const double S_2_f, const int num_obs_all);
RcppExport SEXP _RobustCalibration_Sample_sigma_2_theta_m_no_discrepancy(SEXP paramSEXP, SEXP outputSEXP, SEXP p_thetaSEXP, SEXP XSEXP, SEXP have_meanSEXP, SEXP inv_output_weightsSEXP, SEXP cm_obsSEXP, SEXP S_2_fSEXP, SEXP num_obs_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const int >::type p_theta(p_thetaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type have_mean(have_meanSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type inv_output_weights(inv_output_weightsSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type cm_obs(cm_obsSEXP);
    Rcpp::traits::input_parameter< const double >::type S_2_f(S_2_fSEXP);
    Rcpp::traits::input_parameter< const int >::type num_obs_all(num_obs_allSEXP);
    rcpp_result_gen = Rcpp::wrap(Sample_sigma_2_theta_m_no_discrepancy(param, output, p_theta, X, have_mean, inv_output_weights, cm_obs, S_2_f, num_obs_all));
    return rcpp_result_gen;
END_RCPP
}
// Log_marginal_post_no_discrepancy
double Log_marginal_post_no_discrepancy(const Eigen::VectorXd param, const Eigen::VectorXd output, const int p_theta, Eigen::MatrixXd X, bool have_mean, VectorXd inv_output_weights, const VectorXd cm_obs, const double S_2_f, const int num_obs_all);
RcppExport SEXP _RobustCalibration_Log_marginal_post_no_discrepancy(SEXP paramSEXP, SEXP outputSEXP, SEXP p_thetaSEXP, SEXP XSEXP, SEXP have_meanSEXP, SEXP inv_output_weightsSEXP, SEXP cm_obsSEXP, SEXP S_2_fSEXP, SEXP num_obs_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const int >::type p_theta(p_thetaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type have_mean(have_meanSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type inv_output_weights(inv_output_weightsSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type cm_obs(cm_obsSEXP);
    Rcpp::traits::input_parameter< const double >::type S_2_f(S_2_fSEXP);
    Rcpp::traits::input_parameter< const int >::type num_obs_all(num_obs_allSEXP);
    rcpp_result_gen = Rcpp::wrap(Log_marginal_post_no_discrepancy(param, output, p_theta, X, have_mean, inv_output_weights, cm_obs, S_2_f, num_obs_all));
    return rcpp_result_gen;
END_RCPP
}
// Update_R_inv_y
MatrixXd Update_R_inv_y(VectorXd R_inv_y, List R0, VectorXd beta_delta, String kernel_type, VectorXd alpha, double lambda_z, int num_obs);
RcppExport SEXP _RobustCalibration_Update_R_inv_y(SEXP R_inv_ySEXP, SEXP R0SEXP, SEXP beta_deltaSEXP, SEXP kernel_typeSEXP, SEXP alphaSEXP, SEXP lambda_zSEXP, SEXP num_obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< VectorXd >::type R_inv_y(R_inv_ySEXP);
    Rcpp::traits::input_parameter< List >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type beta_delta(beta_deltaSEXP);
    Rcpp::traits::input_parameter< String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_z(lambda_zSEXP);
    Rcpp::traits::input_parameter< int >::type num_obs(num_obsSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_R_inv_y(R_inv_y, R0, beta_delta, kernel_type, alpha, lambda_z, num_obs));
    return rcpp_result_gen;
END_RCPP
}
// Get_inv_all
List Get_inv_all(const List param, const VectorXd lambda_z, const VectorXi is_SGaSP, const List R0, const List kernel_type, const List alpha_list, const List p_x, const int num_sources);
RcppExport SEXP _RobustCalibration_Get_inv_all(SEXP paramSEXP, SEXP lambda_zSEXP, SEXP is_SGaSPSEXP, SEXP R0SEXP, SEXP kernel_typeSEXP, SEXP alpha_listSEXP, SEXP p_xSEXP, SEXP num_sourcesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type lambda_z(lambda_zSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type is_SGaSP(is_SGaSPSEXP);
    Rcpp::traits::input_parameter< const List >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< const List >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const List >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const List >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< const int >::type num_sources(num_sourcesSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_inv_all(param, lambda_z, is_SGaSP, R0, kernel_type, alpha_list, p_x, num_sources));
    return rcpp_result_gen;
END_RCPP
}
// Sample_delta
VectorXd Sample_delta(const List cov_inv_all, const List tilde_output_cur, const List param, const List p_x, const int num_sources, const int num_obs, const VectorXd rand_norm);
RcppExport SEXP _RobustCalibration_Sample_delta(SEXP cov_inv_allSEXP, SEXP tilde_output_curSEXP, SEXP paramSEXP, SEXP p_xSEXP, SEXP num_sourcesSEXP, SEXP num_obsSEXP, SEXP rand_normSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type cov_inv_all(cov_inv_allSEXP);
    Rcpp::traits::input_parameter< const List >::type tilde_output_cur(tilde_output_curSEXP);
    Rcpp::traits::input_parameter< const List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const List >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< const int >::type num_sources(num_sourcesSEXP);
    Rcpp::traits::input_parameter< const int >::type num_obs(num_obsSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type rand_norm(rand_normSEXP);
    rcpp_result_gen = Rcpp::wrap(Sample_delta(cov_inv_all, tilde_output_cur, param, p_x, num_sources, num_obs, rand_norm));
    return rcpp_result_gen;
END_RCPP
}
// Log_marginal_post_delta
double Log_marginal_post_delta(const Eigen::VectorXd param, Eigen::MatrixXd L, const Eigen::VectorXd delta, const int p_x, const VectorXd CL, const double a, const double b);
RcppExport SEXP _RobustCalibration_Log_marginal_post_delta(SEXP paramSEXP, SEXP LSEXP, SEXP deltaSEXP, SEXP p_xSEXP, SEXP CLSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type L(LSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const int >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type CL(CLSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(Log_marginal_post_delta(param, L, delta, p_x, CL, a, b));
    return rcpp_result_gen;
END_RCPP
}
// Log_profile_lik
List Log_profile_lik(const Eigen::VectorXd param, const String discrepancy_type, const Eigen::VectorXd output, const int p_theta, int p_x, Eigen::MatrixXd X, bool have_mean, const VectorXd cm_obs, const double lambda_z, const List R0, const String kernel_type, const Eigen::VectorXd alpha, const Eigen::VectorXd inv_output_weights, const int num_obs_all, const double S_2_f);
RcppExport SEXP _RobustCalibration_Log_profile_lik(SEXP paramSEXP, SEXP discrepancy_typeSEXP, SEXP outputSEXP, SEXP p_thetaSEXP, SEXP p_xSEXP, SEXP XSEXP, SEXP have_meanSEXP, SEXP cm_obsSEXP, SEXP lambda_zSEXP, SEXP R0SEXP, SEXP kernel_typeSEXP, SEXP alphaSEXP, SEXP inv_output_weightsSEXP, SEXP num_obs_allSEXP, SEXP S_2_fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const String >::type discrepancy_type(discrepancy_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const int >::type p_theta(p_thetaSEXP);
    Rcpp::traits::input_parameter< int >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type have_mean(have_meanSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type cm_obs(cm_obsSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_z(lambda_zSEXP);
    Rcpp::traits::input_parameter< const List >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type inv_output_weights(inv_output_weightsSEXP);
    Rcpp::traits::input_parameter< const int >::type num_obs_all(num_obs_allSEXP);
    Rcpp::traits::input_parameter< const double >::type S_2_f(S_2_fSEXP);
    rcpp_result_gen = Rcpp::wrap(Log_profile_lik(param, discrepancy_type, output, p_theta, p_x, X, have_mean, cm_obs, lambda_z, R0, kernel_type, alpha, inv_output_weights, num_obs_all, S_2_f));
    return rcpp_result_gen;
END_RCPP
}
// Loss_function_no_discrepancy
List Loss_function_no_discrepancy(const Eigen::VectorXd output, const int p_theta, Eigen::MatrixXd X, bool have_mean, const VectorXd cm_obs, const Eigen::VectorXd inv_output_weights, const int num_obs_all, const double S_2_f);
RcppExport SEXP _RobustCalibration_Loss_function_no_discrepancy(SEXP outputSEXP, SEXP p_thetaSEXP, SEXP XSEXP, SEXP have_meanSEXP, SEXP cm_obsSEXP, SEXP inv_output_weightsSEXP, SEXP num_obs_allSEXP, SEXP S_2_fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const int >::type p_theta(p_thetaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type have_mean(have_meanSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type cm_obs(cm_obsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type inv_output_weights(inv_output_weightsSEXP);
    Rcpp::traits::input_parameter< const int >::type num_obs_all(num_obs_allSEXP);
    Rcpp::traits::input_parameter< const double >::type S_2_f(S_2_fSEXP);
    rcpp_result_gen = Rcpp::wrap(Loss_function_no_discrepancy(output, p_theta, X, have_mean, cm_obs, inv_output_weights, num_obs_all, S_2_f));
    return rcpp_result_gen;
END_RCPP
}
