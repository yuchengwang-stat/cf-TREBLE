#include <RcppArmadillo.h>
#include <cmath>
#include <pbv.h>
#include <boost/math/special_functions/beta.hpp>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
//[[Rcpp::depends(RcppArmadillo,pbv,BH)]]
using namespace Rcpp;


double log_PBB_vec_cpp(double p, double alpha, arma::vec Y, arma::vec N) {
  double ap = alpha * p;
  double apc = alpha * (1 - p);
  int n = Y.size();
  double total = 0;
  NumericVector smd;
  double subs = n * R::lbeta(ap, apc);
  
  if (!arma::is_finite(subs)) {
    return 0;
  }
  NumericVector y_plus_ap = wrap(Y + ap);
  NumericVector n_minus_y_plus_apc = wrap(N - Y + apc);
  smd = lbeta(y_plus_ap,n_minus_y_plus_apc);
  for(int i=0;i<n;i++)
  {
    if (arma::is_finite(smd(i))) {
      total += smd(i);
    }
  }
  
  total -= subs;
  return total;
} 


double log_PBB_cpp(double p, double alpha, arma::vec Y, arma::vec N) {
  double ap = alpha * p;
  double apc = alpha * (1 - p);
  int n = Y.size();
  double total = 0;
  double smd = 0;
  double subs = n * R::lbeta(ap, apc);
  
  if (!arma::is_finite(subs)) {
    return 0;
  }
  
  for (int i = 0; i < n; i++) {
    smd = R::lbeta(Y[i] + ap, N[i] - Y[i] + apc);
    if (arma::is_finite(smd)) {
      total += smd;
    }
  }
  
  total -= subs;
  return total;
}


arma::vec log_PBB_1_cpp(double p, double alpha, arma::vec Y, arma::vec N) {
  int n = Y.size();
  arma::vec out(2, arma::fill::zeros); // Initialize with zeros
  double ap = alpha * p;
  double apc = alpha * (1 - p);
  double digamma_ap = R::digamma(ap);
  double digamma_apc = R::digamma(apc);
  double digamma_alpha = R::digamma(alpha);
  double p_digamma_ap = p * digamma_ap;
  double one_minus_p_digamma_apc = (1 - p) * digamma_apc;
  
  double first_term = 0;
  double second_term = 0;
  
  for (int i = 0; i < n; i++) {
    double digamma_y_ap = R::digamma(Y[i] + ap);
    double digamma_n_minus_y_apc = R::digamma(N[i] - Y[i] + apc);
    double digamma_n_alpha = R::digamma(N[i] + alpha);
    
    first_term += digamma_y_ap - digamma_n_minus_y_apc - digamma_ap + digamma_apc;
    second_term += p * digamma_y_ap + (1 - p) * digamma_n_minus_y_apc - digamma_n_alpha - p_digamma_ap - one_minus_p_digamma_apc + digamma_alpha;
  }
  
  out[0] = alpha * first_term;
  out[1] = second_term;
  
  return out;
}

arma::mat22 log_PBB_2_cpp(double p, double alpha, arma::vec Y, arma::vec N) {
  arma::mat22 out;
  int n = N.size();
  
  double ap = alpha * p;
  double apc = alpha * (1 - p);
  double alpha2 = std::pow(alpha, 2);
  double p2 = std::pow(p, 2);
  double pc2 = std::pow((1 - p), 2);
  
  double trigamma_ap = R::trigamma(ap);
  double trigamma_apc = R::trigamma(apc);
  double trigamma_alpha = R::trigamma(alpha);
  
  double digamma_ap = R::digamma(ap);
  double digamma_apc = R::digamma(apc);
  
  double diag1 = 0;
  double diag2 = 0;
  double offdiag = 0;
  
  for (int i = 0; i < n; i++) {
    double trigamma_y_ap = R::trigamma(Y[i] + ap);
    double trigamma_n_minus_y_apc = R::trigamma(N[i] - Y[i] + apc);
    double trigamma_n_alpha = R::trigamma(N[i] + alpha);
    
    double digamma_y_ap = R::digamma(Y[i] + ap);
    double digamma_n_minus_y_apc = R::digamma(N[i] - Y[i] + apc);
    
    // Diagonal elements
    diag1 += alpha2 * (trigamma_y_ap + trigamma_n_minus_y_apc - trigamma_ap - trigamma_apc);
    diag2 += p2 * trigamma_y_ap + pc2 * trigamma_n_minus_y_apc - trigamma_n_alpha - p2 * trigamma_ap - pc2 * trigamma_apc + trigamma_alpha;
    
    // Off-diagonal elements
    offdiag += digamma_y_ap + ap * trigamma_y_ap - digamma_n_minus_y_apc - apc * trigamma_n_minus_y_apc - digamma_ap - ap * trigamma_ap + digamma_apc + apc * trigamma_apc;
  }
  
  out(0, 0) = diag1;
  out(1, 1) = diag2;
  out(0, 1) = offdiag;
  out(1, 0) = offdiag;
  
  return out;
}


double log_prior_cpp(double mu, double sigma, arma::vec delta, arma::vec tau,arma::vec lambda,double nu, double phi)
{
  double n = delta.size();
  arma::vec smd(n);
  for(int i=0;i<n;i++)
  {
    smd[i] = std::log(lambda[i]) + R::dnorm(mu,delta[i],tau[i],1);
  }
  double C = max(smd);
  arma::vec exp_smd_shifted(n);
  for(int i=0;i<n;i++)
  {
    exp_smd_shifted[i] =  std::exp(smd[i] - C);
  }
  double muprior = std::log(sum(exp_smd_shifted)) + C;
  double out = muprior + std::log(2 * sigma) + nu/2 * std::log(phi * nu /2 ) - 
    R::lgammafn(nu/2) - 
    phi * nu /(2*std::pow(sigma,2))-
    2*(1+nu/2)*std::log(sigma);
  return(out);
}
  
arma::vec log_prior_1_cpp(double mu, double sigma, arma::vec delta, arma::vec tau, arma::vec lambda, double nu, double phi) {
  arma::vec out(2);
  int n = delta.n_elem;

  double log_2_sigma = std::log(2 * sigma);
  double half_nu_log_phi_nu_over_2 = nu / 2 * std::log(phi * nu / 2);
  double lgamma_nu_over_2 = R::lgammafn(nu / 2);
  double phi_nu_over_2_sigma_sq = phi * nu / (2 * std::pow(sigma, 2));
  double log_prior_base = log_prior_cpp(mu, sigma, delta, tau, lambda, nu, phi);
  double log_prior_adjustment = log_prior_base - log_2_sigma - half_nu_log_phi_nu_over_2 + lgamma_nu_over_2 + phi_nu_over_2_sigma_sq + 2 * (1 + nu / 2) * std::log(sigma);
  double de = std::exp(log_prior_adjustment);
  
  arma::vec smd = mu - delta; 
  arma::vec abs_smd = arma::abs(smd);
  arma::vec log_lambda_tau = log(lambda) - 2 * log(tau); 
  arma::vec dnorm_values = arma::vec(n);
  
  for (int i = 0; i < n; ++i) {
    dnorm_values[i] = R::dnorm(mu, delta[i], tau[i], 1);
  }
  
  arma::vec numer = log_lambda_tau + log(abs_smd) + dnorm_values;
  double C = arma::max(numer);
  numer -= C;
  
  // Weight contributions by sign of (mu-delta)
  arma::vec contributions = arma::sign(smd) % exp(numer);
  double total = arma::accu(contributions) * -exp(C);
  
  out[0] = total / de;
  out[1] = 1 / sigma + phi * nu * std::pow(sigma, -3) - 2 * (1 + nu / 2) / sigma;
  
  return out;
}




arma::mat22 log_prior_2_cpp(double mu, double sigma, arma::vec delta, arma::vec tau, arma::vec lambda, double nu, double phi){
  int n = delta.size();
  arma::vec tmp_vec(n);
  arma::vec log_tau_sq = 2 * log(tau); 
  double log_2_sigma = std::log(2 * sigma);
  double half_nu_log_phi_nu_over_2 = nu / 2 * std::log(phi * nu / 2);
  double lgamma_nu_over_2 = R::lgammafn(nu / 2);
  double phi_nu_over_2_sigma_sq = phi * nu / (2 * std::pow(sigma, 2));
  double log_prior_base = log_prior_cpp(mu, sigma, delta, tau, lambda, nu, phi) - log_2_sigma - half_nu_log_phi_nu_over_2 + lgamma_nu_over_2 + phi_nu_over_2_sigma_sq + 2 * (1 + nu / 2) * std::log(sigma);
  double de = std::exp(log_prior_base);
  
  for(int i = 0; i < n; i++) {
    tmp_vec[i] = R::dnorm(mu, delta[i], tau[i], 1) + std::log(lambda[i]) - log_tau_sq[i];
  }
  double tmp_C = max(tmp_vec);
  tmp_vec = tmp_vec - tmp_C;
  double tmp_total = std::exp(tmp_C) * sum(exp(tmp_vec));
  
  arma::vec numer(n);
  for(int i = 0; i < n; i++) {
    numer[i] = std::log(lambda[i]) + 2 * std::log(std::abs(mu - delta[i])) - 4 * log(tau[i]) + R::dnorm(mu, delta[i], tau[i], 1);
  }
  double C = max(numer);
  numer = numer - C;
  double total = std::exp(C) * sum(exp(numer));
  
  double out = std::log(total);
  double log_prior1 = log_prior_1_cpp(mu, sigma, delta, tau, lambda, nu, phi)[0] * de;
  out = ((std::exp(out) - tmp_total) * de - std::pow(log_prior1, 2)) / std::pow(de, 2);
  
  arma::mat22 result;
  result.zeros(); 
  result(0, 0) = out;
  result(1, 1) = -1 / std::pow(sigma, 2) - 3 * phi * nu * std::pow(sigma, -4) + 2 * (1 + nu / 2) / std::pow(sigma, 2);
  
  return result;
}


arma::mat22 J_cpp(double mu, double sigma) {
  double sigma_squared = sigma * sigma;
  double denom = 1 + sigma_squared;
  double norm_term = mu / std::sqrt(denom);
  
  double p = R::pnorm(norm_term, 0, 1, true, false);
  double q = pbv::pbv_rcpp_pbvnorm0(norm_term, norm_term, sigma_squared / denom);
  
  double sqrt_denom_inv = 1 / std::sqrt(denom);
  double dpdmu = sqrt_denom_inv * R::dnorm(norm_term, 0, 1, false);
  double dpdsigma = -sigma * mu / (denom / sqrt_denom_inv) * R::dnorm(norm_term, 0, 1, false);
  
  double sqrt_1_plus_sigma_2 = std::sqrt(1 + sigma_squared);
  double sqrt_denom_2_sigma = std::sqrt(denom * (1 + 2 * sigma_squared));
  double partial_mu = mu / sqrt_denom_2_sigma;
  double dnorm_partial_mu = R::dnorm(partial_mu, 0, 1, false);
  
  double dqdmu = 2 * R::dnorm(mu, 0, sqrt_1_plus_sigma_2, false) * R::pnorm(partial_mu, 0, 1, true, false);
  double dqdsigma = 2 * R::dnorm(mu, 0, sqrt_1_plus_sigma_2, false) * 
    ((sigma / sqrt_denom_2_sigma) * dnorm_partial_mu -
    (mu * sigma) / denom * R::pnorm(partial_mu, 0, 1, true, false));
  
  double q_minus_p_squared = q - p * p;
  double p_one_minus_p = p * (1 - p);
  double dalphadmu = (dpdmu - 2 * p * dpdmu) / q_minus_p_squared - p_one_minus_p / std::pow(q_minus_p_squared, 2) * (dqdmu - 2 * p * dpdmu);
  double dalphadsigma = (dpdsigma - 2 * p * dpdsigma) / q_minus_p_squared - p_one_minus_p / std::pow(q_minus_p_squared, 2) * (dqdsigma - 2 * p * dpdsigma);
  
  arma::mat22 J;
  J(0, 0) = dpdmu;
  J(1,0) = dpdsigma;
  J(0,1) = dalphadmu;
  J(1, 1) = dalphadsigma;
  
  return J;
}

List J_grad_cpp(double mu, double sigma) {
  // Common expressions
  double s2 = sigma * sigma;
  double one_plus_s2 = 1.0 + s2;
  double one_plus_2s2 = 1.0 + 2.0 * s2;
  double sqrt_one_plus_s2 = sqrt(one_plus_s2);
  double sqrt_one_plus_2s2 = sqrt(one_plus_2s2);
  double one_plus_s2_inv_sqrt = 1.0 / sqrt_one_plus_s2;
  double one_plus_s2_pow_neg_3_2 = pow(one_plus_s2, -1.5);
  double one_plus_s2_pow_neg_5_2 = pow(one_plus_s2, -2.5);
  double one_plus_s2_pow_neg_7_2 = pow(one_plus_s2, -3.5);
  double sd1 = sqrt_one_plus_s2;
  double sd2 = sqrt_one_plus_2s2;
  double sd1_sd2 = sd1 * sd2;
  double sd1_sq = one_plus_s2;
  double sd2_sq = one_plus_2s2;
  double sd1_fourth = sd1_sq * sd1_sq;
  double sd1_sd2_sq = sd1_sd2 * sd1_sd2;
  double x = mu * one_plus_s2_inv_sqrt;
  double dnorm_x = R::dnorm(x, 0.0, 1.0, false);
  double p = R::pnorm(x, 0.0, 1.0, true, false);
  double rho = s2 / one_plus_s2;
  double q = pbv::pbv_rcpp_pbvnorm0(x, x, rho);
  double dp2dsigma2 = (-mu * one_plus_s2_pow_neg_3_2 + 3.0 * mu * s2 * one_plus_s2_pow_neg_5_2) * dnorm_x -
    one_plus_s2_pow_neg_7_2 * mu * mu * mu * s2 *dnorm_x;
  double dp2dmu2 = -mu * one_plus_s2_pow_neg_3_2 * dnorm_x;
  double dp2dmudsigma = (-sigma * one_plus_s2_pow_neg_3_2 + sigma * mu * mu * one_plus_s2_pow_neg_5_2) * dnorm_x;
  double dpdmu = one_plus_s2_inv_sqrt * dnorm_x;
  double dpdsigma = -sigma * mu * one_plus_s2_pow_neg_3_2 * dnorm_x;
  double z1 = mu / sd1;
  double dnorm_z1 = R::dnorm(z1, 0.0, 1.0, false) / sd1;
  double z2 = mu / (sd1_sd2);
  double pnorm_z2 = R::pnorm(z2, 0.0, 1.0, true, false);
  double dqdmu = 2.0 * dnorm_z1 * pnorm_z2;
  double z3 = mu / sqrt(sd1_sq * sd2_sq);
  double dnorm_z3 = R::dnorm(z3, 0.0, 1.0, false);
  double term1 = sigma / sqrt(sd1_sq * sd2_sq) * dnorm_z3;
  double term2 = (mu * sigma) / sd1_sq * pnorm_z2;
  double dqdsigma = 2.0 * dnorm_z1 * (term1 - term2);
  double dnorm_mu_sd1 = R::dnorm(mu, 0.0, sd1, false);
  double pnorm_mu_over_sd1sd2 = R::pnorm(mu / sd1_sd2, 0.0, 1.0, true, false);
  double dnorm_mu_over_sd1sd2 = R::dnorm(mu / sd1_sd2, 0.0, 1.0, false);
  double mu_over_sd1sd2 = mu/sd1_sd2;
  double dq2dmu2 = -2.0 * mu / sd1_sq * dnorm_mu_sd1 * pnorm_mu_over_sd1sd2 +
    2.0 * dnorm_mu_sd1 / sd1_sd2 * dnorm_mu_over_sd1sd2;
  double coeff1 = -2.0 * mu / sd1_sq * dnorm_mu_sd1;
  double term_a = sigma / sd1_sd2 * dnorm_mu_over_sd1sd2;
  double term_b = mu * sigma / sd1_sq * pnorm_mu_over_sd1sd2;
  double inner_term = term_a - term_b;
  double dq2dmudsigma = coeff1 * inner_term +
    2.0 * dnorm_mu_sd1 * (
        (-sigma / (sd1_sd2_sq) * mu / sd1_sd2 * dnorm_mu_over_sd1sd2)
    - sigma / sd1_sq * pnorm_mu_over_sd1sd2
  - sigma * mu / (sd1_sq * sd1_sd2) * dnorm_mu_over_sd1sd2
    );
  term1 = 2.0 * ( (mu * mu * sigma / sd1_fourth - sigma / (sd1_sq)) * dnorm_mu_sd1 );
  term2 = (sigma / (sd1_sd2) * dnorm_mu_over_sd1sd2) - (mu * sigma / (sd1_sq) * pnorm_mu_over_sd1sd2);
  double first_part = term1 * term2;
  
  // Compute auxiliary terms for the second part
  double denom = std::sqrt(1 + 3 * s2 + 2 * s2 *s2);
  double numerator = sd1_sd2 - 0.5 * sigma * (6 * sigma + 8 * s2*sigma) / denom;
  double C = numerator / (sd1_sd2 * sd1_sd2);
  
  double F = (sigma / (sd1_sd2)) * 0.5 * mu * (6 * sigma + 8 * s2*sigma) / (denom * denom * denom);
  F *= mu_over_sd1sd2 * dnorm_mu_over_sd1sd2;
  
  double G = (mu * sd1_sq - 2 * s2 * mu) / std::pow(sd1, 4);
  double H = G * pnorm_mu_over_sd1sd2;
  
  double I = sigma * mu / (sd1_sq) * 0.5 * mu * (6 * sigma + 8 * s2*sigma) / (denom * denom * denom);
  I *= dnorm_mu_over_sd1sd2;
  
  // Second part of the expression
  double second_part = 2.0 * dnorm_mu_sd1 * ( C * dnorm_mu_over_sd1sd2 + F - H + I );
  
  // Total result
  double dq2dsigma2 = first_part + second_part;
  denom = (q - p * p);
  double denom_sq = denom * denom;
  double denom_cu = denom_sq * denom;
  double denom_qu = denom_cu * denom;
  double dalpha2dmu2 = ((dp2dmu2 - 2 * dpdmu * dpdmu - 2 * p * dp2dmu2) * denom -
                        (dpdmu - 2 * p * dpdmu) * (dqdmu - 2 * p * dpdmu)) / denom_sq -
                        ((dpdmu - 2 * p * dpdmu) * denom_sq - 2 * denom * (dqdmu - 2 * p * dpdmu) * p * (1 - p)) / (denom_qu) * (dqdmu - 2 * p * dpdmu) -
                        (dq2dmu2 - 2 * dpdmu * dpdmu - 2 * p * dp2dmu2) * p * (1 - p) / denom_sq;
                        
  double dalpha2dmudsigma = ((dp2dmudsigma - 2 * dpdmu * dpdsigma - 2 * p * dp2dmudsigma) * denom -
                             (dpdsigma - 2 * p * dpdsigma) * (dqdmu - 2 * p * dpdmu)) / denom_sq -
                             (((dpdmu - 2 * p * dpdmu) * denom_sq - 2 * denom * (dqdmu - 2 * p * dpdmu) * p * (1 - p)) * (dqdsigma - 2 * p * dpdsigma)) / (denom_qu) -
                             (dq2dmudsigma - 2 * dpdsigma * dpdmu - 2 * p * dp2dmudsigma) * p * (1 - p) / denom_sq;
  
  double dalpha2dsigma2 = ((dp2dsigma2 - 2 * dpdsigma * dpdsigma - 2 * p * dp2dsigma2) * denom -
                           (dpdsigma - 2 * p * dpdsigma) * (dqdsigma - 2 * p * dpdsigma)) / denom_sq -
                           ((dpdsigma - 2 * p * dpdsigma) * denom_sq - 2 * denom * (dqdsigma - 2 * p * dpdsigma) * p * (1 - p)) / (denom_qu) * (dqdsigma - 2 * p * dpdsigma) -
                           (dq2dsigma2 - 2 * dpdsigma * dpdsigma - 2 * p * dp2dsigma2) * p * (1 - p) / denom_sq;
  
  // Assemble the matrices
  arma::mat dJdmu(2,2);
  arma::mat dJdsigma(2,2);
  
  dJdmu(0, 0) = dp2dmu2;
  dJdmu(0, 1) = dp2dmudsigma;
  dJdmu(1, 0) = dalpha2dmu2;
  dJdmu(1, 1) = dalpha2dmudsigma;
  
  dJdsigma(0, 0) = dp2dmudsigma;
  dJdsigma(0, 1) = dp2dsigma2;
  dJdsigma(1, 0) = dalpha2dmudsigma;
  dJdsigma(1, 1) = dalpha2dsigma2;
  //return List::create(dalpha2dmu2, dalpha2dsigma2,dalpha2dmudsigma,dq2dmu2,dq2dsigma2,dq2dmudsigma,dpdmu,dpdsigma,dqdmu,dqdsigma,q,p);
  return List::create(dJdmu, dJdsigma);
}

double f_gj_tilde_cpp(double mu, double theta, arma::vec Y, arma::vec N, arma::vec delta, arma::vec tau, arma::vec lambda,double nu,double phi)
{
  double sigma = std::exp(theta);
  double s2 = sigma * sigma;
  double x = mu * pow(1+s2,-0.5);
  double p = R::pnorm(x,0,1,1,0);
  double rho = s2/(1+s2);
  double q = pbv::pbv_rcpp_pbvnorm0(x, x, rho);
  double alpha = p * (1-p) / (q-p * p) -1;
  double out = log_PBB_cpp(p,alpha,Y,N) + log_prior_cpp(mu,sigma,delta,tau,lambda,nu,phi) + theta;
  return out;
}


arma::vec f_gj_1_cpp(double mu, double sigma, arma::vec Y, arma::vec N, arma::vec delta, arma::vec tau, arma::vec lambda,double nu,double phi)
{
  double s2 = sigma * sigma;
  double x = mu * pow(1+s2,-0.5);
  double p = R::pnorm(x,0,1,1,0);
  double rho = s2/(1+s2);
  double q = pbv::pbv_rcpp_pbvnorm0(x, x, rho);
  double alpha = p * (1-p) / (q-p * p) -1;
  arma::mat22 J = J_cpp(mu,sigma);
  arma::vec log_PBB_1 = log_PBB_1_cpp(p,alpha,Y,N);
  arma::vec log_prior_1 = log_prior_1_cpp(mu,sigma,delta,tau,lambda,nu,phi);
  arma::vec out = J * log_PBB_1 +  log_prior_1;
  return out;
}


arma::vec f_gj_1_tilde_cpp(double mu, double theta, arma::vec Y, arma::vec N, arma::vec delta, arma::vec tau, arma::vec lambda,double nu,double phi)
{
  double sigma = std::exp(theta);
  arma::mat tmp_mat = arma::zeros<arma::mat>(2,2);
  tmp_mat(0,0) = 1;
  tmp_mat(1,1) = sigma;
  arma::vec f_gj_1 = f_gj_1_cpp(mu,sigma,Y,N,delta,tau,lambda,nu,phi);
  arma::vec out = tmp_mat * f_gj_1;
  out[1] += 1;
  return out;
  
}


arma::mat22 f_gj_2_tilde_cpp(double mu, double theta, arma::vec Y, arma::vec N, arma::vec delta, arma::vec tau, arma::vec lambda,double nu,double phi)
{
  double avg = 0;
  double sigma = std::exp(theta);
  double s2 = sigma * sigma;
  double x = mu * pow(1+s2,-0.5);
  double p = R::pnorm(x,0,1,1,0);
  double rho = s2/(1+s2);
  double q = pbv::pbv_rcpp_pbvnorm0(x, x, rho);
  double alpha = p * (1-p) / (q-p * p) -1;
  arma::mat22 J = J_cpp(mu,sigma);
  List J_grad = J_grad_cpp(mu,sigma);
  arma::mat22 log_PBB_2 = log_PBB_2_cpp(p,alpha,Y,N);
  arma::vec log_PBB_1 = log_PBB_1_cpp(p,alpha,Y,N);
  arma::mat22 log_prior_2 = log_prior_2_cpp(mu,sigma,delta,tau,lambda,nu,phi);
  arma::mat22 mid;
  arma::mat dJdmu = J_grad[0];
  arma::mat dJdsigma = J_grad[1];
  dJdmu = dJdmu.t();
  dJdsigma = dJdsigma.t();
  mid.col(0) = dJdmu * log_PBB_1;
  mid.col(1) = dJdsigma * log_PBB_1;
  arma::mat22 out = J * log_PBB_2 * J.t() + mid + log_prior_2;
  arma::vec f_gj_1 = f_gj_1_cpp(mu,sigma,Y,N,delta,tau,lambda,nu,phi);
  out(1,0) = sigma * out(1,0);
  out(0,1) = sigma * out(0,1);
  avg = (out(1,0) + out(0,1))/2;
  out(1,0) = avg;
  out(0,1) = avg;
  out(1,1) = s2 * out(1,1) + sigma * f_gj_1[1];
  return out;
}

bool isPositiveDefinite(arma::mat& mat) {
  arma::mat L;
  return chol(L, mat);
}

double target_fn(arma::vec para,arma::vec current_Y,arma::vec current_N,arma::vec delta,arma::vec tau,arma::vec lambda,double nu,double phi)
{
  double out = f_gj_tilde_cpp(para[0],para[1],current_Y,current_N,delta,tau,lambda,nu,phi);
  return out;
}

arma::vec target_fn_grad(arma::vec para,arma::vec current_Y,arma::vec current_N,arma::vec delta,arma::vec tau,arma::vec lambda,double nu,double phi)
{
  arma::vec out = f_gj_1_tilde_cpp(para[0],para[1],current_Y,current_N,delta,tau,lambda,nu,phi);
  return out;
}




arma::mat matrixSquareRootInverse(arma::mat& mat) {
  // Ensure the matrix is symmetric and positive definite
  if(!mat.is_square() || arma::approx_equal(mat, mat.t(), "reldiff", 1e-5) == false) {
    Rcpp::stop("Input must be a symmetric and square matrix.");
  }
  
  // Compute the eigen decomposition
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, mat);
  
  // Compute the square root of eigenvalues (only if they are all positive)
  if(any(eigval <= 0)) {
    Rcpp::stop("Matrix is not positive definite.");
  }
  arma::vec sqrt_eigval = arma::sqrt(eigval);
  
  // Compute the square root matrix
  arma::mat sqrt_mat = eigvec * arma::diagmat(sqrt_eigval) * eigvec.t();
  
  // Invert the square root matrix
  arma::mat sqrt_mat_inv = arma::inv(sqrt_mat);
  
  return sqrt_mat_inv;
}

List optim_rcpp(arma::vec& init_val,arma::vec& current_Y,arma::vec& current_N,arma::vec& delta,arma::vec& tau,arma::vec& lambda,double& nu,double& phi ){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];    
  List control = List::create(Named("fnscale") = -1,
                              Named("maxit") = 1000);
  arma::vec lower(2);
  arma::vec upper(2);
  lower[0] = -2.5;
  lower[1] = R_NegInf;
  upper[0] = 2.5;
  upper[1] = R_PosInf;
  Rcpp::List opt_results = optim(_["par"]    = init_val,
                                 _["fn"]     = Rcpp::InternalFunction(&target_fn),
                                 _["gr"]     = Rcpp::InternalFunction(&target_fn_grad),
                                 _["method"] = "L-BFGS-B",
                                 _["current_Y"] = current_Y,
                                 _["current_N"] = current_N,
                                 _["delta"] = delta,
                                 _["tau"] = tau,
                                 _["lambda"] = lambda,
                                 _["nu"] = nu,
                                 _["phi"] = phi,
                                 _["lower"] = lower,
                                 _["upper"] = upper,
                                 _["control"] = control
                                  );
  arma::vec para_optim = opt_results[0];
  arma::mat H = - f_gj_2_tilde_cpp(para_optim[0],para_optim[1],current_Y,current_N,delta,tau,lambda,nu,phi);
  bool pd = isPositiveDefinite(H);
  while(para_optim(0)>(upper[0]-0.02) | para_optim(0)<(lower[0]+0.02) | !pd)
  {
    lower[0] -= 0.5;
    upper[0] += 0.5;
    if(upper[0] == 5)
    {
      init_val[0] = 0;
    }
    opt_results = optim(_["par"]    = init_val,
                        _["fn"]     = Rcpp::InternalFunction(&target_fn),
                        _["gr"]     = Rcpp::InternalFunction(&target_fn_grad),
                        _["method"] = "L-BFGS-B",
                        _["current_Y"] = current_Y,
                       _["current_N"] = current_N,
                       _["delta"] = delta,
                       _["tau"] = tau,
                       _["lambda"] = lambda,
                       _["nu"] = nu,
                       _["phi"] = phi,
                       _["lower"] = lower,
                       _["upper"] = upper,
                       _["control"] = control);
    para_optim = as<arma::vec>(opt_results[0]);
    arma::mat H = - f_gj_2_tilde_cpp(para_optim[0],para_optim[1],current_Y,current_N,delta,tau,lambda,nu,phi);
    pd = isPositiveDefinite(H);
    //Rcpp::Rcout << H  << std::endl;
    if(upper[0] == 6)
    {
      break;
    }
  }
  H = - f_gj_2_tilde_cpp(para_optim[0],para_optim[1],current_Y,current_N,delta,tau,lambda,nu,phi);
  pd = isPositiveDefinite(H);
  if(!pd)
  {
    opt_results = optim(_["par"]    = init_val,
                        _["fn"]     = Rcpp::InternalFunction(&target_fn),
                        _["gr"]     = Rcpp::InternalFunction(&target_fn_grad),
                        _["method"] = "BFGS",
                        _["current_Y"] = current_Y,
                        _["current_N"] = current_N,
                        _["delta"] = delta,
                        _["tau"] = tau,
                        _["lambda"] = lambda,
                        _["nu"] = nu,
                        _["phi"] = phi,
                        _["control"] = control
    );
    para_optim = as<arma::vec>(opt_results[0]);
  }
  H = - f_gj_2_tilde_cpp(para_optim[0],para_optim[1],current_Y,current_N,delta,tau,lambda,nu,phi);
  //Rcpp::Rcout << H  << std::endl;
  pd = isPositiveDefinite(H);
  if(!pd)
  {
    opt_results = optim(_["par"]    = init_val,
                        _["fn"]     = Rcpp::InternalFunction(&target_fn),
                        _["method"] = "Nelder-Mead",
                        _["current_Y"] = current_Y,
                        _["current_N"] = current_N,
                        _["delta"] = delta,
                        _["tau"] = tau,
                        _["lambda"] = lambda,
                        _["nu"] = nu,
                        _["phi"] = phi,
                        _["control"] = control
    );
    
 }
  return opt_results;
}



arma::vec callTildeFGJ(NumericVector mu,NumericVector theta,NumericVector Y, NumericVector N,
                       arma::vec delta, arma::vec tau, arma::vec lambda, double nu, double phi) {
  Environment env = Environment::global_env();  
  Function tilde_f_gj = env["tilde.f.gj"];  
  return as<arma::vec>(tilde_f_gj(mu, theta, Y, N, delta, tau, lambda, nu, phi));
}

// [[Rcpp::export]]
arma::mat log_lik_laplace_cpp(List leaf_tree, List Y_all, List N_all, arma::vec delta, arma::vec tau, arma::vec lambda,double nu,double phi,arma::mat quad_x, arma::vec quad_w)
{
  double res;
  double C;
  List opt_result;
  const double pi =  3.1415926535;
  double log2pi = log(2*pi);
  arma::vec opt_para;
  double opt_value;
  int n_nodes = leaf_tree.length();
  arma::vec current_node;
  double n_cell_types;
  int current_cell_type;
  arma::vec N_cell_type;
  arma::vec Y_cell_type;
  arma::vec sigma_eval(100);
  arma::vec mu_eval(100);
  double qnormp0 = 0;
  arma::mat H_sqrt_inv;
  arma::mat H;
  arma::vec log_lik(n_nodes);
  arma::vec empty;
  arma::vec current_N;
  arma::vec current_Y;
  arma::mat quad_x_mod;
  int n_quad2 = quad_w.size();
  int tmp_size=0;
  int pos = 0;
  int len = 0;
  double mu;
  double theta;
  double avg;
  arma::vec para;
  arma::mat intermediate(2,n_quad2);
  arma::vec smd(n_quad2);
  arma::vec mu_vec(n_nodes);
  arma::vec theta_vec(n_nodes);
  arma::mat vec_list(n_nodes,3);
  for(int k=0;k<10;k++)
  {
    sigma_eval(k) = k * 10;
  }
  sigma_eval = sigma_eval/80 + 0.2;
  arma::vec sqrt_one_plus_sigma_eval_sq = arma::pow(1+arma::pow(sigma_eval,2),-0.5);
  arma::vec theta_eval = arma::log(sigma_eval);
  for(int i=0;i<n_nodes;i++)
  {
    current_node = as<arma::vec>(leaf_tree[i])-1;
    n_cell_types = current_node.size();
    len = 0;
    for(int m =0;m<n_cell_types;m++)
    {
      current_cell_type = current_node[m];
      len += as<arma::vec>(N_all[current_cell_type]).size();
    }
    current_N.resize(len);
    current_Y.resize(len);
    pos = 0;
    for(int j=0;j<n_cell_types;j++)
    {
      current_cell_type = current_node[j];
      N_cell_type = as<arma::vec>(N_all[current_cell_type]);
      Y_cell_type = as<arma::vec>(Y_all[current_cell_type]);
      tmp_size = Y_cell_type.size();
      current_N.subvec(pos,pos+tmp_size-1) = N_cell_type;
      current_Y.subvec(pos,pos+tmp_size-1) = Y_cell_type;
      pos += tmp_size;
    }
   qnormp0 = R::qnorm((sum(current_Y) + 1.06) / (sum(current_N) + 1.52),0,1,1,0);
   mu_eval = qnormp0 * sqrt_one_plus_sigma_eval_sq;
   double tmp_max = -999999;
   double index_max;
   double tmp;
   for(int l = 0;l<100;l++)
   {
     tmp = f_gj_tilde_cpp(mu_eval[l],theta_eval[l],current_Y,current_N,delta,tau,lambda,nu,phi);
     if(tmp>tmp_max)
     {
       index_max = l;
       tmp_max = tmp;
     }
   }
   arma::vec init_para(2);
   init_para[0] = mu_eval[index_max];
   init_para[1] = theta_eval[index_max];
   opt_result = optim_rcpp(init_para,current_Y,current_N,delta,tau,lambda,nu,phi);
   opt_para = as<arma::vec>(opt_result[0]);
   opt_value = opt_result[1];
   H = - f_gj_2_tilde_cpp(opt_para[0],opt_para[1],current_Y,current_N,delta,tau,lambda,nu,phi);
   avg = (H(0,1) + H(1,0))/2;
   H(0,1) = avg;
   H(1,0) = avg;
   if(!isPositiveDefinite(H))
   {
     for(int c = 0;c<n_nodes;c++)
     {
       log_lik[c] = 999;
       mu_vec[c] = 999;
       theta_vec[c] = 999;
     }
     vec_list.col(0) = log_lik;
     vec_list.col(1) = mu_vec;
     vec_list.col(2) = theta_vec;
     return vec_list;
   }
   H_sqrt_inv = matrixSquareRootInverse(H);
   intermediate = H_sqrt_inv * quad_x;
   intermediate.each_col() += opt_para;
   quad_x_mod = intermediate;
   
   for(int t = 0; t<n_quad2;t++)
   {
     para = quad_x_mod.col(t);
    mu =para(0);
     theta = para(1);
     smd(t) =f_gj_tilde_cpp(mu,theta,current_Y,current_N,delta,tau,lambda,nu,phi);
   }

   smd = quad_w + smd - opt_value - log2pi;
   C = max(smd);
   res = log(sum(exp(smd - C))) + C;
   
   log_lik[i] = opt_value - 0.5*log(arma::det(H)) +log2pi + res;
   mu_vec[i] = opt_para[0];
   theta_vec[i] = opt_para[1];
  }
  vec_list.col(0) = log_lik;
  vec_list.col(1) = mu_vec;
  vec_list.col(2) = theta_vec;
  return vec_list;
}

List estimate_para1(List leaf_tree,List parent_tree, List Y_all, List N_all, arma::vec delta,
                   arma::vec tau, arma::vec lambda,double nu,double phi,
                   arma::mat quad_x, arma::vec quad_w,arma::vec log_lik,arma::vec mu_opt,
                   arma::vec theta_opt,arma::vec prob_vec)
{
  double res;
  double C;
  double res1;
  double res2;
  double C1;
  double C2;
  List opt_result;
  const double pi =  3.1415926535;
  double log2pi = log(2*pi);
  arma::vec opt_para(2);
  double opt_value;
  int n_nodes = leaf_tree.length();
  arma::vec current_node;
  double n_cell_types;
  int current_cell_type;
  arma::vec N_cell_type;
  arma::vec Y_cell_type;
  arma::mat H_sqrt_inv;
  arma::mat H;
  arma::vec empty;
  arma::vec current_N;
  arma::vec current_Y;
  arma::mat quad_x_mod;
  int n_quad2 = quad_w.size();
  int tmp_size=0;
  int pos = 0;
  int len = 0;
  double mu;
  double theta;
  arma::vec para;
  arma::mat intermediate(2,n_quad2);
  arma::vec smd(n_quad2);
  arma::vec smd1(n_quad2);
  arma::vec smd2(n_quad2);
  arma::vec mu_tilde(n_nodes);
  arma::vec theta_tilde(n_nodes);
  arma::vec smd1_sign(n_quad2);
  arma::vec smd2_sign(n_quad2);
  int res1_sign;
  int res2_sign;
  arma::vec weighted_smd;
  arma::uvec current_node1;
  arma::vec mu_out(n_nodes);
  arma::vec theta_out(n_nodes);
  for(int i=0;i<n_nodes;i++)
  {
    current_node = as<arma::vec>(leaf_tree[i])-1;
    n_cell_types = current_node.size();
    len = 0;
    for(int m =0;m<n_cell_types;m++)
    {
      current_cell_type = current_node[m];
      len += as<arma::vec>(N_all[current_cell_type]).size();
    }
    current_N.resize(len);
    current_Y.resize(len);
    pos = 0;
    for(int j=0;j<n_cell_types;j++)
    {
      current_cell_type = current_node[j];
      N_cell_type = as<arma::vec>(N_all[current_cell_type]);
      Y_cell_type = as<arma::vec>(Y_all[current_cell_type]);
      tmp_size = Y_cell_type.size();
      current_N.subvec(pos,pos+tmp_size-1) = N_cell_type;
      current_Y.subvec(pos,pos+tmp_size-1) = Y_cell_type;
      pos += tmp_size;
    }

    H = - f_gj_2_tilde_cpp(mu_opt(i),theta_opt(i),current_Y,current_N,delta,tau,lambda,nu,phi);
    H_sqrt_inv = matrixSquareRootInverse(H);
    intermediate = H_sqrt_inv * quad_x;

    opt_para(0) = mu_opt(i);
    opt_para(1) = theta_opt(i);
    intermediate.each_col() += opt_para;
    quad_x_mod = intermediate;
    for(int t = 0; t<n_quad2;t++)
    {
      para = quad_x_mod.col(t);
      mu =para(0);
      theta = para(1);
      smd(t) = f_gj_tilde_cpp(mu,theta,current_Y,current_N,delta,tau,lambda,nu,phi);
      smd1_sign(t) = arma::sign(mu);
      smd2_sign(t) = arma::sign(theta);
      smd1(t) = smd(t) +  log(abs(mu));
      smd2(t) = smd(t) +  log(abs(theta));
    }
    opt_value = f_gj_tilde_cpp(mu_opt(i),theta_opt(i),current_Y,current_N,delta,tau,lambda,nu,phi);
    //Rcpp::Rcout << smd1  << std::endl;
    smd1 = quad_w + smd1 - opt_value - log2pi;
    smd2 = quad_w + smd2 - opt_value - log2pi;
    C2 = max(smd2);
    C1 = max(smd1);
    //smd1_sign = arma::sign(smd1);
    //smd2_sign = arma::sign(smd2);
    //smd1 = arma::abs(smd1);
    //smd2 = arma::abs(smd2);

    //res2 = sum(smd1_sign %exp(smd2 - C2));
    //res1 =  sum(smd2_sign %exp(smd1 - C1));
    //res1_sign = arma::sign(res1);
    //res2_sign = arma::sign(res2);
    //res1 = abs(res1);
    //res2 = abs(res2);
   
    //res1 = res1_sign * log(res1) + C1;
    //Rcpp::Rcout << res1_sign  << std::endl;
    //res2 = res2_sign * log(res2) + C2;
    res1 = exp(C1) * sum(smd1_sign % exp(smd1 - C1));
    res2 = exp(C2) * sum(smd2_sign % exp(smd2 - C2));
    //res1_sign = arma::sign(res1);
    //res2_sign = arma::sign(res2);
    //res1 = abs(res1);
    //res2 = abs(res2);
    //res2 = sum(smd2_sign % exp(smd2 - C2));
    //res1 = log(res1) + C1;
    //res2 = log(res2) + C2;
    for(int t = 0; t<n_quad2;t++)
    {
      para = quad_x_mod.col(t);
      mu =para(0);
      theta = para(1);
      smd(t) =f_gj_tilde_cpp(mu,theta,current_Y,current_N,delta,tau,lambda,nu,phi);
    }
    smd = quad_w + smd - opt_value - log2pi;
    C = max(smd);
    res = log(sum(exp(smd - C))) + C;
    mu_tilde(i) = res1/exp(res);
    theta_tilde(i) = res2/exp(res);
    
  }
  for(int k=0;k<n_nodes;k++)
  {  
    current_node1 = as<arma::uvec>(parent_tree[k])-1;
    weighted_smd = prob_vec.elem(current_node1) % mu_tilde.elem(current_node1);
    mu_out(k) = sum(weighted_smd);
    weighted_smd = prob_vec.elem(current_node1) % theta_tilde.elem(current_node1);
    theta_out(k) = sum(weighted_smd);
  }
  Rcpp::List vec_list = Rcpp::List::create(
    Rcpp::Named("mu_out")  = mu_out,
    Rcpp::Named("theta_out") = theta_out
  );
  return vec_list;
}

// [[Rcpp::export]]
arma::mat estimate_para(List leaf_tree,List parent_tree, List Y_all, List N_all, arma::vec delta,
                   arma::vec tau, arma::vec lambda,double nu,double phi,
                   arma::mat quad_x, arma::vec quad_w,arma::vec mu_opt,
                   arma::vec theta_opt,arma::vec prob_vec)
{
  double res;
  double C;
  double res1;
  double res2;
  double C1;
  double C2;
  List opt_result;
  const double pi =  3.1415926535;
  double log2pi = log(2*pi);
  arma::vec opt_para(2);
  double opt_value;
  int n_nodes = leaf_tree.length();
  arma::vec current_node;
  double n_cell_types;
  int current_cell_type;
  arma::vec N_cell_type;
  arma::vec Y_cell_type;
  arma::mat H_sqrt_inv;
  arma::mat H;
  arma::vec empty;
  arma::vec current_N;
  arma::vec current_Y;
  arma::mat quad_x_mod;
  int n_quad2 = quad_w.size();
  int tmp_size=0;
  int pos = 0;
  int len = 0;
  double mu;
  double theta;
  arma::vec para;
  arma::mat intermediate(2,n_quad2);
  arma::vec smd(n_quad2);
  arma::vec smd1(n_quad2);
  arma::vec smd2(n_quad2);
  arma::vec mu_tilde(n_nodes);
  arma::vec theta_tilde(n_nodes);
  arma::vec smd1_sign(n_quad2);
  arma::vec smd2_sign(n_quad2);
  int res1_sign;
  int res2_sign;
  arma::vec weighted_smd;
  arma::uvec current_node1;
  arma::vec mu_out(n_nodes);
  arma::vec theta_out(n_nodes);
  for(int i=0;i<n_nodes;i++)
  {
    current_node = as<arma::vec>(leaf_tree[i])-1;
    n_cell_types = current_node.size();
    len = 0;
    for(int m =0;m<n_cell_types;m++)
    {
      current_cell_type = current_node[m];
      len += as<arma::vec>(N_all[current_cell_type]).size();
    }
    current_N.resize(len);
    current_Y.resize(len);
    pos = 0;
    for(int j=0;j<n_cell_types;j++)
    {
      current_cell_type = current_node[j];
      N_cell_type = as<arma::vec>(N_all[current_cell_type]);
      Y_cell_type = as<arma::vec>(Y_all[current_cell_type]);
      tmp_size = Y_cell_type.size();
      current_N.subvec(pos,pos+tmp_size-1) = N_cell_type;
      current_Y.subvec(pos,pos+tmp_size-1) = Y_cell_type;
      pos += tmp_size;
    }
    
    H = - f_gj_2_tilde_cpp(mu_opt(i),theta_opt(i),current_Y,current_N,delta,tau,lambda,nu,phi);
    H_sqrt_inv = matrixSquareRootInverse(H);
    intermediate = H_sqrt_inv * quad_x;
    
    opt_para(0) = mu_opt(i);
    opt_para(1) = theta_opt(i);
    intermediate.each_col() += opt_para;
    quad_x_mod = intermediate;
    for(int t = 0; t<n_quad2;t++)
    {
      para = quad_x_mod.col(t);
      mu =para(0);
      theta = para(1);
      smd(t) = f_gj_tilde_cpp(mu,theta,current_Y,current_N,delta,tau,lambda,nu,phi);
      smd1_sign(t) = arma::sign(mu);
      smd2_sign(t) = arma::sign(1);
      smd1(t) = smd(t) +  log(abs(mu));
      smd2(t) = smd(t) +  2*theta;
    }
    opt_value = f_gj_tilde_cpp(mu_opt(i),theta_opt(i),current_Y,current_N,delta,tau,lambda,nu,phi);
    //Rcpp::Rcout << smd1  << std::endl;
    smd1 = quad_w + smd1 - opt_value - log2pi;
    smd2 = quad_w + smd2 - opt_value - log2pi;
    C2 = max(smd2);
    C1 = max(smd1);
    //smd1_sign = arma::sign(smd1);
    //smd2_sign = arma::sign(smd2);
    //smd1 = arma::abs(smd1);
    //smd2 = arma::abs(smd2);
    
    //res2 = sum(smd1_sign %exp(smd2 - C2));
    //res1 =  sum(smd2_sign %exp(smd1 - C1));
    //res1_sign = arma::sign(res1);
    //res2_sign = arma::sign(res2);
    //res1 = abs(res1);
    //res2 = abs(res2);
    
    //res1 = res1_sign * log(res1) + C1;
    //Rcpp::Rcout << res1_sign  << std::endl;
    //res2 = res2_sign * log(res2) + C2;
    res1 = sum(smd1_sign % exp(smd1 - C1));
    res2 = sum(smd2_sign % exp(smd2 - C2));
    res1_sign = arma::sign(res1);
    res2_sign = arma::sign(res2);
    res1 = abs(res1);
    res2 = abs(res2);
    res1 = log(res1) + C1;
    res2 = log(res2) + C2;
    for(int t = 0; t<n_quad2;t++)
    {
      para = quad_x_mod.col(t);
      mu =para(0);
      theta = para(1);
      smd(t) =f_gj_tilde_cpp(mu,theta,current_Y,current_N,delta,tau,lambda,nu,phi);
    }
    smd = quad_w + smd - opt_value - log2pi;
    C = max(smd);
    res = log(sum(exp(smd - C))) + C;
    mu_tilde(i) = res1_sign * exp( res1-res);
    theta_tilde(i) = res2_sign * exp(res2-res);
    
  }
  for(int k=0;k<n_nodes;k++)
  {  
    current_node1 = as<arma::uvec>(parent_tree[k])-1;
    weighted_smd = prob_vec.elem(current_node1) % mu_tilde.elem(current_node1);
    mu_out(k) = sum(weighted_smd);
    weighted_smd = prob_vec.elem(current_node1) % theta_tilde.elem(current_node1);
    theta_out(k) = sum(weighted_smd);
  }
  arma::mat vec_list(n_nodes,2);
  vec_list.col(0) = mu_out;
  vec_list.col(1) = theta_out;
  return vec_list;
}

// [[Rcpp::export]]
arma::mat find_intervals_cpp(const arma::vec& p, const arma::vec& positions, double increment = 75, double threshold = 0.95, int min_count = 6) {
  arma::uword n = p.n_elem;

  arma::uvec sorted_indices = arma::sort_index(p, "descend");
  arma::vec p_sorted = p.elem(sorted_indices);
  arma::vec positions_sorted_p = positions.elem(sorted_indices);
  
  std::vector<std::vector<double>> intervals;
  arma::vec included_positions = arma::zeros<arma::vec>(n);
  arma::vec positions_vector = positions;
  
  double left_limit = positions_vector.min();
  double right_limit = positions_vector.max();
  
  for (arma::uword i = 0; i < n; ++i) {
  
    if ((i + 1) % 2000 == 0) {
      Rcpp::Rcout << "Current iteration: " << i + 1 << ", Current pos: " << positions_sorted_p[i] << std::endl;
    }
    
    double current_p = p_sorted[i];
    if (current_p < threshold) break;  
    
    double current_pos = positions_sorted_p[i];
    arma::uvec idx_in_positions = arma::find(positions_vector == current_pos);
    if (idx_in_positions.is_empty()) continue;
    arma::uword idx = idx_in_positions[0];
    if (included_positions[idx] == 1) continue; 
    
    double start_pos = current_pos;
    double end_pos = current_pos;
    
    bool left_stop = false;
    bool right_stop = false;
    while (true) {
      if (!left_stop) {
        double new_start = start_pos - increment;
        if (new_start < left_limit) {
          new_start = left_limit;
          left_stop = true;
        }
        arma::uvec interval_indices = arma::find(positions_vector >= new_start && positions_vector <= end_pos);
        arma::uvec interval_indices_new = arma::find(positions_vector >= new_start && positions_vector < start_pos);
        if (interval_indices_new.is_empty()) {
          left_stop = true;
        } else {
          double avg_p = arma::mean(p.elem(interval_indices));
          if (avg_p > threshold) {
            start_pos = new_start;
          } else {
            left_stop = true;
          }
        }
      }
      if (!right_stop) {
        double new_end = end_pos + increment;
        if (new_end > right_limit) {
          new_end = right_limit;
          right_stop = true;
        }
        arma::uvec interval_indices = arma::find(positions_vector >= start_pos && positions_vector <= new_end);
        arma::uvec interval_indices_new = arma::find(positions_vector > end_pos && positions_vector <= new_end);
        if (interval_indices_new.is_empty()) {
          right_stop = true;
        } else {
          double avg_p = arma::mean(p.elem(interval_indices));
          if (avg_p > threshold) {
            end_pos = new_end;
          } else {
            right_stop = true;
          }
        }
      }
      if (left_stop && right_stop) break;
    }
    arma::uvec interval_indices = arma::find(positions_vector >= start_pos && positions_vector <= end_pos);
    if (interval_indices.n_elem >= (arma::uword)min_count) {
      intervals.push_back({start_pos, end_pos, current_p});
      included_positions.elem(interval_indices).ones();
    }
  }
  if (!intervals.empty()) {
    arma::mat intervals_mat(intervals.size(), 3);
    for (size_t i = 0; i < intervals.size(); ++i) {
      intervals_mat(i, 0) = intervals[i][0];
      intervals_mat(i, 1) = intervals[i][1];
      intervals_mat(i, 2) = intervals[i][2];
    }
    arma::uvec sort_indices = arma::sort_index(intervals_mat.col(0));
    intervals_mat = intervals_mat.rows(sort_indices);
    arma::mat merged_intervals;
    merged_intervals.set_size(0, 3);
    merged_intervals.insert_rows(0, intervals_mat.row(0));
    for (arma::uword i = 1; i < intervals_mat.n_rows; ++i) {
      arma::rowvec last_interval = merged_intervals.row(merged_intervals.n_rows - 1);
      arma::rowvec current_interval = intervals_mat.row(i);
      if (last_interval[1] >= current_interval[0]) {
        last_interval[1] = std::max(last_interval[1], current_interval[1]);
        last_interval[2] = std::max(last_interval[2], current_interval[2]);
        merged_intervals.row(merged_intervals.n_rows - 1) = last_interval;
      } else {
        merged_intervals.insert_rows(merged_intervals.n_rows, current_interval);
      }
    }
    return merged_intervals;
  } else {
    return arma::mat(0, 3);
  }
}


void writeDataToFileBuffered(const std::string& filename, 
                             const arma::mat& data, 
                             bool append = true, 
                             const std::string& delimiter = ",") {
  std::ios_base::openmode mode = std::ios_base::out;
  if (append) {
    mode |= std::ios_base::app;
  } else {
    mode |= std::ios_base::trunc;
  }
  std::ofstream outfile(filename.c_str(), mode);
  if (!outfile.is_open()) {
    Rcpp::stop("Error: " + filename);
  }
  
  std::ostringstream buffer;
  for (size_t row = 0; row < data.n_rows; ++row) {
    for (size_t col = 0; col < data.n_cols; ++col) {
      buffer << data(row, col);
      if (col != data.n_cols - 1) {
        buffer << delimiter; 
      }
    }
    buffer << "\n"; 
  }
  outfile << buffer.str();
  outfile.close();
}

// [[Rcpp::export]]
void find_pairs(DataFrame df, DataFrame ref, std::string filename, int min_read = 20)
{
  StringVector chr = df["V1"];
  StringVector fragment = df["V3"];
  arma::vec index = df["V2"];
  arma::vec rep = df["V4"];
  StringVector ref_chr = ref["chr"];
  arma::vec ref_pos = ref["pos"];
  arma::vec ref_score = ref["score"];
  arma::vec ref_idx = ref["idx"];
  arma::vec ref_read = ref["read"];
  arma::mat out(500000,6);
  int count = 0;
  int ref_length = ref_chr.size();
  int df_length = fragment.size();
  int i = 0;
  int j = 0;
  int k = 0;
  int current_ref_pos;
  std::string current_fragment;
  int current_fragment_length;
  int CpG_1;
  int CpG_2;
  double mu_1;
  double mu_2;
  int distance;
  int replicate;
  int current_index;
  arma::rowvec added_row(6);
  while(true)
  {
    if(i % 200000 == 0)
    {
      Rcpp::Rcout << i  << std::endl;
    }
    if(i == df_length - 1)
    { 
      arma::uvec rowIndices = arma::regspace<arma::uvec>(0, count-1);
      writeDataToFileBuffered(filename, out.rows(rowIndices));
      break;
    }
    current_fragment = fragment[i];
    current_index = index[i] - 1;
    current_fragment_length = current_fragment.size();
    if(current_fragment_length == 1)
    {
      i += 1;
      continue;
    }
    j = 0;
    k = 1;
    while(true)
    {
      if((j == current_fragment_length - 1) | (ref_read[j + current_index] < min_read))
      {
        i += 1;
        break;
      }
      current_ref_pos = ref_pos[j + current_index];
      while(true)
      {
        if(k == current_fragment_length)
        {
          break;
        }
        if((abs(ref_pos[k + current_index] - current_ref_pos)<=200) & (ref_read[k + current_index] >= min_read))
        {
          CpG_1 = (current_fragment[j] == 'C') ? 1 : 0;
          CpG_2 = (current_fragment[k] == 'C') ? 1 : 0;
          mu_1 = R::qnorm(ref_score[current_index + j], 0.0, 1.0, 1, 0);
          mu_2 = R::qnorm(ref_score[current_index + k], 0.0, 1.0, 1, 0);
          //mu_1 = ref_score[current_index + j];
          //mu_2 = ref_score[current_index + k];
          replicate = rep[i];
          distance = ref_pos[k + current_index] - current_ref_pos;
          added_row[0] = CpG_1;
          added_row[1] = CpG_2;
          added_row[2] = mu_1;
          added_row[3] = mu_2;
          added_row[4] = distance;
          added_row[5] = replicate;
          out.row(count) = added_row;
          count += 1;
          if(count == 500000)
          {
            writeDataToFileBuffered(filename, out);
            count = 0;
          }
          k += 1;
        }else{
          break;
        }
          
      }
      j += 1;
      k = j + 1;
    }
  }
  
}

// [[Rcpp::export]]
void find_pairs_subj(List df, std::string filename)
{
  arma::vec pos = df["pos"];
  arma::vec mean = df["mean"];
  arma::vec sd = df["sd"];
  arma::mat mu = df["mu"];
  double m = static_cast<int>(mu.n_cols);
  m -= 1;
  int n = pos.size();
  int i = 0;
  int j = 1;
  int current_pos;
  arma::mat out(500000,2);
  int count = 0;
  int distance;
  double corr;
  arma::rowvec added_row(2);
  while(true)
  {
    if(i % 200000 == 0)
    {
      Rcpp::Rcout << i  << std::endl;
    }
    if(i == n - 1)
    {
      arma::uvec rowIndices = arma::regspace<arma::uvec>(0, count-1);
      writeDataToFileBuffered(filename, out.rows(rowIndices));
      break;
    }
    current_pos = pos[i];
    if(abs(pos[j] - current_pos )<=200 & j <= n - 1)
    {
      distance = pos[j] - current_pos;
      corr = dot((mu.row(i)-mean[i])/sd[i], (mu.row(j)-mean[j])/sd[j])/m;
      added_row[0] = distance;
      added_row[1] = corr;
      out.row(count) = added_row;
      count += 1;
      if(count == 500000)
      {
        writeDataToFileBuffered(filename, out);
        count = 0;
      }
      j += 1;
    }
    else{
      i += 1;
      j = i + 1;
    }
  }
}

double pmvn_rcpp_arma(arma::vec lower, arma::vec upper, arma::vec mean, arma::mat sigma) {
  NumericVector lower_r = NumericVector(lower.begin(),lower.end());
  NumericVector upper_r = NumericVector(upper.begin(),upper.end());
  NumericVector mean_r = NumericVector(mean.begin(),mean.end());
  NumericMatrix sigma_r = as<NumericMatrix>(wrap(sigma));
  Environment mvtnorm = Environment::namespace_env("mvtnorm");
  Function pmvnorm = mvtnorm["pmvnorm"];
  NumericVector result = pmvnorm(_["lower"] = lower_r,
                                 _["upper"] = upper_r,
                                 _["mean"] = mean_r,
                                 _["sigma"] = sigma_r,
                                 _["abseps"] = 1e-3);
  double cdf = result[0];
  return cdf;
}

arma::mat calculate_sigma_cpp(arma::mat distance,
                              double rho1 = 0.5058050, double omega1 = 0.5209102,
                              double theta1 = 163.5032626, double theta2 = -28.5323214)
{
  int n = distance.n_rows;
  arma::mat sigma_mat(n,n);
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<=i;j++)
    {
      if(i == j)
      {
        sigma_mat(i,j) = 1;
      }
      else{
        sigma_mat(i,j) = rho1*(omega1*exp(-(distance(i,j)/theta1))+(1-omega1)*exp(-pow(distance(i,j)/theta2,2)));
        sigma_mat(j,i) = sigma_mat(i,j);
      }
    }
  }
  return sigma_mat;
}


// arma::mat sigma_mat()

// [[Rcpp::export]]
List find_overlap_fragments(arma::ivec frag_start, std::vector<std::string> frag_seq,
                            arma::ivec frag_count,
                            arma::ivec ref_start, arma::ivec ref_end) {
  int n_frag = frag_start.n_elem;
  int n_ref = ref_start.n_elem;
  
  std::vector<int> res_start_idx;
  std::vector<int> res_end_idx;
  std::vector<std::string> res_seq;
  std::vector<int> res_count; 
  std::vector<int> DMRs_idx; 
  int ref_idx = 0;
  for(int i = 0; i < n_frag; i++) {
    int frag_s = frag_start[i];
    std::string seq = frag_seq[i];
    int count = frag_count[i];
    int frag_len = seq.size();
    int frag_e = frag_s + frag_len - 1;
    while(ref_idx < n_ref && ref_end[ref_idx] < frag_s) {
      ref_idx++;
    }

    int tmp_ref_idx = ref_idx;

    while(tmp_ref_idx < n_ref && ref_start[tmp_ref_idx] <= frag_e) {
      int overlap_start = std::max(frag_s, ref_start[tmp_ref_idx]);
      int overlap_end = std::min(frag_e, ref_end[tmp_ref_idx]);
      int seq_start = overlap_start - frag_s;
      int seq_end = overlap_end - frag_s;
      std::string overlap_seq = seq.substr(seq_start, seq_end - seq_start + 1);
      res_start_idx.push_back(overlap_start);
      res_end_idx.push_back(overlap_end);
      res_seq.push_back(overlap_seq);
      res_count.push_back(count); 
      DMRs_idx.push_back(tmp_ref_idx);
      tmp_ref_idx++;
    } 
  }
  
  return List::create(
    _["start_idx"] = wrap(res_start_idx),
    _["end_idx"] = wrap(res_end_idx),
    _["seq"] = res_seq,
    _["count"] = wrap(res_count),
    _["DMR_idx"] = wrap(DMRs_idx)
  );
}

// [[Rcpp::export]]
arma::uvec findIndicesInRange(const arma::vec& sorted_vec, double start, double end) {
  auto lower = std::lower_bound(sorted_vec.begin(), sorted_vec.end(), start);
  auto upper = std::upper_bound(sorted_vec.begin(), sorted_vec.end(), end);
  int start_idx = lower - sorted_vec.begin();
  int end_idx = upper - sorted_vec.begin() - 1;
  if (start_idx > end_idx) {
    return arma::uvec(); 
  }
  arma::uvec indices(end_idx - start_idx + 1);
  for (int i = start_idx; i <= end_idx; ++i) {
    indices[i - start_idx] = i;
  }
  
  return indices;
} 

// [[Rcpp::export]]
arma::vec calculate_1_fragment_lik_cpp(std::string fragment, arma::vec pos, arma::mat mu,
                                    arma::mat sigma2,arma::mat corr_para_mat)
{
  int n = fragment.size();
  arma::vec lower(n);
  arma::vec upper(n);
  for(int i = 0;i<n;i++)
  { 
    if(fragment[i] == 'C')
    {
      lower(i) = 0.0;
      upper(i) = 999.0;
    }else if(fragment[i] == 'T'){
      lower(i) = -999.0;
      upper(i) = 0.0;
    }else{
      lower(i) = -999.0;
      upper(i) = 999.0;
    }
  } 
  arma::mat sigma_mat(n,n);
  arma::mat dist_mat(n,n);
  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      dist_mat(i, j) = std::abs(pos(i) - pos(j));
      dist_mat(j, i) = dist_mat(i, j);
    }
  }
  
  int m = mu.n_cols;
  arma::vec likelihood(m);
  for(int i = 0;i<m;i++)
  {
    arma::vec tmp_corr_para = corr_para_mat.col(i);
    arma::mat corr_mat = calculate_sigma_cpp(dist_mat,tmp_corr_para(0),tmp_corr_para(1),
                                             tmp_corr_para(2),tmp_corr_para(3));
    arma::vec tmp_sigma2 = sigma2.col(i);
    arma::vec tmp_mu = mu.col(i);
    sigma_mat = corr_mat;
    //for(int j = 0;j <n;j++)
    //{
    //  sigma_mat(j,j) = corr_mat(j,j) + tmp_sigma2(j);
    //}
    likelihood(i) = pmvn_rcpp_arma(lower,upper,tmp_mu,sigma_mat);
    
  }

  return(likelihood);
}

// [[Rcpp::export]]
List calculate_fragment_lik_cpp(DataFrame overlap_fragments,
                                     DataFrame reference, arma::vec index,
                                     arma::mat mu_mat, arma::mat sigma2_mat,List corr_mat)
{
  arma::vec position = reference["pos"];
  arma::vec start_idx = overlap_fragments["start_idx"];
  arma::vec end_idx = overlap_fragments["end_idx"];
  arma::uvec DMR_idx = overlap_fragments["DMR_idx"];
  int m = mu_mat.n_cols;
  int n = end_idx.size();
  std::vector<std::string> all_fragments = overlap_fragments["seq"];
  arma::vec count = overlap_fragments["count"];
  arma::mat out(m,n);
  for(int i =0;i<n;i++)
  {
    std::string fragment = all_fragments[i];
//    Rcpp::Rcout << fragment  << std::endl;
    arma::uvec current_idx = findIndicesInRange(index,start_idx(i),end_idx(i));
    arma::mat mu = mu_mat.rows(current_idx);
    arma::mat sigma2 = sigma2_mat.rows(current_idx);
    arma::vec pos = position(current_idx);
    arma::mat corr_para_mat = corr_mat(DMR_idx(i));
    out.col(i) = calculate_1_fragment_lik_cpp(fragment,pos,mu,sigma2,corr_para_mat);
    if(i % 1000 == 0)
    {
      Rcpp::Rcout << i  << std::endl;
    }
    //Rcpp::Rcout << out.col(i)   << std::endl;
  }
  return List::create(
    _["likelihood"] = out,
    _["replicates"] = count
  );
}

// [[Rcpp::export]]
NumericMatrix find_pairs_cpp(DataFrame df, DataFrame ref, int min_read = 20) {
  std::vector<std::string> fragment = df["seq"];
  arma::vec index = df["start_idx"];
  arma::vec rep = df["count"];
  
  StringVector ref_chr = ref["chr"];
  arma::vec ref_pos = ref["pos"];
  arma::vec ref_score = ref["score"];
  arma::vec ref_idx = ref["idx"];
  arma::vec ref_read = ref["read"];
  
  int ref_length = ref_chr.size();
  int df_length = fragment.size();
  bool all_single_length = true;
  for (int i = 0; i < df_length; i++) {
    if ((int)fragment[i].size() > 1) {
      all_single_length = false;
      break;
    }
  }

  if (all_single_length) {
    NumericMatrix empty_result(0, 6);
    return empty_result;
  }
  
  std::unordered_map<int, int> idx_map;
  idx_map.reserve(ref_length);
  for (int r = 0; r < ref_length; r++) {
    idx_map[(int)ref_idx[r]] = r;
  }
  
  arma::mat out(500000, 6);
  int count = 0;
  
  for (int i = 0; i < df_length; i++) {
    std::string current_fragment = fragment[i];
    int current_fragment_length = (int)current_fragment.size();
    
    if (current_fragment_length <= 1) {
      continue;  
    }
    
    int idx_val = (int)index[i];
    auto it = idx_map.find(idx_val);
    if (it == idx_map.end()) {
      continue; 
    }
    
    int current_index = it->second;
    int j = 0, k = 1;
    while (j < current_fragment_length) {
      if ((j + current_index) >= ref_length || (j + current_index) < 0) {
        break;  
      }
      
      if (ref_read[j + current_index] < min_read) {
        j += 1;
        k = j + 1;
        continue; 
      }
      
      int current_ref_pos = (int)ref_pos[j + current_index];
      while (k < current_fragment_length) {
        if ((k + current_index) >= ref_length || (k + current_index) < 0) {
          break;  
        }
        
        if (current_fragment[j] == '.' || current_fragment[k] == '.') {
          k += 1;
          continue;  
        }
        
        if (std::abs((int)ref_pos[k + current_index] - current_ref_pos) <= 200 &&
            ref_read[k + current_index] >= min_read) {
          int CpG_1 = (current_fragment[j] == 'C') ? 1 : 0;
          int CpG_2 = (current_fragment[k] == 'C') ? 1 : 0;
          double mu_1 = ref_score[current_index + j];
          double mu_2 = ref_score[current_index + k];
          int distance = (int)ref_pos[k + current_index] - current_ref_pos;
          int replicate = (int)rep[i];
          
          if (count >= (int)out.n_rows) {
            arma::mat temp(out.n_rows * 2, 6);
            temp.rows(0, count - 1) = out.rows(0, count - 1);
            out = temp;
          }
          
          out(count, 0) = CpG_1;
          out(count, 1) = CpG_2;
          out(count, 2) = mu_1;
          out(count, 3) = mu_2;
          out(count, 4) = distance;
          out(count, 5) = replicate;
          count += 1;
        }
        k += 1;
      }
      j += 1;
      k = j + 1;
    }
  }
  
  arma::mat final_out = out.rows(0, count - 1);
  return wrap(final_out);
}
