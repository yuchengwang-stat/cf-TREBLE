// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using arma::mat;
using arma::vec;
using arma::rowvec;
using arma::accu;

// [[Rcpp::export]]
Rcpp::List my_celfie_rcpp(const vec& X,          // length G: methylated reads
                          const vec& D_X,        // length G: total reads
                          const mat& W,          // G × CT: reference methylated reads
                          const mat& D_W,        // G × CT: reference total reads
                          int   n_iter = 10000,  // max iterations (default)
                          double thres = 1e-5)   // convergence threshold 
{
  const std::size_t G  = X.n_elem;
  const std::size_t CT = W.n_cols;
  
  // ---------------- Initialization ----------------
  vec mu(CT, arma::fill::ones);
  mu /= static_cast<double>(CT);                 
  
  mat phi = W / D_W;                             
  mat phi_unmeth = 1.0 - phi;
  
  // workspaces
  mat U1(G, CT);
  mat U0(G, CT);
  vec U1_sum(G);
  vec U0_sum(G);
  
  int n_converge = 0;
  
  for (int iter = 0; iter < n_iter; ++iter) {
    ++n_converge;

    rowvec mu_row = mu.t();
    U1 = phi;            U1.each_row() %= mu_row;   
    U0 = phi_unmeth;     U0.each_row() %= mu_row;   
    
    U1_sum = arma::sum(U1, 1);
    U0_sum = arma::sum(U0, 1);
    
    vec tmp1 = X / U1_sum;
    vec tmp0 = (D_X - X) / U0_sum;
    
    U1.each_col() %= tmp1;
    U0.each_col() %= tmp0;
    mat U = U1 + U0;
    
    vec mu_old = mu;
    mu = arma::sum(U, 0).t();
    mu /= accu(mu);

    phi = (U1 + W) / (U + D_W);
    phi_unmeth = 1.0 - phi;

    if (accu(arma::abs(mu - mu_old)) < thres) {
      return Rcpp::List::create(Rcpp::Named("cellular_frac") = mu,
                                Rcpp::Named("iterations")    = n_converge);
    }
  }
  
  // reached max iterations
  return Rcpp::List::create(Rcpp::Named("cellular_frac") = mu,
                            Rcpp::Named("iterations")    = n_converge);
}

/*
 Usage in R:
 Rcpp::sourceCpp("my_celfie_rcpp.cpp")
 res <- my_celfie_rcpp(X, D_X, W, D_W)
 res$cellular_frac
 res$iterations
 */
