#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List fx(const arma::mat& EM,
              const arma::cube& TM){
  //arma::mat EM = Rcpp::as<arma::mat>(em);
  //arma::cube TM = Rcpp::as<arma::cube>(tm);

  short int nm = EM.n_rows, ns = EM.n_cols;
  arma::mat FM(ns, nm + 1);
  arma::rowvec e, f;
  arma::vec s(nm);

  FM.fill(1.0/ns);

  for(int i = 0; i < nm; ++i){

    f = arma::trans(FM.col(i))*TM.slice(i);
    e = EM.row(i);

    for(int j = 0; j < ns; ++j){
      FM(j,i+1) = f(j)*e(j);
    }

    s(i) = sum(FM.col(i+1));
    FM.col(i+1) = FM.col(i+1)/s(i); // normalise using scaling

  }

  return Rcpp::List::create(Rcpp::Named("FM") = FM,
                            Rcpp::Named("s") = s);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bx(const arma::vec& s,
                       const arma::mat& EM,
                       const arma::cube& TM){
  // arma::mat EM = Rcpp::as<arma::mat>(em);
  //arma::cube TM = Rcpp::as<arma::cube>(tm);
  //arma::vec SV = Rcpp::as<arma::vec>(sv);

  short int nm = EM.n_rows, ns = EM.n_cols;
  arma::mat B(nm, ns);
  arma::rowvec b, e;
  arma::vec be(ns);

  B.fill(1.0/ns);

  for(int i = nm - 1; i > 0; --i){
    b = B.row(i);
    e = EM.row(i);

    for(int j = 0; j < ns; ++j){
      be(j) = b(j)*e(j);
    }

    B.row(i-1) = (arma::trans(be)*TM.slice(i))/s(i);
  }

  return Rcpp::wrap(arma::trans(B));
}
