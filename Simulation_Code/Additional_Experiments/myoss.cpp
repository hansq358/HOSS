#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <queue>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
arma::vec lev(const arma::mat& X) {
  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, X, "left");
  arma::vec row_crossprods(U.n_rows);
  
  for (arma::uword i = 0; i < U.n_rows; ++i) {
    arma::rowvec row = U.row(i);
    row_crossprods(i) = arma::dot(row, row);
  }
  
  arma::vec result = row_crossprods / X.n_cols;
  
  return result;
}



// [[Rcpp::export]]
arma::vec bottom_k(arma::vec x, unsigned int k) {
  arma::vec x2 = x; // save a copy of x
  arma::vec ind(k); // save the indexes of the smallest k numbers
  std::nth_element(x.begin(), x.begin() + k - 1, x.end()); // std::greater<double>());
  for(int ii=0, i=0; i<int(x.n_elem) && ii<int(k); i++){
    if(x2[i] <= x[k-1])  ind[ii++] = i;  // +1 for R
  }
  return ind;
}

// [[Rcpp::export]]
arma::vec top_k(arma::vec x, unsigned int k) {
  return bottom_k(-x,k);
}
  
// [[Rcpp::export]]
arma::vec Dscr_cpp(arma::mat X, arma::vec xa, arma::mat y, double ya, double tPow) {
  // set tPow=2 for D2 and tPow=4 for D4
  // y is a row vector
  int n=X.n_rows;
  int p=X.n_cols;
  arma::vec B = zeros<vec>(n);
  for(int i=0; i<n; i++){
    B(i) = pow(accu(X.row(i)==y)+p-xa(i)/2-ya/2,tPow); // current used
    //B(i) = pow(p-accu(abs(X.row(i)-y))/2,tPow); // not work
  }
  return B;
}



// [[Rcpp::export]]
arma::uvec OAJ2_cpp(arma::mat x, int k, double tPow=2){
  int n=x.n_rows;
  arma::uvec candi=linspace<uvec>(1,n,n);
  arma::uvec ind=linspace<uvec>(1,k,k);
  arma::vec L=sum(pow(x,2),1);
  arma::vec xa=L;
  uword mm=L.index_max();
  ind(0)=candi(mm);
  candi.shed_row(mm);
  L.shed_row(mm);
  
  arma::mat sx=sign(x);
  double r=log(n/k)/log(k);
  for(int i=1; i<k; i++){
    if(i==1)
      L=Dscr_cpp(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);
    else
      L=L+Dscr_cpp(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);
    
    
    mm=L.index_min();
    ind(i)=candi(mm);
    candi.shed_row(mm);
    L.shed_row(mm);
    
    int nc=floor(n/pow(i,r));
    //Rcout << ind(i) << std::endl;
    //double nc=n/pow(i,r)/L.n_elem;
    if((i>1) & (L.n_elem>double(nc))){
      //arma::uvec tt=as<arma::uvec>(bottom_k(L,nc));
      //Rcout << bottom_k(L,nc) << std::endl;
      arma::uvec tt=arma::conv_to<arma::uvec>::from(bottom_k(L,nc));
      L=L.elem(tt);
      candi=candi.elem(tt);
    }
  }
  return ind;
}



















