/*
####################################################
## Author: Bruce
## Date : Nov 12,2017
## Description: This script implements QR and Sweep
####################################################
 
###########################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not change your working directory
## anywhere inside of your code. If you do, I will be unable 
## to grade your work since R will attempt to change my 
## working directory to one that does not exist.
###########################################################
 
 */ 

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


/* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
Sign function for later use 
~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
double signC(double d){
  return d<0?-1:d>0? 1:0;
}

int justi(int m,int n){
  if(m==n){
    m=m-1;
  }
  return (m);
}

mat subsitute(mat m1,mat m2, int i, int n){
  m1.submat(i,i,n-1,n-1)=m2;
  return m1;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Problem 1: QR decomposition 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */  

// [[Rcpp::export()]]
List myQRC(const mat A){ 
  
  /*
  Perform QR decomposition on the matrix A
  Input: 
  A, an n x m matrix (mat)
  
#############################################
## FILL IN THE BODY OF THIS FUNCTION BELOW ##
#############################################
  
  */
  mat R=A;
  int n = R.n_rows;
  int m = R.n_cols;  
  mat Q = mat(n,n,fill::eye);
  int k= justi(m,n);
  
  for(int i=0;i<k;i++){
    mat H=mat(n,n,fill::eye);
    vec a=R(span(i,n-1),i);
    double xs=signC(a[0])*norm(a,2)*(-1);
    vec a2=zeros<vec>(a.n_elem);
    a2[0]=xs;
    vec u=(a-a2)/norm((a-a2),2);
    mat H2=mat(a.n_elem,a.n_elem,fill::eye);
    H2=H2-2*(u*u.t());
    H=subsitute(H,H2,i,n);
    R=H*R;
    Q=H*Q;
  } 
  // Function should output a List 'output', with 
  // Q.transpose and R
  // Q is an orthogonal n x n matrix
  // R is an upper triangular n x m matrix
  // Q and R satisfy the equation: A = Q %*% R
  List output;
  output["Q"] = Q.t();
  output["R"] = R;
  return(output);
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
Problem 2: Sweep operator 
~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
mat mySweepC(const mat A, int m){
  /*
  Perform a SWEEP operation on A with the pivot element A[m,m].
  
  A: a square matrix (mat).
  m: the pivot element is A[m, m]. 
  Returns a swept matrix B (which is m by m).
  
  Note the "const" in front of mat A; this is so you
  don't accidentally change A inside your code.
  
#############################################
## FILL IN THE BODY OF THIS FUNCTION BELOW ##
#############################################
  
*/ 
  
  mat B = A;
  int n = B.n_rows;
  
  for(int k=0;k<m;k++){
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        if((i!=k) & (j!=k)){
          B(i,j)=B(i,j)-B(i,k)*B(k,j)/B(k,k);
        }
      }
    }
    for(int i=0;i<n;i++){
      if(i!=k){
        B(i,k)=B(i,k)/B(k,k);
      }
    }
    for(int j=0;j<n;j++){
      if(j!=k){
        B(k,j)=B(k,j)/B(k,k);
      }
    }
    B(k,k)=-1/B(k,k);
  }
  
// Return swept matrix B
  return(B);
}
