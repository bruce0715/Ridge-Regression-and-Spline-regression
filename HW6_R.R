1#############################################################
## Author: Bruce
## Date : Nov 12, 2017
## Description: This script implements ridge regression as 
## well as piecewise linear spline regression.
#############################################################
 
#############################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. You can add examples at the
## end of the script (in the "Optional examples" section) to 
## double-check your work, but MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not use the function "setwd" anywhere
## in your code. If you do, I will be unable to grade your 
## work since R will attempt to change my working directory
## to one that does not exist.
#############################################################

## Source your Rcpp file (put in the name of your 
## Rcpp file)
library(Rcpp)

#Please put the Cpp file in working directory to source it
sourceCpp('HW6_CPP.cpp')

##################################
## Function 1: QR decomposition ##
##################################

myQR <- function(A){
  
  ## Perform QR decomposition on the matrix A
  ## Input: 
  ## A, an n x m matrix
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################  
  n=dim(A)[1]
  m=dim(A)[2]
  R=A
  Q=diag(rep(1,n))
  k=ifelse(m==n,m-1,m)
  
  for(i in 1:k){
    H=diag(rep(1,n))
    a=R[i:nrow(R),i]
    xs=-sign(a[1])*norm(a,"2")
    a2=rep(0,length(a))
    a2[1]=xs
    u=(a-a2)/norm(a-a2,"2")
    H2=diag(rep(1,length(a)))
    H2=H2-(2*u%*%t(u))
    H[i:nrow(H),i:ncol(H)]=H2
    
    R=H%*%R
    Q=H%*%Q
  }
  
  
  
  ## Function should output a list with Q.transpose and R
  ## Q is an orthogonal n x n matrix
  ## R is an upper triangular n x m matrix
  ## Q and R satisfy the equation: A = Q %*% R
  return(list("Q" = t(Q), "R" = R))
  
}


#################################
## Function 2: Sweep operation ##
#################################

mySweep <- function(A, m){
  
  # Perform a SWEEP operation on A with the pivot element A[m,m].
  # 
  # A: a square matrix.
  # m: the pivot element is A[m, m].
  # Returns a swept matrix.
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################
  n<-dim(A)[1]
  B=A
  
  for (k in 1:m){
    for (i in 1:n){
      for (j in 1:n){
        if (i!=k&j!=k){
          B[i,j]=B[i,j]-B[i,k]*B[k,j]/B[k,k]
        }
      }
    }
    for (i in 1:n){
      if(i!=k){
        B[i,k]=B[i,k]/B[k,k]
      }
    }
    for (j in 1:n){
      if(j!=k){
        B[k,j]=B[k,j]/B[k,k]
      }
    }
    B[k,k]=-1/B[k,k]
  }
  
  ## The function outputs the matrix B
  A=B
  return(A)
  
}


##################################
## Function 3: Ridge regression by Sweep ##
##################################

myRidge <- function(X, Y, lambda){
  
  # Perform ridge regression of Y on X.
  # 
  # X: an n x p matrix of explanatory variables.
  # Y: an n vector of dependent variables. Y can also be a 
  # matrix, as long as the function works.
  # lambda: regularization parameter (lambda >= 0)
  # Returns beta, the ridge regression solution.
  
  ##################################
  ## FILL IN THIS SECTION OF CODE ##
  ##################################
  n=dim(X)[1]
  p=dim(X)[2]
  Z=cbind(rep(1,n),X,Y)
  A=t(Z)%*%Z
  D=diag(rep(lambda,(p+2)))
  D[1,1]=0 #don't add for beta_0
  D[(p+2),(p+2)]=0 #don't add for y
  A=A+D
  S=mySweepC(A,p+1)
  #S=mySweep(A,(p+1))
  
  #if use mySweep or mySweepC
  beta_ridge=S[1:(p+1),(p+2)]
  
  ## Function should output the vector beta_ridge, the 
  ## solution to the ridge regression problem. beta_ridge
  ## should have p + 1 elements.
  return(beta_ridge)
}

##################################
## Function 3.1: Ridge regression QR ##
##################################

# myRidgeQR <- function(X, Y, lambda){
#   
#   # Perform ridge regression of Y on X.
#   # 
#   # X: an n x p matrix of explanatory variables.
#   # Y: an n vector of dependent variables. Y can also be a 
#   # matrix, as long as the function works.
#   # lambda: regularization parameter (lambda >= 0)
#   # Returns beta, the ridge regression solution.
#   ##################################
#   ## FILL IN THIS SECTION OF CODE ##
#   ##################################
#   n=dim(X)[1] 
#   p=dim(X)[2]
#   
#   X2=rbind(X,diag(rep(sqrt(lambda),p)))
#   Y2=matrix(c(Y,rep(0,p)),(n+p),1)
#   
#   Z=cbind(c(rep(1,n),rep(0,p)),X2,Y2)
#   #R=myQR(Z)$R
#   R=myQRC(Z)$R
#   R1=R[1:(p+1),1:(p+1)]
#   Y1=R[1:(p+1),p+2]
#   beta_ridge=solve(R1,Y1)
#   
#   ## Function should output the vector beta_ridge, the 
#   ## solution to the ridge regression problem. beta_ridge
#   ## should have p + 1 elements.
#   return(beta_ridge)
# }


####################################################
## Function 4: Piecewise linear spline regression ##
####################################################


mySpline <- function(x, Y, lambda, p = 100){
  
  # Perform spline regression of Y on X.
  # 
  # x: An n x 1 vector or matrix of explanatory variables.
  # Y: An n x 1 vector of dependent variables. Y can also be a 
  # matrix, as long as the function works.
  # lambda: regularization parameter (lambda >= 0)
  # p: Number of cuts to make to the x-axis.
  
  ##################################
  ## FILL IN THIS SECTION OF CODE ##
  ##################################
  n=dim(x)[1]
  
  X=matrix(x,nrow=n)
  for (k in (1:(p-1))/p)
    X=cbind(X,(x>k)*(x-k))
  beta_spline=myRidge(X,Y,lambda)
  #beta_spline=myRidgeQR(X,Y,lambda)
  
  Yhat=cbind(rep(1,n),X)%*%beta_spline
  
  
  ## Function should a list containing two elements:
  ## The first element of the list is the spline regression
  ## beta vector, which should be p + 1 dimensional (here, 
  ## p is the number of cuts we made to the x-axis).
  ## The second element is y.hat, the predicted Y values
  ## using the spline regression beta vector. This 
  ## can be a numeric vector or matrix.
  output <- list(beta_spline = beta_spline, predicted_y = Yhat)
  return(output)
}

#########test performance############
 n=20
 p=100
 sigma=0.1
 lambda=1.
 set.seed(90)
 x=runif(n)
 x=sort(x)
 set.seed(90)
 Y=x^2+rnorm(n)*sigma
 X=matrix(x,nrow=n)

 ptm <- Sys.time()

output=mySpline(X,Y,lambda=lambda,p=p)

#Sys.time()-ptm














# 
# 
# ######plot graph###########
# lambda=0
# temp_output=mySpline(X,Y,lambda,p=p)
# plot(X,Y,main='lambda=0',xlab='X values',ylab='Y values')
# abline(myRidge(X,Y,lambda)[1],myRidge(X,Y,lambda)[2],lwd=2)
# #par(new=TRUE)
# lines(X,temp_output$predicted_y,type='l',col='red',lwd=2)
# legend('topleft',legend = c('Ridge','Spline'),col=c('black','red'),lty=1,cex=0.8)
# 
# lambda=0.01
# temp_output=mySpline(X,Y,lambda,p=p)
# plot(X,Y,main='lambda=0.01',xlab='X values',ylab='Y values')
# abline(myRidge(X,Y,lambda)[1],myRidge(X,Y,lambda)[2],lwd=2)
# #par(new=TRUE)
# lines(X,temp_output$predicted_y,type='l',col='red',lwd=2)
# legend('topleft',legend = c('Ridge','Spline'),col=c('black','red'),lty=1,cex=0.8)
# 
# lambda=1
# temp_output=mySpline(X,Y,lambda,p=p)
# plot(X,Y,main='lambda=1',xlab='X values',ylab='Y values')
# abline(myRidge(X,Y,lambda)[1],myRidge(X,Y,lambda)[2],lwd=2)
# #par(new=TRUE)
# lines(X,temp_output$predicted_y,type='l',col='red',lwd=2)
# legend('topleft',legend = c('Ridge','Spline'),col=c('black','red'),lty=1,cex=0.8)
# 
# lambda=5
# temp_output=mySpline(X,Y,lambda,p=p)
# plot(X,Y,main='lambda=5',xlab='X values',ylab='Y values')
# abline(myRidge(X,Y,lambda)[1],myRidge(X,Y,lambda)[2],lwd=2)
# #par(new=TRUE)
# lines(X,temp_output$predicted_y,type='l',col='red',lwd=2)
# legend('topleft',legend = c('Ridge','Spline'),col=c('black','red'),lty=1,cex=0.8)
# 
# 
# #######plot error##############
# #create test data
# n_t=20
# p_t=100
# sigma_t=0.1
# lambda_t=1.
# set.seed(1)
# x_t=runif(n_t)
# x_t=sort(x_t)
# set.seed(1)
# Y_t=x_t^2+rnorm(n_t)*sigma_t
# X_t=matrix(x_t,nrow=n_t)
# 
# 
# ##for ridge error###
# 
# #training error
# train_err=c()
# for(lambda in seq(0,5,0.005)){
#   beta=myRidge(X,Y,lambda)
#   err=sum((Y-(beta[2]*X+beta[1]))^2)/n
#   train_err=c(train_err,err)
# }
# 
# #testing error
# test_err=c()
# for(lambda in seq(0,5,0.005)){
#   beta2=myRidge(X,Y,lambda)
#   err2=sum((Y_t-(beta2[2]*X_t+beta2[1]))^2)/n_t
#   test_err=c(test_err,err2)
# }
# 
# plot(seq(0,5,0.005),test_err,type='l',ylim=c(0,0.075),xlab='Lambda Values',ylab='Avg. Errors',main="Training VS Test Error (Ridge Regression)",col='green',lwd=2)
# lines(seq(0,5,0.005),train_err,type='l',col='red',lwd=2)
# legend('bottomright',legend = c('Testing','Training'),col=c('green','red'),lty=1,cex=0.8)
# 
# 
# ####for spline error######
# #training error
# p=100
# 
# s_train_err=c()
# for(lambda in seq(0,1,0.02)){
#   s_output=mySpline(X,Y,lambda,p)
#   s_err=sum((Y-s_output$predicted_y)^2)/n
#   s_train_err=c(s_train_err,s_err)
# }
# 
# #testing error
# s_test_err=c()
# X_t2=X_t
# for (k in (1:(p-1))/p){
#   X_t2=cbind(X_t2,(X_t>k)*(X_t-k))
# }
# 
# 
# for(lambda in seq(0,1,0.02)){
#   
#   out_put_t=mySpline(X,Y,lambda,p)
#   s_beta_t=out_put_t$beta_spline
#   
#   Yhat_t=cbind(rep(1,n_t),X_t2)%*%matrix(s_beta_t,length(s_beta_t),1)
#   s_err=sum((Y_t-Yhat_t)^2)/n_t
#   s_test_err=c(s_test_err,s_err)
# }
# 
# 
# 
# plot(seq(0.02,1,0.02),s_test_err[-1],type='l',ylim=c(0.005,0.01),xlab='Lambda Values',ylab='Avg.Errors',main="Training VS Test Error (Spline Regression)",col='green',lwd=2)
# lines(seq(0.02,1,0.02),s_train_err[-1],type='l',ylim=c(0,0.007),col='red',lwd=2)
# legend('bottomright',legend = c('Testing','Training'),col=c('green','red'),lty=1,cex=0.8)
