#rm(list=ls(all=T))
library(inline)

txt<-'
#include<Eigen/Dense>
using namespace Rcpp;
using namespace Eigen;

MatrixXd m1=as<MatrixXd> (xs);
ColPivHouseholderQR<MatrixXd> dec(m1);
MatrixXd q1=dec.householderQ();
MatrixXd p1=dec.colsPermutation();
//MatrixXd r1=dec.matrixR().triangularView<Upper>()*p1.inverse();
MatrixXd r1=dec.matrixR().triangularView<Upper>();

return List::create(Named("Q")=q1, Named("R")=r1, Named("P")=p1);
'

qr1<-cxxfunction(signature(xs="numeric"),body=txt,plugin="RcppEigen")

#m1<-matrix(rnorm(9),nr=3)
#q1<-qr(m1)
#q2<-qr1(m1)
