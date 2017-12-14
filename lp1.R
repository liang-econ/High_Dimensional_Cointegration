rm(list=ls(all=T))
library(glmnet)
library(Matrix)
library(lsgl)
library(LDRTools)

source("qr.R")
####################################################
#kk=3

rank0<-1
dm0<-2

pi0<-matrix(c(-1,1,-0.5,0.5),nr=2)
pi1<-pi0+diag(dm0)
print(abs(eigen(pi1)$val))

bb1<-diag(c(0.4,0.4))
bb2<-matrix(0,nr=2,nc=2)
bb3<-diag(c(0.4,0.4))

mc1<-cbind(bb3,bb2-bb3,bb1-bb2,-(pi1+bb1))*-1
tmc1<-diag(3*dm0)
tmc2<-cbind(matrix(0,nr=3*dm0,nc=dm0),tmc1)
mc2<-rbind(tmc2,mc1)
abs(eigen(mc2)$val)


dn1<-10
#n1<-1000+dn1
#nn1<-n1-kk
gm1<-3
sdm0<-Matrix(0,dm0,dm0,sparse=T)
for (i in 1:dm0){sdm0[i,i]=1}


cmm1<-function(x){
  cm1<-rowMeans(x)
  mm1<-matrix(rep(c(cm1),ncol(x)),nr=nrow(x),nc=ncol(x),byrow=F)
  return(x-mm1)
}

dg1<-function(n1,kk){
cov1<-matrix(0,dm0,dm0)
w1<-diag(sqrt(c(1.25,0.75)))%*%matrix(rnorm(dm0*n1),dm0)
y1<-w1
dy1<-w1
for (i in 5:n1){
    y1[,i]<-y1[,i]+(pi1+bb1)%*%y1[,i-1]-(bb1-bb2)%*%y1[,i-2]-(bb2-bb3)%*%y1[,i-3]-bb3%*%y1[,i-4]
    dy1[,i]<-y1[,i]-y1[,i-1]
}

y1<-y1[,-(1:dn1)]
dy1<-dy1[,-(1:dn1)]

ly1<-y1[,(kk):(ncol(y1)-1)]
dy0<-lg1(dy1,y1,kk)
ldy0<-dy0[-(1:dm0),]
dy0<-dy0[1:dm0,]

return(list(dy0,ly1,ldy0))
}


lag1<-function(gm1=3){
nn1<-ncol(dy0)
mm1<-diag(nn1)-t(ly1)%*%solve(ly1%*%t(ly1)+log(nn1))%*%ly1

tmdy0<-dy0%*%mm1
tmldy0<-ldy0%*%mm1

mdy0<-cmm1(tmdy0)
mldy0<-cmm1(tmldy0)

hb1<-(mdy0%*%t(mldy0)/nn1)%*%solve(mldy0%*%t(mldy0)/nn1+log(nn1)/nn1)
mhb1<-matrix(hb1,nr=dm0*dm0,byrow=F)
wb1<-(apply(abs(mhb1),2,max))^(-gm1)

cy0<-c(mdy0)
kx1<-t(mldy0)%x%sdm0
gid1<-rep(1:kk,each=dm0*dm0)
fit1<-lsgl::fit(kx1,cy0,intercept=F,grouping=gid1,lambda=0.01,alpha=0,groupWeights=wb1)

bic1<-0
for (i in 1:length(fit1$beta)){
    err1<-matrix(cy0-kx1%*%fit1$beta[[i]],nr=dm0)
    sig1<-err1%*%t(err1)/ncol(err1)
    bic1[i]<-log(det(sig1))+length(fit1$beta[[i]]@x)*log(nn1)/nn1
}

return(list(fit1,bic1))
}
############################################
#kk=3
lg1<-function(dy1,y1,kk=3){
if (kk>0){
dy0<-dy1[,-(1:kk)]
for (j in 1:kk){
    tdy1<-dy1[,(kk+1-j):(ncol(dy1)-j)]
    dy0<-rbind(dy0,tdy1)
    rm(tdy1)
}
}
return(dy0)
}
##############################################
rank1<-function(gm1=3){
nn1<-ncol(dy0)
mm1<-diag(nn1)-t(ldy0)%*%solve(ldy0%*%t(ldy0)+log(nn1))%*%ldy0

tmdy0<-dy0%*%mm1
tmly1<-ly1%*%mm1

mdy0<-cmm1(tmdy0)
mly1<-cmm1(tmly1)


hp1<-(mdy0%*%t(mly1)/nn1)%*%solve(mly1%*%t(mly1)/nn1)
dhp2<-qr1(t(hp1))

wr1<-sqrt(rowSums(dhp2$R^2))^(-gm1)
pmdy0<-t(dhp2$P)%*%mdy0
qmly1<-t(dhp2$Q)%*%mly1
ss1<-(dhp2$Q)

cy0<-c(pmdy0)
#kx1<-t(qmly1)%x%diag(dm0)
kx1<-t(qmly1)%x%sdm0
gid1<-rep(1:dm0,each=dm0)

fit1<-lsgl::fit(kx1,cy0,intercept=F,grouping=gid1,lambda=0.01,alpha=0,groupWeights=wr1)

bic1<-0
for (i in 1:length(fit1$beta)){
    err1<-matrix(cy0-kx1%*%fit1$beta[[i]],nr=dm0)
    sig1<-err1%*%t(err1)/ncol(err1)
    bic1[i]<-log(det(sig1))+length(fit1$beta[[i]]@x)*log(nn1)/nn1
}
return(list(fit1,bic1,ss1))
#return(fit1)
}

set.seed(2017)
sn1<-c(100,200,300,400)
np1<-5000
kk<-5
#df=10, rho=0.0,0.2,0.4,0.6
#df=15, rho=0.0,0.2,0.4,0.6

rk1<-matrix(0,nr=np1,nc=length(sn1))
pl1<-matrix(0,nr=np1,nc=length(sn1))


for (j in 1:length(sn1)){
for (i in 1:np1){#change here
print(c(i,j, rnorm(1)))
l1<-dg1(n1=sn1[j],kk);dy0<-l1[[1]];ly1<-l1[[2]];ldy0<-l1[[3]];nn1<-ncol(dy0)
rf1<-rank1(gm1=3)
d1<-30
l1<-which.min(rf1[[2]][-(1:d1)])+d1
rk1[i,j]<-length(rf1[[1]]$beta[l1][[1]]@x)/dm0

pf1<-lag1(gm1=4)
l1<-which.min(pf1[[2]][-(1:d1)])+d1
#pl1[i,j]<-length(pf1[[1]]$beta[l1][[1]]@x)/(dm0*dm0)
x1<-pf1[[1]]$beta[l1][[1]]
mx1<-matrix(x1,nr=dm0*dm0,byrow=F)
cmx1<-colSums(mx1!=0)
id1<-max(which(cmx1!=0))
pl1[i,j]<-id1

print(c(rk1[i,j],pl1[i,j]))

rm(l1,dy0,ly1,ldy0,rf1,pf1,x1,mx1,cmx1)
}
save(rk1,pl1,file="liao1g3g4.RData")
}
save(rk1,pl1,file="liao1g3g4.RData")



rk1<-matrix(0,nr=np1,nc=length(sn1))
pl1<-matrix(0,nr=np1,nc=length(sn1))


for (j in 1:length(sn1)){
for (i in 1:np1){#change here
print(c(i,j, rnorm(1)))
l1<-dg1(n1=sn1[j],kk);dy0<-l1[[1]];ly1<-l1[[2]];ldy0<-l1[[3]];nn1<-ncol(dy0)
rf1<-rank1(gm1=3)
d1<-30
l1<-which.min(rf1[[2]][-(1:d1)])+d1
rk1[i,j]<-length(rf1[[1]]$beta[l1][[1]]@x)/dm0

pf1<-lag1(gm1=5)
l1<-which.min(pf1[[2]][-(1:d1)])+d1
#pl1[i,j]<-length(pf1[[1]]$beta[l1][[1]]@x)/(dm0*dm0)
x1<-pf1[[1]]$beta[l1][[1]]
mx1<-matrix(x1,nr=dm0*dm0,byrow=F)
cmx1<-colSums(mx1!=0)
id1<-max(which(cmx1!=0))
pl1[i,j]<-id1

print(c(rk1[i,j],pl1[i,j]))

rm(l1,dy0,ly1,ldy0,rf1,pf1,x1,mx1,cmx1)
}
save(rk1,pl1,file="liao1g3g5.RData")
}
save(rk1,pl1,file="liao1g3g5.RData")

q()

