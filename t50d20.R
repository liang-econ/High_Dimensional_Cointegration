rm(list=ls(all=T))
library(glmnet)
library(Matrix)
library(lsgl)
set.seed(2017)
source("qr.R")
####################################################
kk=3
rank0<-10
dm0<-50
m1<-matrix(rnorm(rank0*dm0),dm0)
m2<-matrix(rnorm(rank0*dm0),dm0)

q1<-qr.Q(qr(m1))
q2<-qr.Q(qr(m2))
#print(t(q2)%*%q1)
#ld0<-diag(c(2,2,2,1,1,rep(0,15)))
ld0<-diag(c(2,2,2,1,1,rep(1,5),rep(0,40)))
#ld0<-diag(c(2,1,1,1,1,rep(1,20),rep(0,5)))

a1<-ld0%*%q1-q2
b1<-q2
#print(eigen(t(a1)%*%a1)$val)
pi0<-a1%*%t(b1)
pi1<-pi0+diag(dm0)
#print(abs(eigen(pi1)$val))

q1<-a1
q2<-b1

#bb1<-diag(runif(dm0,-0.5,0.5))
bb1<-diag(rep(0,dm0))

tm1<-matrix(0,dm0*2,dm0*2)
tm1[(dm0+1):(dm0*2),1:dm0]=-bb1
tm1[1:dm0,(dm0+1):(dm0*2)]=diag(dm0)
tm1[(dm0+1):(dm0*2),(dm0+1):(dm0*2)]=pi1+bb1

va1<-abs(eigen(tm1)$val)
sum(round(va1,2)==1.00)
sum(round(va1,2)<1.00)

dn1<-10
n1<-1000+dn1
nn1<-n1-kk
gm1<-3
sdm0<-Matrix(0,dm0,dm0,sparse=T)
for (i in 1:dm0){sdm0[i,i]=1}


cmm1<-function(x){
  cm1<-rowMeans(x)
  mm1<-matrix(rep(c(cm1),ncol(x)),nr=nrow(x),nc=ncol(x),byrow=F)
  return(x-mm1)
}

dg1<-function(n1,kk){
w1<-matrix(rt(dm0*n1,df=20),dm0,n1)
y1<-w1
dy1<-w1
for (i in 3:n1){
    y1[,i]<-y1[,i]+(pi1+bb1)%*%y1[,i-1]-bb1%*%y1[,i-2]
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


lag1<-function(){
nn1<-ncol(dy0)
mm1<-diag(nn1)-t(ly1)%*%solve(ly1%*%t(ly1))%*%ly1

tmdy0<-dy0%*%mm1
tmldy0<-ldy0%*%mm1

mdy0<-cmm1(tmdy0)
mldy0<-cmm1(tmldy0)

hb1<-(mdy0%*%t(mldy0)/nn1)%*%solve(mldy0%*%t(mldy0)/nn1)
mhb1<-matrix(hb1,nr=dm0*dm0,byrow=F)
wb1<-(apply(abs(mhb1),2,max))^{-gm1}

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
kk=3
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
rank1<-function(){
nn1<-ncol(dy0)

mdy0<-cmm1(dy0)
mly1<-cmm1(ly1)

hp1<-(mdy0%*%t(mly1)/nn1)%*%solve(mly1%*%t(mly1)/nn1+1/nn1)
#dhp1<-qr(t(hp1))
#qr1<-qr.R(dhp1)
#qq1<-qr.Q(dhp1)
dhp2<-qr1(t(hp1))

wr1<-abs(diag(dhp2$R))^(-gm1)
pmdy0<-t(dhp2$P)%*%mdy0
qmly1<-t(dhp2$Q)%*%mly1

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
return(list(fit1,bic1))
#return(fit1)
}


#sn1<-c(300,500,800,1000)
#sn1<-c(1000,2000,3000,4000,5000)
sn1<-seq(1000,3000,by=500)

np1<-100
rk1<-matrix(0,nr=np1,nc=length(sn1))
pl1<-matrix(0,nr=np1,nc=length(sn1))

for (j in 1:length(sn1)){
for (i in 1:np1){#change np1 here
print(c(i,j,rnorm(1)))
l1<-dg1(n1=sn1[j],3);dy0<-l1[[1]];ly1<-l1[[2]];ldy0<-l1[[3]];nn1<-ncol(dy0)
rf1<-rank1()
d1<-30
l1<-which.min(rf1[[2]][-(1:d1)])+d1
rk1[i,j]<-length(rf1[[1]]$beta[l1][[1]]@x)/dm0
rm(l1,dy0,ly1,ldy0,rf1)
}
save(rk1,file="t50d20_rank.RData")

}

save(rk1,file="t50d20_rank.RData")

q()


