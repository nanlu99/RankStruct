source("dpFunc.R")
source("dpFuncContinue.R")
source("gpFunc.R")

set.seed(123)
library(Rsolnp)
library(pracma)

##parameters##
h=0.02
nitem=10
TT=100 #100, 200, 500, 1000
m=TT*3 #number of pairwise comparison
mmcpt=8 #number of candidate change points
mmept=5 #number of estimation points we use between two candidates #18,27,36
mmgpt=10 #number of estimation points used for group in a interval
mm0=(mmcpt+1)*mmept #all points used for estimation
ttlen=m*nitem
tt=m*nitem*(0:(m-1))/(m-1)
S=matrix(nrow=m,ncol=nitem) #score at the observed time
RES=array(dim=c(nitem,nitem,m)) #pairwise comparison result for (i,j,m)
#Porg=array(dim=c(nitem,nitem,mm0)) #original P & P0 
#P0org=array(dim=c(nitem,nitem,mm0))
T=diag(rep(1,nitem))
for(i in 1:(nitem-1)){
T[i,i+1]=-1}
T[nitem,]=rep(1,nitem)
Tinv=solve(T)


##generate pi^*(t) phaseI
a0=c(rep(2/nitem,nitem*0.3),rep(1/nitem,nitem*0.3))
for(i in 1:(nitem*0.3)){
S[1:TT,i]=a0[i]+0.3/nitem*sin(18*pi*tt[1:TT]/ttlen)
}
for(i in (nitem*0.3+1):(nitem*0.6)){
S[1:TT,i]=a0[i]-0.2/nitem*sin(18*pi*tt[1:TT]/ttlen)
}
b0=(rep(1,TT)-rowSums(S[1:TT,1:(nitem*0.6)]))/(nitem*0.4)
for(i in (nitem*0.6+1):nitem){
S[1:TT,i]=b0
}
##generate pi^*(t) phaseII
a0=c(rep(1.5/nitem,nitem*0.5))
#c(rep(1.2/nitem,nitem*0.5))
for(i in 1:(nitem*0.5)){
S[(TT+1):(2*TT),i]=a0[i]+0.2/nitem*sin(18*pi*tt[(TT+1):(2*TT)]/ttlen)
#a0[i]+0.7/nitem*sin(18*pi*tt[(TT+1):(2*TT)]/ttlen)
}
b0=(rep(1,TT)-rowSums(S[(TT+1):(2*TT),1:(nitem*0.5)]))/(nitem*0.5)
for(i in (nitem*0.5+1):nitem){
S[(TT+1):(2*TT),i]=b0
}
##generate pi^*(t) phaseIII
a0=c(rep(2/nitem,nitem*0.3),rep(1/nitem,nitem*0.3))
for(i in 1:(nitem*0.3)){
S[(2*TT+1):(3*TT),i]=a0[i]+0.3/nitem*sin(18*pi*tt[(2*TT+1):(3*TT)]/ttlen)
}
for(i in (nitem*0.3+1):(nitem*0.6)){
S[(2*TT+1):(3*TT),i]=a0[i]-0.2/nitem*sin(18*pi*tt[(2*TT+1):(3*TT)]/ttlen)
}
b0=(rep(1,TT)-rowSums(S[(2*TT+1):(3*TT),1:(nitem*0.6)]))/(nitem*0.4)
for(i in (nitem*0.6+1):nitem){
S[(2*TT+1):(3*TT),i]=b0
}
S[,c(3,6)]<-S[,c(6,3)]
cpt=c(1/3,2/3)##true change point
plot((1:m)/m, S[,1],ylim=c(0,0.5))
lines((1:m)/m, S[,4])
lines((1:m)/m, S[,7])


##run##
B=500 #500
RE=matrix(0,ncol=8,nrow=B)
for(bct in 1:B){
print("No.")
print(bct)
RES=generate(S)
gammals=c(0.002*(1:5),0.02*(1:5),0.2*(1:5)) #),0.1*(1:10) *(1:3),0.1*(1:3))
lambdals=c(0.02*(1:5),0.2*(1:5),2*(1:5)) #1:3 #sort(c(1e-5,5e-5*(1:2),5e-4*(1:2),5e-3*(1:2),0.05*(1:2),0.5*(1:2)),decreasing=TRUE) #

##dp and 10-fold CV
CVre=array(dim=c(10,length(gammals),length(lambdals)))
for(CVct in 1:10){
CVrmls=CVct+(0:(m/10-1))*10 #the removed list of CV

RESCVtrain=RES[,,-CVrmls]
ttCVtrain=tt[-CVrmls]
PCVtrain=PFunc(RESCVtrain,ttCVtrain) #calculate observed P
PCVtrainLess=PFuncLess(RESCVtrain,ttCVtrain)

RESCVtest=RES[,,CVrmls]
ttCVtest=tt[CVrmls]
PCVtest=PFunc(RESCVtest,ttCVtest)

dpre=dpFunc(RESCVtrain,ttCVtrain,PCVtrain,PCVtrainLess,gammals,lambdals) #change point result of dp
dpre
rfre=nlk4CV(PCVtest,dpre)
CVre[CVct,,]=rfre  #[[1]]
}

meanCVre=apply(CVre,c(2,3),mean)
pos=which(meanCVre==min(meanCVre,na.rm=T),arr.ind=TRUE)
gammaCV=gammals[pos[ceiling(nrow(pos)/2),1]]  #get the CV gamma
lambdaCV=lambdals[pos[ceiling(nrow(pos)/2),2]]
print("##########CVpara:")
print(gammaCV)
print(lambdaCV)
P=PFunc(RES,tt)
PLess=PFuncLess(RES,tt)
dpregamma=dpFunc(RES,tt,P,PLess,gammaCV,lambdaCV)
cpesgamma=dpregamma[[1]][[1]][[1]]/(mmcpt+1)
print("change point!!!!!!!!!!!!!!!!")
print(cpesgamma)


##gp
gpes=GPFunc(RES)/(mmcpt+1)
print(gpes)


##judge
if(length(cpesgamma)>length(cpt)){
RE[bct,4]=1
}else if(length(cpesgamma)==length(cpt)){
RE[bct,3]=1
}else 
{RE[bct,2]=1}
RE[bct,1]=hausdorff_dist(cpt, cpesgamma)

if(length(gpes)>length(cpt)){
RE[bct,8]=1
}else if(length(gpes)==length(cpt)){
RE[bct,7]=1
}else 
{RE[bct,6]=1}
RE[bct,5]=hausdorff_dist(cpt, gpes)

print("#################################")
print(RE[bct,])
}

RE
colMeans(RE)

save(RE, file = paste0("res/Para1.RData"))


res=matrix(nrow=0,ncol=8)
res<-rbind(res,colMeans(RE))
colnames(res)<-rep(c("Hausdorff","hat_p<p","hat_p=p","hat_p>p"),2)
rownames(res)<-TT*h
res
write.csv(res,"res/Para1.csv")