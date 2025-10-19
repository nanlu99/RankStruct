rm(list=ls())
library(glmnet)
library(sparsegl)


#######################
##generation function##
#######################
generate=function(Sx){
S=Sx
##generate y
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
ds=S[,j]/(S[,j]+S[,i])
RES[i,j,]=as.numeric(runif(m)<ds)
RES[j,i,]=1-RES[i,j,] 
}}
return(RES)
}



########################
##calculate P function##for likelihood calculation
########################
PFunc = function(REStmp,tttmp){
Porg=array(dim=c(nitem,nitem,mm0)) #original P

#estimate P
#t0=m*nitem*( sort(rep(((1:mmept)/(mmept+1))/(mmcpt+1),each=mmcpt+1)+((0:(mmcpt))/(mmcpt+1))) )
t0=m*nitem*( rep(((0:(mmcpt))/(mmcpt+1)),each=mmept)+((1:mmept)/(mmept+1))/(mmcpt+1) )

for(k in 1:mm0){
ix=t0[k]
phi=dnorm(abs(tttmp-ix)/ttlen,sd = h)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(nrow = nitem,ncol=nitem)
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((REStmp[i,j,])*phi)
A[i,j]=sRESphi/sphin
A[j,i]=(sphi-sRESphi)/sphin
#P0org[i,j,k]=sRESphi #for another likelihood
#P0org[j,i,k]=sphi-sRESphi
}}
A[is.na(A)]=0
A=A+diag(1-rowSums(A))
Porg[,,k]=A
#P0org[,,k]=P0org[,,k]/sum(P0org[,,k],na.rm=TRUE)
}

return(Porg)
}



########################
##calculate P function##for grouping from L to R
########################
PFuncLess = function(REStmp,tttmp){ 
Porg=list()
for(PL in 0:mmcpt){
PorgtmpL=list()
for(PR in (PL+1):(mmcpt+1)){
Porgtmp=array(dim=c(nitem,nitem,mmgpt)) #original P
lstmp=which( (tttmp/(m*nitem))>=(PL/(mmcpt+1)) & 
(tttmp/(m*nitem))<=(PR/(mmcpt+1)) )
ttLR=tttmp[lstmp]
RESLR=REStmp[,,lstmp]
htmp=h #(PR-PL)*h

#t0=m*nitem*( sort(rep(((1:mmept)/(mmept+1))/(mmcpt+1),each=mmcpt+1)+((0:(mmcpt))/(mmcpt+1))) )
#t0=m*nitem*( rep(((PL:(PR-1))/(mmcpt+1)),each=mmept)+((1:mmept)/(mmept+1))/(mmcpt+1) )
t0=m*nitem*( PL+(1:mmgpt)/(1+mmgpt)*(PR-PL) )/(mmcpt+1)

for(k in 1:mmgpt){
#print(k)
ix=t0[k]
phi=dnorm(abs(ttLR-ix)/ttlen,sd = htmp)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(nrow = nitem,ncol=nitem)
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((RESLR[i,j,])*phi)
A[i,j]=sRESphi/sphin
A[j,i]=(sphi-sRESphi)/sphin
#P0org[i,j,k]=sRESphi #for another likelihood
#P0org[j,i,k]=sphi-sRESphi
}}
A[is.na(A)]=0
A=A+diag(1-rowSums(A))
Porgtmp[,,k]=A
#P0org[,,k]=P0org[,,k]/sum(P0org[,,k],na.rm=TRUE)
}
PorgtmpL[[PR]]=Porgtmp}
Porg[[PL+1]]=PorgtmpL}
return(Porg)
}





#######################
##dynamic programming##
#######################
dpFunc=function(RESpt,ttpt,Pob,PobLess,gammals,lambdals){
#estimation result for all gamma -- 
#1.change point integer; 234.nog,pi in each interval
temp<-matrix(vector("list"),nrow=mmcpt+1,ncol=mmcpt+1)
temp1<-matrix(vector("list"),nrow=mmcpt+1,ncol=mmcpt+1)
#temp2<-matrix(vector("list"),nrow=mmcpt+1,ncol=mmcpt+1)
##core circulate: caculate dp
r_cand=1:(mmcpt+1)
l_cand=0:mmcpt
dp=matrix(Inf,nrow=mmcpt+1,ncol=mmcpt+1)
dp1=matrix(Inf,nrow=mmcpt+1,ncol=mmcpt+1)
dp2=matrix(Inf,nrow=mmcpt+1,ncol=mmcpt+1)
for(r in r_cand){
for(l in l_cand){
if(l>=r){break}
inter=m*nitem*( rep(((l:(r-1))/(mmcpt+1)),each=mmept)+((1:mmept)/(mmept+1))/(mmcpt+1) )
#(l*mmept+1):(r*mmept)  #l:r
PobLR=Pob[,,(l*mmept+1):(r*mmept)] #
#interLess=m*nitem*( l+(1:mmgpt)/(1+mmgpt)*(r-l) )/(mmcpt+1)
PobLRLess=PobLess[[l+1]][[r]]

estimate=opt_group(RESpt,ttpt,PobLRLess,PobLR,inter,l,r) #group for Pob on inter/mm0
#print(estimate)
dp1[(l+1),r]=estimate[[1]] #neg loglike
dp2[(l+1),r]=estimate[[2]]*estimate[[3]] #total group number
temp[[(l+1),r]]=estimate[[5]] #nog
temp1[[(l+1),r]]=estimate[[4]] #estimated pi
#temp2[[(l+1),r]]=estimate[[7]] #ordered es coff
}}


cp_est_all<-list()
gammact=1
for(gamma in gammals){
print("gamma")
print(gamma)

cp_est_lambda=list()
#lambdact=1
for(lambda1ct in 1:length(lambdals)){ #,0.03,0.05,0.1,0.02
lambda1=lambdals[lambda1ct]
if(lambda1<=gamma){
next}
print(lambda1)
dp=dp1+lambda1*dp2 #0.05+lambda1*dp2

p=rep(-1,mmcpt+1) #c(-1,rep(1,mm-1)) #record how to get the minimum value for a certain point
b=dp[1,1:(mmcpt+1)] #c(0,rep(Inf,mmcpt+1)) #minimum dp-value from left end to a certain point

for(r in r_cand[-1]){
for(l in l_cand[-1]){
if(l>=r){break}
bb=b[l]+gamma+dp[l+1,r]
if(bb<b[r]){
b[r]=bb
p[r]=l}}}
#p

k=mmcpt+1
cp_est<-vector()
dtmp=c(mmcpt+1)
while(k>0){
d=p[k]
if(d>0){
cp_est<-c(cp_est,d) }
k=d
#print(d)
dtmp=c(dtmp,d)
}
dtmp=sort(dtmp)
#dtmp
cp_est=sort(cp_est) #change point estimation -- integer
cp_Est<-cp_est/(mmcpt+1) #change point estimation -- [0,1]
print(cp_Est)
dtmp=c(0,cp_est,mmcpt+1)

gatmp=list()
gatmp[[1]]=cp_est
gatmpA=list()
gatmpB=list()
gatmpC=list()
for(i in 1:(length(dtmp)-1)){
gatmpA[[i]]=temp[[(dtmp[i]+1),dtmp[i+1]]] #nog,pi in each interval
gatmpB[[i]]=temp1[[(dtmp[i]+1),dtmp[i+1]]]
#gatmpC[[i]]=temp2[[(dtmp[i]+1),dtmp[i+1]]]
}
gatmp[[2]]=gatmpA
gatmp[[3]]=gatmpB
#gatmp[[4]]=gatmpC
cp_est_lambda[[lambda1ct]]<-gatmp

#lambdact=lambdact+1
}
cp_est_all[[gammact]]<-cp_est_lambda
gammact=gammact+1
}#}
return(cp_est_all) 
}




########################
#Functions for DP method
########################
##group optimazation for a given interval## 
#1.neg loglike; 2.nog; 3.inter_len; 4.hatpi; 5.gp res;
#6. estimated order; 7. estimated ordered coff
#group for Pob on inter/mm0
opt_group=function(RESpt,ttpt,PobLRLess,PobLR,interval,l,r){
#print(interval)
inter_len=length(interval)
adglre=ADGLCV(RESpt,ttpt,PobLRLess,interval,l,r)
hatpi=adglre[[1]]
lkvalue=lk(PobLR,hatpi,inter_len) #negative loglike in interval for hatpi
unique_rows <- unique(signif(t(hatpi), digits = 3))
nog <- nrow(unique(unique_rows)) #number of group
REtemp=list()
REtemp[[1]]=lkvalue
REtemp[[2]]=nog
REtemp[[3]]=inter_len
REtemp[[4]]=hatpi
REtemp[[5]]=adglre[[2]] #group result
#REtemp[[6]]=adglre[[3]] #estimated order
#REtemp[[7]]=adglre[[4]] #estimated coff
return(REtemp)
}


##likelihood function
lk=function(PobLR,x,inter_len){
x=t(x)
z=0
for(k in 1:inter_len){
for(i in 1:nitem){
for(j in 1:nitem){
if(i!=j){
z=z+PobLR[i,j,k]*log(x[j,k]/(x[i,k]+x[j,k]))
#z=z+P0[i,j,k]*log(x[j,k]/(x[i,k]+x[j,k]))
}}}
}
return(-z) ##return negative log likelihood
}




########################
##adaptive group lasso##
########################
ADGLCV=function(RESpt,ttpt,Pobtrain,interval,l,r){ 
#1.hatpi; 2.hatgp
P=Pobtrain #[,,interval]

##estimate order
A=P[,,1] #order at the first point of the interval
p<-eigen(t(A))
x1=Re(p$vectors[,1])/sum(Re(p$vectors[,1])) 
oddy=sort(x1,index.return = TRUE)$ix

Pod=P
inter_len=mmgpt
for(k in 1:inter_len){
Pod[,,k]=P[oddy,oddy,k]
}

##calculate XX and YY
XX=matrix(0,nrow=nitem*inter_len,ncol=(nitem-1)*inter_len)
Theta=rep(-1,nitem*inter_len)
YY=rep(-1,nitem*inter_len)
P1=array(dim=c(nitem,nitem,inter_len))
II=diag(rep(1,nitem))
for(k in 1:inter_len){
P1[,,k]=(t(Pod[,,k])-II)%*%Tinv}
for(i in 1:inter_len){
YY[((i-1)*nitem+1):(i*nitem)]=(II-t(Pod[,,i]))%*%rep(1/nitem,nitem)}
for(i in 1:inter_len){
#XX[((i-1)*nitem+1):(i*nitem),((i-1)*(nitem-1)+1):(i*(nitem-1))]=P1[,1:(nitem-1),i]
XX[((i-1)*nitem+1):(i*nitem),seq(i,inter_len*(nitem-1),by=inter_len)]=P1[,1:(nitem-1),i]}

#get the pre-estimator without penalty
coefs <- glm.fit(XX, YY, family = gaussian())$coefficients
#adaptive group lasso
#index=c(rep(1:(nitem-1),mm))
index=rep(1:(nitem-1),each=inter_len)
ind_weights <- 1 / abs(coefs)
grp_weights <- numeric(length(unique(index)))
for (i in unique(index)) {
  grp_weights[i] <-  1 / sqrt(sum(coefs[which(index == i)]^2))
}

fit2<-sparsegl(XX,YY,
     group = index, family = "gaussian", 
     pf_group = grp_weights,  intercept = FALSE,asparse = 0,standardize=FALSE)

      lams2=fit2$lambda
      cr2=rep(-1,length(fit2$lambda))
      iii=1      
      for(lam2 in lams2){
        thetatemp=fit2$beta[,iii]
        cr2[iii]=CR5(YY,XX,thetatemp,coefs,inter_len) 
        iii=iii+1
      }
      lam=lams2[which.min(cr2)]


      hatthetab=c(t(matrix(fit2$beta[,which.min(cr2)],inter_len)))
      #which.min(cr2)
      #matrix(hatthetab,5)

##give pi estimation using opt results
hatpib=matrix(ncol=nitem,nrow=inter_len)
for(k in 1:inter_len){ #RIGHT!!!already process above
hatpib[k,]=Tinv%*%c(hatthetab[(1+(k-1)*(nitem-1)):(k*(nitem-1))],0)+rep(1/nitem,nitem)
}
hatpi=signif(hatpib,digits=5) #trunc(hatpib*1e10)/1e10
pihato=matrix(nrow=inter_len,ncol=nitem)
for(i in 1:inter_len){
pihato[i,oddy]=hatpi[i,]}

hatpi=signif(hatpi,digits=3)
cluster=rep(1,nitem)
for(i in 2:(nitem)){
if(identical(hatpi[,i-1],hatpi[,i])){
cluster[i]=cluster[i-1]
}else{cluster[i]=cluster[i-1]+1}
}
gre=rep(-1,nitem)
gre[oddy]=cluster
gre=max(gre)-gre+1

piref=KRCref(RESpt,ttpt,interval,l,r,gre)
REtemp=list()
REtemp[[1]]=piref
REtemp[[2]]=gre
#REtemp[[3]]=oddy
#REtemp[[4]]=fit2$beta[,which.min(cr2)]
return(REtemp)} 


CR5=function(YY,XX,theta,coefs,inter_len){
  thetam=matrix(theta,ncol=nitem-1,nrow=inter_len)
  #th=c(t(matrix(theta,inter_len)))
  norlogc=log(mean((YY-XX%*%theta)^2)+0.1*var(YY)) #0.1
  thetanorm=apply(thetam,2,function(x) norm(x,type='2'))
  dfa=sum(thetanorm>1e-10)
  if(dfa==0)
    return(Inf)
coefsm=matrix(coefs,ncol=nitem-1,nrow=inter_len,byrow=F) #should be byrow=F
coefsnorm=apply(coefsm,2,function(x) norm(x,type='2'))
dfb=sum(thetanorm/coefsnorm*(inter_len-1))
dff=dfa+dfb
  zz=(nitem*inter_len)*norlogc+(log(nitem*inter_len))*dff #+2*0*log(choose((nitem*mm),(ceiling(dff)))) #2
  return(zz)
}


