##simu function using sparsegl
rm(list=ls())
library(glmnet)
library(sparsegl)


#########################
##functions
##generation function##
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

##BIC for dynamic group
##adaptive group lasso##
ADGL=function(RESx){
RES=RESx

##estimate order
ix=eps0
phi=dnorm(abs(tt-ix)/ttlen,sd = h)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(nrow = nitem,ncol=nitem)
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((RES[i,j,])*phi)
A[i,j]=sRESphi/sphin
A[j,i]=(sphi-sRESphi)/sphin
}}
A[is.na(A)]=0
A=A+diag(1-rowSums(A))
p<-eigen(t(A))
x1=Re(p$vectors[,1])/sum(Re(p$vectors[,1])) 
pidy=x1
oddy=sort(x1,index.return = TRUE)$ix


##estimate P
P=array(dim=c(nitem,nitem,mm))
Pod=P
t0=tttr/ttlentr*ttlen #m*nitem*(0:(mm-1))/mm
for(k in 1:mm){
ix=t0[k]
phi=dnorm(abs(tt-ix)/ttlen,sd = h)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(nrow = nitem,ncol=nitem)
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((RES[i,j,])*phi)
A[i,j]=sRESphi/sphin
A[j,i]=(sphi-sRESphi)/sphin
}}
A[is.na(A)]=0
A=A+diag(1-rowSums(A))
P[,,k]=A}

for(k in 1:mm){
Pod[,,k]=P[oddy,oddy,k]
}

##calculate XX and YY
XX=matrix(0,nrow=nitem*mm,ncol=(nitem-1)*mm)
Theta=rep(-1,nitem*mm)
YY=rep(-1,nitem*mm)
P1=array(dim=c(nitem,nitem,mm))
II=diag(rep(1,nitem))
for(k in 1:mm){
P1[,,k]=(t(Pod[,,k])-II)%*%Tinv}
for(i in 1:mm){
YY[((i-1)*nitem+1):(i*nitem)]=(II-t(Pod[,,i]))%*%rep(1/nitem,nitem)}
for(i in 1:mm){
#XX[((i-1)*nitem+1):(i*nitem),((i-1)*(nitem-1)+1):(i*(nitem-1))]=P1[,1:(nitem-1),i]
XX[((i-1)*nitem+1):(i*nitem),seq(i,mm*(nitem-1),by=mm)]=P1[,1:(nitem-1),i]}

#get the pre-estimator without penalty
coefs <- glm.fit(XX, YY, family = gaussian())$coefficients
#adaptive group lasso
#index=c(rep(1:(nitem-1),mm))
index=rep(1:(nitem-1),each=mm)
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
        cr2[iii]=CR5(YY,XX,thetatemp,coefs) 
        iii=iii+1
      }
      lam=lams2[which.min(cr2)]
      hatthetab=c(t(matrix(fit2$beta[,which.min(cr2)],mm)))
      #which.min(cr2)
      #matrix(hatthetab,5)


##give pi estimation using opt results
hatpib=matrix(ncol=nitem,nrow=mm)
for(k in 1:mm){ #RIGHT!!!already process above
hatpib[k,]=Tinv%*%c(hatthetab[(1+(k-1)*(nitem-1)):(k*(nitem-1))],0)+rep(1/nitem,nitem)
}
hatpi=signif(hatpib,digits=5)
pihato=matrix(nrow=mm,ncol=nitem)
for(i in 1:mm){
pihato[i,oddy]=hatpi[i,]}
return(pihato)}


CR5=function(YY,XX,theta,coefs){
  thetam=matrix(theta,ncol=nitem-1,nrow=mm)
  #th=c(t(matrix(theta,mm)))
  norlogc=log(mean((YY-XX%*%theta)^2)+0.1*var(YY)) #0.1
  thetanorm=apply(thetam,2,function(x) norm(x,type='2'))
  dfa=sum(thetanorm>1e-10)
  if(dfa==0)
    return(Inf)
coefsm=matrix(coefs,ncol=nitem-1,nrow=mm,byrow=F) #should be byrow=F
coefsnorm=apply(coefsm,2,function(x) norm(x,type='2'))
dfb=sum(thetanorm/coefsnorm*(mm-1))
dff=dfa+dfb
  zz=(nitem*mm)*norlogc+(log(nitem*mm))*dff #+2*0*log(choose((nitem*mm),(ceiling(dff)))) #2
  return(zz)
}



##static group##
SG=function(RESx){
RES=RESx
##estimate P
Pst=matrix(nrow = nitem,ncol=nitem)
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
Pst[i,j]=mean(RES[i,j,])/nitem
Pst[j,i]=1/nitem-Pst[i,j]
}}
Pst[is.na(Pst)]=0
Pst=Pst+diag(1-rowSums(Pst))

##give an order
p<-eigen(t(Pst))
x1=Re(p$vectors[,1])/sum(Re(p$vectors[,1])) 
pist=x1
odst=sort(x1,index.return = TRUE)$ix
Pstod=Pst[odst,odst]

##calculate XX and YY
II=diag(rep(1,nitem))
XXst=matrix(0,nrow=nitem,ncol=(nitem-1))
YYst=rep(-1,nitem)
YYst=(II-t(Pstod))%*%rep(1/nitem,nitem)
XXst=(t(Pstod)-II)%*%Tinv[,1:(nitem-1)]

##optimization
#get the pre-estimator without penalty
library(pracma)
coefst=lm(YYst~XXst+0)$coefficients

##using SRC
cv.fit<-glmnet(x=XXst,y=YYst,intercept =FALSE,penalty.factor = 1/abs(coefst))
cv.fit$beta
#plot(cv.fit, label = TRUE)
lams=cv.fit$lambda
hb=rep(-1,length(cv.fit$lambda))
iii=1      
for(lam in lams){
thetatemp=coef(cv.fit,s=lam)[-1]
hb[iii]=Hbic(YYst,XXst,thetatemp) 
iii=iii+1
}
lam=lams[which.min(hb)]
hatthetast=coef(cv.fit,s=lam)[-1]

##give pi estimation using opt results
hatpist1=Tinv%*%c(hatthetast,0)+rep(1/nitem,nitem)
hatpist=matrix(rep(t(hatpist1),mm),nrow=mm,ncol=nitem,byrow = T)
hatpi=signif(hatpist,digits=5)

pihato=matrix(nrow=mm,ncol=nitem)
for(i in 1:mm){
pihato[i,odst]=hatpi[i,]}
return(pihato)}


Hbic=function(YYst,XXst,hatthetast){
  ms=sum(abs(hatthetast)>1e-5&abs(hatthetast)!=Inf)
  if(ms==0)
    return(Inf)
  return(nitem*log(mean((XXst%*%hatthetast-YYst)^2)+0.2*var(YYst))+log(nitem)*ms+2*0.5*log(choose(nitem,ms)))
}

##KRC##
KRC=function(RESx){
RES=RESx
KRCpi=matrix(ncol=nitem,nrow=mm)
##estimate KRCpi
t0=tttr/ttlentr*ttlen #m*nitem*(0:(mm-1))/mm
for(k in 1:mm){
ix=t0[k]
phi=dnorm(abs(tt-ix)/ttlen,sd = h)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(nrow = nitem,ncol=nitem)
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((RES[i,j,])*phi)
A[i,j]=sRESphi/sphin
A[j,i]=(sphi-sRESphi)/sphin
}}
A[is.na(A)]=0
A=A+diag(1-rowSums(A))
p<-eigen(t(A))
x1=Re(p$vectors[,1])/sum(Re(p$vectors[,1])) 
KRCpi[k,]=x1}
hatpi=signif(KRCpi,digits=5)
return(hatpi)}



##judge no order##
judgeNO=function(hatpix){
#estimation of original pi
pihato=matrix(nrow=mm,ncol=nitem)
for(i in 1:mm){
pihato[i,]=hatpix[i,]}

#Kendall tau cor
Kd=rep(-1,mm)
for(i in 1:mm){
Kd[i]=cor(pihato[i,],Strgp[i,],method='kendall')
}
d1=mean(Kd) 

#MSE
mse=rep(-1,mm)
for(i in 1:mm){
mse[i]=norm(pihato[i,]-Str[i,],'2')/norm(Str[i,],'2')
}
d2=mean(mse) 

#sensitivity and specificity
sen=0
nsen=0
spe=0
nspe=0
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){          
if( ((i %in% (1:(nitem*0.3)))&(j %in% (1:(nitem*0.3)))) ||
 ((i %in% ((nitem*0.3+1):(nitem*0.6)))&(j %in% ((nitem*0.3+1):(nitem*0.6)))) ||
 ((i %in% ((nitem*0.6+1):nitem))&(j %in% ((nitem*0.6+1):nitem)))
){
#if(norm(pihato[,i]-pihato[,j],'2')/mm>1e-10){
#print(i)
#print(j)}
sen=sen+(norm(pihato[,i]-pihato[,j],'2')/mm<1e-10) #1e-10
nsen=nsen+1}
else{
spe=spe+(norm(pihato[,i]-pihato[,j],'2')/mm>1e-10)
nspe=nspe+1}
}}
sen/nsen
spe/nspe
return(c(d1,d2,sen/nsen,spe/nspe))}


##KRC for refit##
KRCref=function(RESx, hatpi){
##give the group estimation
gp=rep(1,nitem)
if(norm(hatpi[,2]-hatpi[,1],'2')/mm>=1e-10)
{gp[2]=2}
for(i in 3:(nitem)){ 
temp=apply(hatpi[,1:(i-1)], 2, 
function(col) norm(hatpi[,i]-col,'2')/mm)
temp0=min(temp)
if(temp0<1e-10) #1e-10
{gp[i]=max(gp[which(temp==temp0)])}
else{gp[i]=max(gp)+1}
}
##estimate P for the grouped items
RES=RESx
KRCpigp=matrix(ncol=nitem,nrow=mm)
t0=tttr/ttlentr*ttlen #m*nitem*(0:(mm-1))/mm
for(k in 1:mm){
ix=t0[k]
phi=dnorm(abs(tt-ix)/ttlen,sd = h)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(0,nrow = max(gp),ncol=max(gp))
A1=matrix(0,nrow = max(gp),ncol=max(gp))
A2=matrix(0,nrow = max(gp),ncol=max(gp))
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((RES[i,j,])*phi)
A1[gp[i],gp[j]]=A1[gp[i],gp[j]]+sRESphi
A2[gp[i],gp[j]]=A2[gp[i],gp[j]]+sphin
A1[gp[j],gp[i]]=A1[gp[j],gp[i]]+(sphi-sRESphi)
A2[gp[j],gp[i]]=A2[gp[j],gp[i]]+sphin
}}
A=A1/A2
A[is.na(A)]=0
A=A+diag(1-rowSums(A))
p<-eigen(t(A))
x1=Re(p$vectors[,1])/sum(Re(p$vectors[,1])) 
x1
for(i in 1:nitem){
KRCpigp[k,i]=x1[gp[i]]}
KRCpigp=KRCpigp/rowSums(KRCpigp)
KRCpigp
}
hatpi=signif(KRCpigp,digits=5)
return(hatpi)}


##static group for refit##
SGref=function(RESx, hatpist){
##give the group estimation
gp=rep(1,nitem)
if(norm(hatpist[,2]-hatpist[,1],'2')/mm>=1e-10)
{gp[2]=2}
for(i in 3:(nitem)){ 
temp=apply(hatpist[,1:(i-1)], 2, 
function(col) norm(hatpist[,i]-col,'2')/mm)
temp0=min(temp)
if(temp0<1e-10) #1e-10
{gp[i]=max(gp[which(temp==temp0)])}
else{gp[i]=max(gp)+1}
}
##estimate P for the grouped items
RES=RESx
SGpigp=rep(0,nitem)
SGhatpigp=matrix(ncol=nitem,nrow=mm)
Pstgp=matrix(0,nrow = max(gp),ncol=max(gp))
phi=rep(1,m)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(0,nrow = max(gp),ncol=max(gp))
A1=matrix(0,nrow = max(gp),ncol=max(gp))
A2=matrix(0,nrow = max(gp),ncol=max(gp))
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((RES[i,j,])*phi)
A1[gp[i],gp[j]]=A1[gp[i],gp[j]]+sRESphi
A2[gp[i],gp[j]]=A2[gp[i],gp[j]]+sphin
A1[gp[j],gp[i]]=A1[gp[j],gp[i]]+(sphi-sRESphi)
A2[gp[j],gp[i]]=A2[gp[j],gp[i]]+sphin
}}
A=A1/A2
A[is.na(A)]=0
A=A+diag(1-rowSums(A))
p<-eigen(t(A))
x1=Re(p$vectors[,1])/sum(Re(p$vectors[,1])) 
for(i in 1:nitem){
SGpigp[i]=x1[gp[i]]}
SGpigp=SGpigp/sum(SGpigp)
for(i in 1:mm){
SGhatpigp[i,]=SGpigp}
SGhatpigp=signif(SGhatpigp,digits=5)
return(SGhatpigp)
}




