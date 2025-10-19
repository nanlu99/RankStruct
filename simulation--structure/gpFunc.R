GPFunc=function(RES){
gpeach=matrix(nrow=(mmcpt+1),ncol=nitem)

for(i0 in 0:mmcpt){
intgp=(i0+(1:mmgpt)/(mmgpt+1))/(mmcpt+1)
lsgp=which( (tt/(m*nitem))>=(i0/(mmcpt+1)) & 
(tt/(m*nitem))<=((i0+1)/(mmcpt+1)) )

ttgp=tt[lsgp]
RESgp=RES[,,lsgp]
hatpi=ADGL(RESgp,ttgp,intgp) 

gp=rep(1,nitem)
if(norm(hatpi[,2]-hatpi[,1],'2')/(mmept)>=1e-10)
{gp[2]=2}
for(i in 3:(nitem)){ 
temp=apply(hatpi[,1:(i-1)], 2, 
function(col) norm(hatpi[,i]-col,'2')/(mmept))
temp0=min(temp)
if(temp0<1e-10) #1e-10
{gp[i]=max(gp[which(temp==temp0)])}
else{gp[i]=max(gp)+1}
}
print(gp)
gpeach[i0+1,]=gp
}

gpcpes=NULL
for(i0 in 1:(mmcpt)){
if(!identical(gpeach[i0,],gpeach[i0+1,])){
gpcpes=c(gpcpes,i0)}}
return(gpcpes)
}



##adaptive group lasso##
ADGL=function(RESx,ttgp,intgp){
RES=RESx

##estimate order
ix=intgp[1]
phi=dnorm(abs(ttgp-ix)/ttlen,sd = h)
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
P=array(dim=c(nitem,nitem,mmgpt))
Pod=P
t0=m*nitem*intgp
for(k in 1:mmgpt){
ix=t0[k]
phi=dnorm(abs(ttgp-ix)/ttlen,sd = h)
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

for(k in 1:mmgpt){
Pod[,,k]=P[oddy,oddy,k]
}

##calculate XX and YY
XX=matrix(0,nrow=nitem*mmgpt,ncol=(nitem-1)*mmgpt)
Theta=rep(-1,nitem*mmgpt)
YY=rep(-1,nitem*mmgpt)
P1=array(dim=c(nitem,nitem,mmgpt))
II=diag(rep(1,nitem))
for(k in 1:mmgpt){
P1[,,k]=(t(Pod[,,k])-II)%*%Tinv}
for(i in 1:mmgpt){
YY[((i-1)*nitem+1):(i*nitem)]=(II-t(Pod[,,i]))%*%rep(1/nitem,nitem)}
for(i in 1:mmgpt){
#XX[((i-1)*nitem+1):(i*nitem),((i-1)*(nitem-1)+1):(i*(nitem-1))]=P1[,1:(nitem-1),i]
XX[((i-1)*nitem+1):(i*nitem),seq(i,mmgpt*(nitem-1),by=mmgpt)]=P1[,1:(nitem-1),i]}

#get the pre-estimator without penalty
coefs <- glm.fit(XX, YY, family = gaussian())$coefficients
#adaptive group lasso
#index=c(rep(1:(nitem-1),mm))
index=rep(1:(nitem-1),each=mmgpt)
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
        cr2[iii]=CR5gp(YY,XX,thetatemp,coefs) 
        iii=iii+1
      }
      lam=lams2[which.min(cr2)]
      hatthetab=c(t(matrix(fit2$beta[,which.min(cr2)],mmgpt)))
      #which.min(cr2)
      #matrix(hatthetab,5)


##give pi estimation using opt results
hatpib=matrix(ncol=nitem,nrow=mmgpt)
for(k in 1:mmgpt){
hatpib[k,]=Tinv%*%c(hatthetab[(1+(k-1)*(nitem-1)):(k*(nitem-1))],0)+rep(1/nitem,nitem)
}
hatpi=signif(hatpib,digits=5) #hatpi=trunc(hatpib*1e10)/1e10
pihato=matrix(nrow=mmgpt,ncol=nitem)
for(i in 1:mmgpt){
pihato[i,oddy]=hatpi[i,]}
return(pihato)}


CR5gp=function(YY,XX,theta,coefs){
  thetam=matrix(theta,ncol=nitem-1,nrow=mmgpt)
  #th=c(t(matrix(theta,mm)))
  norlogc=log(mean((YY-XX%*%theta)^2)+0.1*var(YY)) #0.1
  thetanorm=apply(thetam,2,function(x) norm(x,type='2'))
  dfa=sum(thetanorm>1e-10)
  if(dfa==0)
    return(Inf)
coefsm=matrix(coefs,ncol=nitem-1,nrow=mmgpt,byrow=F) #should be byrow=F
coefsnorm=apply(coefsm,2,function(x) norm(x,type='2'))
dfb=sum(thetanorm/coefsnorm*(mmgpt-1))
dff=dfa+dfb
  zz=(nitem*mmgpt)*norlogc+(log(nitem*mmgpt))*dff #+2*0*log(choose((nitem*mm),(ceiling(dff)))) #2
  return(zz)
}


