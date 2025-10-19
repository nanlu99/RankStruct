##KRC for refit##
#change time points
KRCref=function(RESx,ttpt,intervalRef,l,r,gp){
##estimate P for the grouped items
lstmp=which( (ttpt/(m*nitem))>=(l/(mmcpt+1)) & 
(ttpt/(m*nitem))<=(r/(mmcpt+1)) )
ttlr=ttpt[lstmp]
RESlr=RESx[,,lstmp]

t0=intervalRef
KRCpigp=matrix(ncol=nitem,nrow=length(t0))
for(k in 1:length(t0)){
ix=t0[k]
phi=dnorm(abs(ttlr-ix)/ttlen,sd = h)
sphi=sum(phi) 
sphin=sphi*nitem
A=matrix(0,nrow = max(gp),ncol=max(gp))
A1=matrix(0,nrow = max(gp),ncol=max(gp))
A2=matrix(0,nrow = max(gp),ncol=max(gp))
for(i in 1:(nitem-1)){
for(j in (i+1):nitem){
sRESphi=sum((RESlr[i,j,])*phi)
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




#################################################
##given hat pi, CV
##calculate the negative likelihood for CV
#################################################
nlk4CV=function(PobtestOut,ptre){ #parameters: RES & partition result
#return: 1.CV error for each gamma

resgamma=matrix(nrow=length(gammals),ncol=length(lambdals)) 
#neg likelihood result for each para
#refgamma=list() #refit result for each gamma
for(jj in 1:length(gammals)){
for(jj2 in 1:length(lambdals)){
if(lambdals[jj2]<=gammals[jj]){
next}
#hatpirf=matrix(ncol=nitem,nrow=mm) #refit hat pi 
gmre=ptre[[jj]][[jj2]] #result for each gamma
ttpt=c(0,gmre[[1]]/(mmcpt+1),1) #time partition
ttptmm=floor(c(0,gmre[[1]],mmcpt+1)) #time partition using mmcpt
rsdO=0
for(j in 1:(length(ttpt)-1)){
#RESint=ceiling(dim(REStest)[3]*ttpt[j]+1e-3):floor(dim(REStest)[3]*ttpt[j+1]) #result interval
#gp=gmre[[j+1]]
intRefit=(ttptmm[j]*mmept+1):(ttptmm[j+1]*mmept) 
#interval

PobLRCV=PobtestOut[,,intRefit]
inter_lenCV=length(intRefit)
lkvalue=lk(PobLRCV,gmre[[3]][[j]],inter_lenCV) #

rsdO=rsdO+lkvalue
}
resgamma[jj,jj2]=rsdO
}}

return(resgamma) #,ptre[[which.min(resgamma)]][1],ptre[[which.min(resgamma)]][2]))
}


