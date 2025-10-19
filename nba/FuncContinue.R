##KRC for refit##
#change time points
KRCref=function(ADtmp,intervalRef,l,r,gp){
##estimate P for the grouped items
lstmp=which( (ADtmp$time)>=ttsep[l+1] & 
(ADtmp$time)<=ttsep[r+1] )
ADlr=ADtmp[lstmp,]

t0=intervalRef
ttlr=unique(ADlr$time)
KRCpigp=matrix(ncol=nitem,nrow=length(t0))
for(k in 1:length(t0)){
ix=t0[k]
phi=dnorm(abs(ttlr-ix),sd = h)
  P=matrix(0,nrow = max(gp),ncol=max(gp))
  N=P
  
  for(i in 1:nrow(ADlr)){
    P[gp[ADlr$A[i]],gp[ADlr$B[i]]]=P[gp[ADlr$A[i]],gp[ADlr$B[i]]]+as.numeric(1-ADlr$AisWin[i])*phi[which(ttlr==ADlr$time[i])]
    N[gp[ADlr$A[i]],gp[ADlr$B[i]]]=N[gp[ADlr$A[i]],gp[ADlr$B[i]]]+phi[which(ttlr==ADlr$time[i])]
    P[gp[ADlr$B[i]],gp[ADlr$A[i]]]=P[gp[ADlr$B[i]],gp[ADlr$A[i]]]+(as.numeric(ADlr$AisWin[i]))*phi[which(ttlr==ADlr$time[i])]
    N[gp[ADlr$B[i]],gp[ADlr$A[i]]]=N[gp[ADlr$B[i]],gp[ADlr$A[i]]]+phi[which(ttlr==ADlr$time[i])]
  }
  
  P=P/N
  diag(P)=0
nitemtmp=length(unique(c(ADlr$A,ADlr$B)))
if(nitemtmp!=nitem){
print("wrong!!!")}
  P=P/nitem
  P=P+diag(1-rowSums(P))
P[is.na(P)]=0
p<-eigen(t(P))
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
ttptmm=floor(c(0,gmre[[1]],mmcpt+1)) #time partition using mmcpt
rsdO=0
for(j in 1:(length(ttsep[gmre[[1]]])+1)){
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


