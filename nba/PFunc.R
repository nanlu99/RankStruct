PFunc=function(ADob){
Porg=array(dim=c(nitem,nitem,mm0)) #original P
t0=t0save
ttob=unique(ADob$time)
for(k in 1:mm0){
ix=t0[k]
  phi=dnorm(abs(ttob-ix),sd = h)
  P=matrix(0,nrow = nitem,ncol=nitem)
  N=P
  
  for(i in 1:nrow(ADob)){
    P[ADob$A[i],ADob$B[i]]=P[ADob$A[i],ADob$B[i]]+as.numeric(1-ADob$AisWin[i])*phi[which(ttob==ADob$time[i])]
    N[ADob$A[i],ADob$B[i]]=N[ADob$A[i],ADob$B[i]]+phi[which(ttob==ADob$time[i])]
    P[ADob$B[i],ADob$A[i]]=P[ADob$B[i],ADob$A[i]]+(as.numeric(ADob$AisWin[i]))*phi[which(ttob==ADob$time[i])]
    N[ADob$B[i],ADob$A[i]]=N[ADob$B[i],ADob$A[i]]+phi[which(ttob==ADob$time[i])]
  }
  
  P=P/N
  diag(P)=0
nitemtmp=length(unique(c(ADob$A,ADob$B)))
if(nitemtmp!=nitem){
print("wrong!!!")}
  P=P/nitem
  P=P+diag(1-rowSums(P))
P[is.na(P)]=0
#  a=eigen(t(P)) #no panel P
#  p_hat=Re(a$vectors[,1])/sum(Re(a$vectors[,1]))
#  Mp[gg,]=p_hat
  Porg[,,k]=P
print(ix)
}
return(Porg)
}


########################
PFuncLess = function(ADob1){
Porg=list()
for(PL in 0:mmcpt){
PorgtmpL=list()
for(PR in (PL+1):(mmcpt+1)){
Porgtmp=array(dim=c(nitem,nitem,mmgpt)) #original P
lstmp=which( (ADob1$time)>=ttsep[PL+1] & 
(ADob1$time)<=ttsep[PR+1] )
#ttLR=tttmp[lstmp]
ADLR=ADob1[lstmp,]
htmp=h #(PR-PL)*

#t0=m*nitem*( sort(rep(((1:mmept)/(mmept+1))/(mmcpt+1),each=mmcpt+1)+((0:(mmcpt))/(mmcpt+1))) )
#t0=m*nitem*( rep(((PL:(PR-1))/(mmcpt+1)),each=mmept)+((1:mmept)/(mmept+1))/(mmcpt+1) )
t0=ttsep[PL+1]+(1:mmgpt)/(1+mmgpt)*(ttsep[PR+1]-ttsep[PL+1])
ttLR=unique(ADLR$time)

for(k in 1:mmgpt){
#print(k)
ix=t0[k]
phi=dnorm(abs(ttLR-ix),sd = htmp)
  P=matrix(0,nrow = nitem,ncol=nitem)
  N=P
  
  for(i in 1:nrow(ADLR)){
    P[ADLR$A[i],ADLR$B[i]]=P[ADLR$A[i],ADLR$B[i]]+as.numeric(1-ADLR$AisWin[i])*phi[which(ttLR==ADLR$time[i])]
    N[ADLR$A[i],ADLR$B[i]]=N[ADLR$A[i],ADLR$B[i]]+phi[which(ttLR==ADLR$time[i])]
    P[ADLR$B[i],ADLR$A[i]]=P[ADLR$B[i],ADLR$A[i]]+(as.numeric(ADLR$AisWin[i]))*phi[which(ttLR==ADLR$time[i])]
    N[ADLR$B[i],ADLR$A[i]]=N[ADLR$B[i],ADLR$A[i]]+phi[which(ttLR==ADLR$time[i])]
  }
  
  P=P/N
  diag(P)=0
nitemtmp=length(unique(c(ADLR$A,ADLR$B)))
if(nitemtmp!=nitem){
print("wrong!!!")}
  P=P/nitem
  P=P+diag(1-rowSums(P))
P[is.na(P)]=0
Porgtmp[,,k]=P
}
PorgtmpL[[PR]]=Porgtmp
print(PR)}
print("left")
print(PL)
Porg[[PL+1]]=PorgtmpL}
return(Porg)
}
