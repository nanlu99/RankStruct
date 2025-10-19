source("Func.R")
source("FuncContinue.R")
source("PFunc.R")


library(tidyverse)
library(Rsolnp)
library(pracma)


AD<-read.csv("data/0919nba.csv")
AD=AD[3:13291,]
AD=AD%>%select(A,B,AisWin,time) #sample
seasonSep=1/11*(0:11)
tradeSep=c(0.053411355,0.144048805,0.236429283,0.332295817,
0.417704183,0.508341633,0.598979084,0.689616534,0.781997012,
0.869148406,0.959785857)
seasonPick=(2015:2019)-2009+1


tt=sort(unique(AD$time))
h=1/11
nitem=30
ttsep=sort(c(seasonSep[seasonPick[1]:length(seasonSep)],
tradeSep[seasonPick]))
AD<-AD[(AD$time>min(ttsep))&(AD$time<max(ttsep)),]
mmcpt=length(ttsep)-2 #number of candidate change points
mmept=5 #number of estimation points we use between two candidates #18,27,36
mmgpt=10 #number of estimation points used for group in a interval
mm0=(mmcpt+1)*mmept #all points used for estimation
#ttlen=m*nitem
#tt=m*nitem*(0:(m-1))/(m-1)
t0save=NULL
for(i in 1:(length(ttsep)-1)){
t0tmp=ttsep[i]+((1:mmept)/(mmept+1))*(ttsep[i+1]-ttsep[i])
t0save=c(t0save,t0tmp)
}

#PP=array(dim=c(nitem,nitem,mm)) #P in KRC
#Mp=matrix(ncol=nitem,nrow=length(lix)) #no panel pi estimation 
T=diag(rep(1,nitem))
for(i in 1:(nitem-1)){
T[i,i+1]=-1}
T[nitem,]=rep(1,nitem)
Tinv=solve(T)



gammals=c(0.002*(1:5),0.02*(1:5),0.2*(1:5)) #1e-4,1e-2,0.02,0.05,
lambdals=c(0.02*(1:5),0.2*(1:5),2*(1:5)) #c(0.02*(1:5),0.2*(1:5),2*(1:5)) #1:3 
#gamma=1e-4

##dp and 10-fold CV
CVre=array(dim=c(10,length(gammals),length(lambdals)))
for(CVct in 1:10){
CVrmtmp0=CVct+(0:(length(tt)/10))*10 #the removed list of CV
CVrmtmp=CVrmtmp0[which(CVrmtmp0<=length(tt))]
CVrmls=which(AD$time %in% tt[CVrmtmp])

ADtrain=AD[-CVrmls,]
#ttCVtrain=tt[-CVrmls]
PCVtrain=PFunc(ADtrain) #calculate observed P
PCVtrainLess=PFuncLess(ADtrain)

ADtest=AD[CVrmls,]
#ttCVtest=tt[CVrmls]
PCVtest=PFunc(ADtest)

dpre=dpFunc(ADtrain,PCVtrain,PCVtrainLess,gammals,lambdals) #change point result of dp
dpre
rfre=nlk4CV(PCVtest,dpre)
CVre[CVct,,]=rfre  #[[1]]
print("number of CV:")
print(CVct)
}

meanCVre=apply(CVre,c(2,3),mean)
pos=which(meanCVre==min(meanCVre,na.rm=T),arr.ind=TRUE)
gammaCV=gammals[pos[ceiling(nrow(pos)/2),1]]  #get the CV gamma
lambdaCV=lambdals[pos[ceiling(nrow(pos)/2),2]]
print("##########CVpara:")
print(gammaCV)
print(lambdaCV)
P=PFunc(AD)
PLess=PFuncLess(AD)
dpregamma=dpFunc(AD,P,PLess,gammaCV,lambdaCV)
cpesgamma=ttsep[dpregamma[[1]][[1]][[1]]+1] #/(mmcpt+1)
print("change point!!!!!!!!!!!!!!!!")
print(cpesgamma)


####################
####################

#xx=rep(((0:(mmcpt))/(mmcpt+1)),each=mmept)+((1:mmept)/(mmept+1))/(mmcpt+1)
xx=t0save
yy<-matrix(nrow=0,ncol=nitem)
for(i in 1:(length(cpesgamma)+1)){
yytmp=dpregamma[[1]][[1]][[3]][[i]] #[3:10,]
yy<-rbind(yy,yytmp)
}

plot(1,1,xlim=c(min(xx),max(xx)),ylim=c(0,max(yy)),type='n',
las=1,xlab='t',ylab=expression(pi),cex.axis=1.5,cex.lab=2.5)
for(i in 1:ncol(yy)){
lines(xx,yy[,i],lty=5,col = i+1,lwd=4)
} 
linexxls=cpesgamma
for(linexx in linexxls){
abline(v=linexx,lty=2)
#abline(v=(linexx-(1/(mm+1))),lty=2)
#rect(xleft=(linexx-(1/(mm+1))),xright=linexx,ybottom=0,ytop=max(yy),
#density = 10,border=NA)
}

save.image(file = "output/nba1519.RData")
