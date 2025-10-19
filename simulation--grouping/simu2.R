source("simuFunc.R")

##simu result--para@2
set.seed(123)
library(pracma)
##parameters
h=0.05
B=500 #10
eps0=0.001 #points for od estimation
eps=0.01 #pertubation of score
Lnitem=c(10,20,50)
Lmm=30 #10*(1:4) #number of estimation ponints we use
Lm=c(10,20,30,50,100,200) #50+50*(1:3) #number of pairwise comparison
mjdg=array(-100,dim=c(length(Lnitem),length(Lm),20))


for(cnitem in 1:length(Lnitem)){
nitem=Lnitem[cnitem]
T=diag(rep(1,nitem))
for(i in 1:(nitem-1)){
T[i,i+1]=-1}
T[nitem,]=rep(1,nitem)
Tinv=solve(T)

#for(cmm in 1:length(Lmm)){
mm=Lmm #[cmm]
for(cm in 1:length(Lm)){
m=Lm[cm]
print("nitem")
print(nitem)
print("mm")
print(mm)
print("m")
print(m)
ttlen=m*nitem
tt=m*nitem*seq(0,1,length.out=m)
S=matrix(nrow=m,ncol=nitem)
RES=array(dim=c(nitem,nitem,m))


##generate pi^*(t)
a0=c(rep(1.9/nitem,nitem*0.3),rep(0.1/nitem,nitem*0.3))
for(i in 1:(nitem*0.1)){
S[,i]=a0[i]+0.5/nitem*sin(3*pi*tt/ttlen)+eps/nitem
}
for(i in (nitem*0.1+1):(nitem*0.2)){
S[,i]=a0[i]+0.5/nitem*sin(3*pi*tt/ttlen)
}
for(i in (nitem*0.2+1):(nitem*0.3)){
S[,i]=a0[i]+0.5/nitem*sin(3*pi*tt/ttlen)-eps/nitem
}

for(i in (nitem*0.3+1):(nitem*0.4)){
S[,i]=a0[i]+0.6/nitem*atan(pi*tt/ttlen)+eps/nitem
}
for(i in (nitem*0.4+1):(nitem*0.5)){
S[,i]=a0[i]+0.6/nitem*atan(pi*tt/ttlen)
}
for(i in (nitem*0.5+1):(nitem*0.6)){
S[,i]=a0[i]+0.6/nitem*atan(pi*tt/ttlen)-eps/nitem
}

b0=(rep(1,m)-rowSums(S[,1:(nitem*0.6)]))/(nitem*0.4)
for(i in (nitem*0.6+1):(nitem*0.8)){
S[,i]=b0+eps/nitem
}
for(i in (nitem*0.8+1):nitem){
S[,i]=b0-eps/nitem
}

#true score
tttr=mm*nitem*(1:(mm))/(mm+1)
ttlentr=mm*nitem
Str=matrix(nrow=mm,ncol=nitem) 

for(i in 1:(nitem*0.1)){
Str[,i]=a0[i]+0.5/nitem*sin(3*pi*tttr/ttlentr)+eps/nitem
}
for(i in (nitem*0.1+1):(nitem*0.2)){
Str[,i]=a0[i]+0.5/nitem*sin(3*pi*tttr/ttlentr)
}
for(i in (nitem*0.2+1):(nitem*0.3)){
Str[,i]=a0[i]+0.5/nitem*sin(3*pi*tttr/ttlentr)-eps/nitem
}

for(i in (nitem*0.3+1):(nitem*0.4)){
Str[,i]=a0[i]+0.6/nitem*atan(pi*tttr/ttlentr)+eps/nitem
}
for(i in (nitem*0.4+1):(nitem*0.5)){
Str[,i]=a0[i]+0.6/nitem*atan(pi*tttr/ttlentr)
}
for(i in (nitem*0.5+1):(nitem*0.6)){
Str[,i]=a0[i]+0.6/nitem*atan(pi*tttr/ttlentr)-eps/nitem
}

b0=(rep(1,mm)-rowSums(Str[,1:(nitem*0.6)]))/(nitem*0.4)
for(i in (nitem*0.6+1):(nitem*0.8)){
Str[,i]=b0+eps/nitem
}
for(i in (nitem*0.8+1):nitem){
Str[,i]=b0-eps/nitem
}


#true score grouped
Strgp=matrix(nrow=mm,ncol=nitem) 

for(i in 1:(nitem*0.1)){
Strgp[,i]=a0[i]+0.5/nitem*sin(3*pi*tttr/ttlentr)
}
for(i in (nitem*0.1+1):(nitem*0.2)){
Strgp[,i]=a0[i]+0.5/nitem*sin(3*pi*tttr/ttlentr)
}
for(i in (nitem*0.2+1):(nitem*0.3)){
Strgp[,i]=a0[i]+0.5/nitem*sin(3*pi*tttr/ttlentr)
}

for(i in (nitem*0.3+1):(nitem*0.4)){
Strgp[,i]=a0[i]+0.6/nitem*atan(pi*tttr/ttlentr)
}
for(i in (nitem*0.4+1):(nitem*0.5)){
Strgp[,i]=a0[i]+0.6/nitem*atan(pi*tttr/ttlentr)
}
for(i in (nitem*0.5+1):(nitem*0.6)){
Strgp[,i]=a0[i]+0.6/nitem*atan(pi*tttr/ttlentr)
}

b0=(rep(1,mm)-rowSums(Str[,1:(nitem*0.6)]))/(nitem*0.4)
for(i in (nitem*0.6+1):(nitem*0.8)){
Strgp[,i]=b0
}
for(i in (nitem*0.8+1):nitem){
Strgp[,i]=b0
}



##
jdg=matrix(-1,nrow=B,ncol=20)
for(b in 1:B){
print(b)
RES=generate(S)
hatpi=ADGL(RES)
hatpirf=KRCref(RES,hatpi)
hatpist=SG(RES)
hatpistrf=SGref(RES,hatpist)


#Str-hatpi
#Str-hatpist
#judgeNO(hatpi)
#judgeNO(hatpist)

KRCpi=KRC(RES)

jdg[b,]=c(judgeNO(hatpi),judgeNO(hatpirf),judgeNO(hatpist),
judgeNO(hatpistrf),judgeNO(KRCpi))
#print(b)
print(jdg[b,])
}
#jdg
mjdg[cnitem,cm,]=colMeans(jdg,na.rm=T)
print(mjdg[cnitem,cm,])
}}
#}


#############
#mjdgre<-matrix(nrow=0,ncol=20)

#for(cnitem in 1:length(Lnitem)){
mjdgtmp=mjdg #[cnitem,,,]
mjdg2=array(-100,dim=c(length(Lm),20,length(Lnitem)))
for(kk in 1:length(Lnitem)){
mjdg2[,,kk]=mjdgtmp[kk,,]}
mjdg3=apply(mjdg2,2,c)
mjdg3
#mjdgre<-rbind(mjdgre,mjdg3)
#}
# write.csv(mjdg3,file = "res/para2-mjdg3.csv")



## process results
res <- mjdg3[,c(1+4*(0:4),2+4*(0:4),3,11,4,12)]
res
cnames <- c("GKRC","GKRC(rf)","GRC","GRC(rf)","KRC")

colnames(res) <- c(rep(cnames,2),rep(c("GKRC","GRC"),2))
rownames(res) <- rep(h*Lm,3)

res <- res[c(-1,-2),c(-3,-8)] #res[c(9,10,14,15),c(-3,-8)] 
colnames(res)[c(3,7)]<-rep("GRC",2)
res

save(mjdg3, res, file = paste0("res/para2.RData"))
write.csv(res,"res/para2.csv")