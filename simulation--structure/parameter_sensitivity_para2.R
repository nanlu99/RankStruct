source("dpFunc.R")
source("dpFuncContinue.R")

## check the sensitivity of the results with respect to the hyperparameters
## para-2
set.seed(123)
library(Rsolnp)
library(pracma)


##parameters##
h=0.02 #0.05
nitem=10
TT=100 #50 #50, 100,200,500,1000
m=TT*2 #number of pairwise comparison
mmcpt=9 #number of candidate change points
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
a0=c(rep(2.5/nitem,nitem*0.2),rep(0.625/nitem,nitem*0.8))
for(i in 1:(nitem*0.2)){
S[1:TT,i]=a0[i]+0.6/nitem*sin(7*pi*tt[1:TT]/ttlen)
}
for(i in (nitem*0.2+1):(nitem)){
S[1:TT,i]=a0[i]-0.15/nitem*sin(7*pi*tt[1:TT]/ttlen)
}
##generate pi^*(t) phaseII
a0=c(rep(1.5/nitem,nitem*0.2),rep(0.65/nitem,nitem*0.2))
for(i in 1:(nitem*0.2)){
S[(TT+1):(2*TT),i]=a0[i]-0.4/nitem*sin(15*pi*tt[(TT+1):(2*TT)]/ttlen)
}-for(i in (nitem*0.2+1):(nitem*0.4)){
S[(TT+1):(2*TT),i]=a0[i]+2.5/nitem*(tt[(TT+1):(2*TT)]/ttlen-1/2)^(1/10)
}
b0=(rep(1,TT)-rowSums(S[(TT+1):(2*TT),1:(nitem*0.4)]))/(nitem*0.6)
for(i in (nitem*0.4+1):nitem){
S[(TT+1):(2*TT),i]=b0
}
#S[,c(3,6)]<-S[,c(6,3)]
cpt=c(1/2)##true change point
plot((1:m)/m, S[,1],ylim=c(0,0.5))
lines((1:m)/m, S[,3])
lines((1:m)/m, S[,5])

##run##
B=50 #500 #500
RE=matrix(0,ncol=8,nrow=B)
gammals=c(0.002*(1:5),0.02*(1:5),0.2*(1:5)) #),0.1*(1:10) *(1:3),0.1*(1:3))
lambdals=c(0.02*(1:5),0.2*(1:5),2*(1:5)) #1:3 #sort(c(1e-5,5e-5*(1:2),5e-4*(1:2),5e-3*(1:2),0.05*(1:2),0.5*(1:2)),decreasing=TRUE) #

re_sensitivity = array(dim=c(B,length(gammals),length(lambdals)))
for(bct in 1:B){
print("No.")
print(bct)
RES=generate(S)


P=PFunc(RES,tt)
PLess=PFuncLess(RES,tt)
dpregamma=dpFunc(RES,tt,P,PLess,gammals,lambdals)

for(ii in 1:(length(gammals))){
for(jj in 1:(length(lambdals))){
re_sensitivity[bct, ii, jj] = length(dpregamma[[ii]][[jj]][[1]])
cpestmp=dpregamma[[ii]][[jj]][[1]]/(mmcpt+1)
print("change point!!!!!!!!!!!!!!!!")
print(cpestmp)
print(re_sensitivity[bct, ii, jj])
}}

print("#################################")
}

mean_re_sensitivity <- apply(re_sensitivity, c(2, 3), mean)

rownames(mean_re_sensitivity) <- gammals
colnames(mean_re_sensitivity) <- lambdals

# Print the resulting matrix with row and column names
print(mean_re_sensitivity)

file_name <- "res/mean_re_sensitivity_para2.csv"

# Save the matrix to a CSV file
write.csv(mean_re_sensitivity, file = file_name)

# Confirm that the file has been saved
absolute_path <- normalizePath(file_name)

# Print the absolute path
cat("The matrix has been saved to:", absolute_path)
