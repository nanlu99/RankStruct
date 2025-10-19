library(ggplot2)

# --- Plot Text Size Settings ---
legend_text_size <- 25
legend_title_size <- 25
axis_title_size <- 35
axis_text_size <- 30
plot_title_size <- 25
# -----------------------------

load("output/nba1519.RData")

Mp=matrix(ncol=nitem,nrow=length(t0save)) #no panel pi estimation 

ADob=AD
Porg=array(dim=c(nitem,nitem,mm0)) #original P
t0=t0save
ttob=unique(ADob$time)
gg=1
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
  a=eigen(t(P)) #no panel P
  p_hat=Re(a$vectors[,1])/sum(Re(a$vectors[,1]))
  Mp[gg,]=p_hat
gg=gg+1
  Porg[,,k]=P
print(ix)
}

dfplot<-data.frame(plott0=rep(t0,each=nitem),y=as.vector(t(Mp)),group=team)
ggplot(dfplot,aes(x=plott0,y=y,group=group,color=factor(group)))+
geom_line(linewidth= 2,alpha=0.5)+
scale_x_continuous("Season",breaks=ttsep[c(TRUE,FALSE)],labels=partitionls[c(TRUE,FALSE)])+
  labs(y = expression(hat(pi)),color="Team") +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
#        panel.background = element_rect(fill = "transparent",colour = NA),
 #       plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
#        legend.position = c(.05,.935),
#        legend.box.background = element_rect(color="black"),
        legend.text = element_text(size = legend_text_size),
        legend.title=element_text(size = legend_title_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        axis.text.x = element_text(size = axis_text_size),
        axis.text.y = element_text(size = axis_text_size),
        plot.title = element_text(size = plot_title_size, face = "bold")) 

ggsave("plt/noPanelPlot.pdf",units='cm',height=20,width=30)
