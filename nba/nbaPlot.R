##from 0# --- Plot Text Size Settings ---
legend_text_size <- 30
axis_text_size <- 20
axis_title_size <- 35
title_size <- 50
legend_key_height_cm <- 4
# -----------------------------aPlot
library(dplyr)
#install.packages("gridExtra")
library(janitor)
library(ggplot2)
library(gridExtra)
library(grid)

# --- Plot Text Size Settings ---
legend_text_size <- 40
axis_text_size <- 25
axis_title_size <- 42
# -----------------------------

load("output/nba1519.RData")

team=c("ATL","BOS","CLE","NOP","CHI","DAL","DEN","GSW","HOU",
"LAC","LAL","MIA","MIL","MIN","BKN","NYK","ORL","IND","PHI",
"PHX","POR","SAC","SAS","OKC","TOR","UTA","MEM","WAS","DET","CHA")

partitionls=c()
yearPick=c(min(seasonPick-1),seasonPick)+2008
for(i in yearPick){
partitionls=c(partitionls,as.character(i),paste(i,"TradeDDL",sep=""))
}
partitionls=partitionls[-length(partitionls)]

titlels=c()
posi=match(cpesgamma,ttsep)
posi=c(1,posi,length(ttsep))
for(posict in 1:(length(posi)-1)){
titletmp=paste(partitionls[posi[posict]],"-",
partitionls[posi[posict+1]],sep="")
titlels=c(titlels,titletmp)
}


for(i in 1:(length(cpesgamma)+1)){
print(cpesgamma[i])
print((which(ttsep==cpesgamma[i])-1)/2+min(yearPick))
titletmp=(which(ttsep==cpesgamma[i])-1)/2+min(yearPick)

cpestmp=c(ttsep[1],cpesgamma,ttsep[length(ttsep)])

tttmp0=cpestmp[i] #ttsep[which(ttsep==cpesgamma[i])]
tttmp1=cpestmp[i+1] #ttsep[which(ttsep==cpesgamma[i+1])]
ADtmp<-AD[which(AD$time>tttmp0 & AD$time<tttmp1),c(1,2)]

pairres=as.data.frame.matrix(table(ADtmp,useNA="ifany"))
pairall=pairres+t(pairres)
mat=pairres/pairall
mat=data.matrix(mat)
mat_df<-as.data.frame(as.table(mat))
mat_df$odVar1=team[mat_df$Var1]
mat_df$odVar2=team[mat_df$Var2]
#min(pairall==t(pairall))

plot1<-ggplot(mat_df,aes(x=odVar1,y=odVar2,fill=Freq))+
scale_x_discrete("Team A",expand=c(0,0))+#,labels=team)+
scale_y_discrete("Team B",expand=c(0,0))+#,labels=team)+
geom_tile()+
scale_fill_gradient(low="white",high="blue")+
theme(legend.title=element_blank(),
	legend.text=element_text(size=legend_text_size),
	axis.text.y=element_text(size=axis_text_size),
    axis.text.x=element_text(size=axis_text_size, angle = 45, hjust = 1),
	axis.title=element_text(size=axis_title_size),
    legend.key.height = unit(legend_key_height_cm, "cm"))



gptmp=dpregamma[[1]][[1]][[2]][[i]]
gpnumber=cumsum(as.data.frame(table(gptmp))$Freq)[-length(table(gptmp))]+0.5
numbertmp=as.data.frame(sort(gptmp,index.return = TRUE)$ix)
teamsort=team[sort(gptmp,index.return = TRUE)$ix]
numbertmp$nb=1:30

for(jj in 1:nrow(mat_df)){
mat_df$Var1new[jj]=which(numbertmp[,1]==mat_df$Var1[jj])
mat_df$Var2new[jj]=which(numbertmp[,1]==mat_df$Var2[jj])
}


plot2<-ggplot(mat_df,aes(x=Var1new,y=Var2new,fill=Freq))+
scale_x_continuous("Team A",expand=c(0,0),breaks=1:30,labels=teamsort)+
scale_y_continuous("Team B",expand=c(0,0),breaks=1:30,labels=teamsort)+
geom_tile()+
scale_fill_gradient(low="white",high="blue")+
geom_hline(yintercept=gpnumber,linewidth=1)+
geom_vline(xintercept=gpnumber,linewidth=1)+
theme(legend.title=element_blank(),
	legend.text=element_text(size=legend_text_size),
	axis.text.y=element_text(size=axis_text_size),
    axis.text.x=element_text(size=axis_text_size, angle = 45, hjust = 1),
	axis.title=element_text(size=axis_title_size),
    legend.key.height = unit(legend_key_height_cm, "cm"))

g=grid.arrange(plot1,plot2,ncol=2,top=textGrob(titlels[i],gp=gpar(fontsize=title_size)))
ggsave(paste("plt/nbaPlot",i,".pdf",sep=""), plot=g, units='cm', height=35, width=75)
}