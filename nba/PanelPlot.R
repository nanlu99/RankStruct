library(ggplot2)

# --- Plot Text Size Settings ---
legend_text_size <- 25
legend_title_size <- 25
axis_title_size <- 35
axis_text_size <- 30
plot_title_size <- 25
# -----------------------------

load("output/nba1519.RData")

partitionls=c()
yearPick=c(min(seasonPick-1),seasonPick)+2008
for(i in yearPick){
partitionls=c(partitionls,as.character(i),paste(i,"TradeDDL",sep=""))
}
partitionls=partitionls[-length(partitionls)]



dfplot<-data.frame(plott0=rep(t0save,each=nitem),y=as.vector(t(yy)),group=team)
ggplot(dfplot,aes(x=plott0,y=y,group=group,color=factor(group)))+
geom_line(linewidth= 2,alpha=0.5)+
scale_x_continuous("Season",breaks=ttsep[c(TRUE,FALSE)],labels=partitionls[c(TRUE,FALSE)])+
  labs(y = expression(hat(pi)),color="Team") +
geom_vline(xintercept=cpesgamma,linetype=3,linewidth=1.5)+  
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

ggsave("plt/PanelPlot.pdf",units='cm',height=20,width=30)
