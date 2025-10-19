library(ggplot2)

xxa <- seq(0, 1/3-0.0001, 0.0001)
xxb <- seq(1/3+0.0001, 2/3-0.0001, 0.0001)
xxc <- seq(2/3+0.0001, 1, 0.0001)
xx <- c(xxa, xxb, xxc)

yy1a <- 0.2 + 0.03 * sin(18 * pi * xxa)
yy2a <- 0.1 - 0.02 * sin(18 * pi * xxa)
yy3a <- 0.025 - 0.0075 * sin(18 * pi * xxa)

yy1b <- 0.15 + 0.02 * sin(18 * pi * xxb)
yy2b <- 0.05 - 0.02 * sin(18 * pi * xxb)

yy1c <- 0.2 + 0.03 * sin(18 * pi * xxc)
yy2c <- 0.1 - 0.02 * sin(18 * pi * xxc)
yy3c <- 0.025 - 0.0075 * sin(18 * pi * xxc)

#3*yy1a+3*yy2a+4*yy3a
#5*yy1b+5*yy2b

yy1_3 <- c(yy1a, yy1b, yy1c)
yy4_5 <- c(yy2a, yy1b, yy2c)
yy6 <- c(yy2a, yy2b, yy2c)
yy7_10 <- c(yy3a, yy2b, yy3c)


#####################

data <- data.frame(xx, yy1_3, yy4_5, yy6, yy7_10)

ggplot(data, aes(x = xx)) +
  geom_line(aes(y = yy1_3, linetype = "Items1", color = "Items1", alpha="Items1"), size = 2) +
  geom_line(aes(y = yy4_5, linetype = "Items4", color = "Items4"), size = 2, alpha=0.9) +
  geom_line(aes(y = yy6, linetype = "Item6", color = "Item6", alpha="Item6"), size = 2, alpha=0.75) +
  geom_line(aes(y = yy7_10, linetype = "Items7", color = "Items7", alpha="Items7"), size =2, alpha=1) +
  labs(x = "t", y = expression(paste(pi, "*"))) +
  #theme_minimal() +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
#        legend.position = c(.05,.935),
#        legend.box.background = element_rect(color="black"),
legend.text = element_text(family = 'serif',size = 23),
legend.title=element_blank(),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
	  legend.key.size=unit(5,"lines")) +
  scale_linetype_manual(name = "", values =c(Items1 = "solid",Items4 = "dotdash",Item6 = "dashed",Items7 = "dotted"),
                        labels = c("Items 1-3", "Items 4-5", "Item 6", "Items 7-10")) +
  scale_color_manual(name = "", values =c(Items1 = "gold3",Items4 = "#3574b6",Item6 = "#7CA878",Items7 = "red3"),
                        labels = c("Items 1-3", "Items 4-5", "Item 6", "Items 7-10")) +
  scale_alpha_manual(name = "", values =c(Items1 = 0.4,Items4 = 0.9,Item6 = 0.75,Items7 = 1),
                        labels = c("Items 1-3", "Items 4-5", "Item 6", "Items 7-10")) +
  guides(linetype = guide_legend(override.aes = list(size = 1.5)))


#####################
#####################
xxa <- seq(0, 1/2-0.0001, 0.0001)
xxb <- seq(1/2+0.0001, 1, 0.0001)
xx <- c(xxa, xxb)

yy1a <- 0.25 + 0.06 * sin(7 * pi * xxa)
yy2a <- 0.0625 - 0.015 * sin(7 * pi * xxa)

yy1b <- 0.15 - 0.04 * sin(15 * pi * xxb)
yy2b <- 0.065 + 0.25 * ((xxb-0.5)^0.2)
yy3b <- 0.095 + 0.08/6 * sin(15 * pi * xxb)-0.5/6*((xxb-0.5)^0.2)


yy1_2 <- c(yy1a, yy1b)
yy3_4 <- c(yy2a, yy2b)
yy6_10 <- c(yy2a, yy3b)

#2*yy1a+8*yy2a
#2*yy1b+2*yy2b+6*yy3b



#####################

data <- data.frame(xx, yy1_2, yy3_4, yy6_10)

ggplot(data, aes(x = xx)) +
  geom_line(aes(y = yy1_2, linetype = "Items1", color = "Items1", alpha="Items1"), size = 2) +
  geom_line(aes(y = yy3_4, linetype = "Items3", color = "Items3", alpha="Items3"), size = 2) +
  geom_line(aes(y = yy6_10, linetype = "Item6", color = "Item6", alpha="Item6"), size = 2) +
  labs(x = "t", y = expression(paste(pi, "*"))) +
  #theme_minimal() +
  theme_bw() +#
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),#
#        legend.position = c(.05,.935),#
#        legend.box.background = element_rect(color="black"),
legend.text = element_text(family = 'serif',size = 23),#
legend.title=element_blank(),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
	  legend.key.size=unit(5,"lines")) +
  scale_linetype_manual(name = "", values =c(Items1 = "dashed",Items3 = "solid",Item6 = "dotted"),
                        labels = c("Items 1-2", "Items 3-4", "Items 6-10")) +
  scale_color_manual(name = "", values =c(Items1 = "gold3",Items3 = "skyblue",Item6 = "red3"),
                        labels = c("Items 1-2", "Items 3-4", "Items 6-10")) +
  scale_alpha_manual(name = "", values =c(Items1 = 1,Items3 = 1,Item6 = 1),
                        labels = c("Items 1-2", "Items 3-4", "Items 6-10")) +
  guides(linetype = guide_legend(override.aes = list(size = 1.5)))



