###############################
# User-adjustable parameters  #
###############################
# Text / scaling (mirrors simulation--structure/simuPlot.R)
BASE_CEX    <- 1.5   # overall scaling via par(cex=...)
AXIS_CEX    <- 2.5   # axis tick label size
LAB_CEX     <- 3     # axis label size
LEGEND_CEX  <- 2.2   # legend text size
TITLE_CEX   <- 3.0   # (reserved)
YLABEL_LINE <- 0     # line offset for y-axis label
FONT_FAMILY <- "Times"  # set to desired font family (e.g., 'Times', 'Helvetica', 'Arial'); NULL for default

# Line widths
LWD_GROUP1       <- 4   # first group (unique style)
LWD_GROUP_OTHERS <- 5   # other groups / bands
LWD_LEGEND       <- 3   # legend key width

# Output device size
PDF_WIDTH  <- 8
PDF_HEIGHT <- 7.5

# Colors
COL_GROUP1 <- '#3574b6'
COL_GROUP2 <- '#c74d12'
COL_GROUP3 <- 'darkorchid4'

# Sequence resolution (can change for smoother lines)
SEQ_BY <- 0.001   # previously 0.008; set smaller for smoother curves

# Whether to save to PDF (set FALSE if you just want interactive plots)
SAVE_PDF <- TRUE

###############################
# Setup                       #
###############################
op <- par(cex = BASE_CEX, family = FONT_FAMILY)
if (!dir.exists('plt')) dir.create('plt', recursive = TRUE)
setlwd <- LWD_GROUP_OTHERS  # keep legacy variable name used below

###############################
# Plot 1 (simu1)              #
###############################
if (SAVE_PDF) pdf('plt/simu1.pdf', width = PDF_WIDTH, height = PDF_HEIGHT, family = FONT_FAMILY)

xx  <- seq(0, 1, SEQ_BY)
yy1 <- 0.2 + 0.03*sin(6*pi*xx)
yy1p <- yy1 + 0.001
yy1n <- yy1 - 0.001
yy2 <- 0.1 - 0.02*sin(6*pi*xx)
yy2p <- yy2 + 0.001
yy2n <- yy2 - 0.001
yy3  <- (rep(1,length(xx)) - 3*yy1 - 3*yy2)/4
yy3p <- yy3 + 0.001
yy3n <- yy3 - 0.001

plot(1,1,xlim=c(0,1),ylim=c(0,0.3),type='n',
     las=1,xlab='t',ylab='', cex.axis=AXIS_CEX, cex.lab=LAB_CEX, yaxt='n')
# Following mirrors the original style: use the +/- bands instead of central line
lines(xx, yy1p, lty=5, col=COL_GROUP1, lwd=setlwd)
lines(xx, yy1n, lty=5, col=COL_GROUP1, lwd=setlwd)
lines(xx, yy2,  lty=1, col=COL_GROUP2, lwd=setlwd)
lines(xx, yy2p, lty=1, col=COL_GROUP2, lwd=setlwd)
lines(xx, yy2n, lty=1, col=COL_GROUP2, lwd=setlwd)
lines(xx, yy3p, lty=3, col=COL_GROUP3, lwd=setlwd)
lines(xx, yy3n, lty=3, col=COL_GROUP3, lwd=setlwd)
title(ylab=expression(paste(pi,"*")), line=YLABEL_LINE, cex.lab=LAB_CEX)

legend(x='topleft',
       legend=c('Group 1','Group 2','Group 3'),
       lty=c(5,1,3),
       lwd=c(LWD_GROUP1, LWD_GROUP_OTHERS, LWD_GROUP_OTHERS),
       col=c(COL_GROUP1, COL_GROUP2, COL_GROUP3),
       bty='n',
       cex=LEGEND_CEX)
if (SAVE_PDF) dev.off()

###############################
# Plot 2 (simu2)              #
###############################
if (SAVE_PDF) pdf('plt/simu2.pdf', width = PDF_WIDTH, height = PDF_HEIGHT, family = FONT_FAMILY)

pert <- 0.001
xx  <- seq(0,1,SEQ_BY)
yy1 <- 0.19 + 0.05*sin(3*pi*xx)
yy1p <- yy1 + pert
yy1n <- yy1 - pert
yy2 <- 0.01 + 0.06*atan(pi*xx)
yy2p <- yy2 + pert
yy2n <- yy2 - pert
yy3p <- (rep(1,length(xx)) - 3*yy1 - 3*yy2)/4 + pert
yy3n <- (rep(1,length(xx)) - 3*yy1 - 3*yy2)/4 - pert

plot(1,1,xlim=c(0,1),ylim=c(0,0.32),type='n',
     las=1,xlab='t',ylab='', cex.axis=AXIS_CEX, cex.lab=LAB_CEX, yaxt='n')
lines(xx, yy1,  lty=5, col=COL_GROUP1, lwd=setlwd)
lines(xx, yy1p, lty=5, col=COL_GROUP1, lwd=setlwd)
lines(xx, yy1n, lty=5, col=COL_GROUP1, lwd=setlwd)
lines(xx, yy2,  lty=1, col=COL_GROUP2, lwd=setlwd)
lines(xx, yy2p, lty=1, col=COL_GROUP2, lwd=setlwd)
lines(xx, yy2n, lty=1, col=COL_GROUP2, lwd=setlwd)
lines(xx, yy3p, lty=3, col=COL_GROUP3, lwd=setlwd)
lines(xx, yy3n, lty=3, col=COL_GROUP3, lwd=setlwd)
title(ylab=expression(paste(pi,"*")), line=YLABEL_LINE, cex.lab=LAB_CEX)

legend(x='topleft',
       legend=c('Group 1','Group 2','Group 3'),
       lty=c(5,1,3),
       lwd=c(LWD_GROUP1, LWD_GROUP_OTHERS, LWD_GROUP_OTHERS),
       col=c(COL_GROUP1, COL_GROUP2, COL_GROUP3),
       bty='n',
       cex=LEGEND_CEX)
if (SAVE_PDF) dev.off()

# Restore original par if needed later
par(op)