This folder is for the group method simulation part

procedure to reproduce the results:

1. Rscript simu1.R > runout/output-simu1.log 2> runout/error-simu1.log
2. Rscript simu2.R > runout/output-simu2.log 2> runout/error-simu2.log
3. Rscript simuPlot.R



files
- simu1.R: main code for simulation parameter 1
- simu2.R: main code for simulation parameter 2
- simuFunc.R: functions for simulations
- simuPlot.R: plot for simulation parameters
