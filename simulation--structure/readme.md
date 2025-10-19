This folder is for the structure detection method simulation part:
procedure to reproduce the results:

1. Rscript Para1.R > runout/output-Para1.log 2> runout/error-Para1.log
2. Rscript Para2.R > runout/output-Para2.log 2> runout/error-Para2.log
3. Rscript simuPlot.R
4. Rscript parameter_sensitivity.R > runout/output-sen1.log 2> runout/error-sen1.log
5. Rscript parameter_sensitivity_para2.R > runout/output-sen2.log 2> runout/error-sen2.log




files
- Para1.R, Para2.R: main files for parameter settings 1 and 2
- parameter_sensitivity.R, parameter_sensitivity_para2.R: main files for sensitivity analysis
- dpFunc.R, dpFuncContinue.R: functions
- gpFunc.R: functions for grouping
- simuPlot.R: plots