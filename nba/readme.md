
1. Rscript NBAmain.R > runout/output-nba.log 2> runout/error-nba.log
2. Rscript nbaPlot.R
3. Rscript PanelPlot.R
4. Rscript noPanelPlot.R



files:
- NBAmain.R: main for NBA analysis
- nbaPlot.R: plot the heatmap
- PanelPlot.R: plot the score estimation with panelty
- noPanelPlot.R: plot the score estimation without panelty
- PFunc.R, Func.R, FuncContinue.R: functions