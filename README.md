# RankStruct

This repository contains the code accompanying the paper **"Preference-Based Dynamic Ranking Structure Recognition"**.  
It includes three main modules with reproducible analyses and plotting utilities:

- `nba/` — NBA data analysis and visualization  
- `simulation--grouping/` — grouping simulation experiments  
- `simulation--structure/` — structure-detection simulation and sensitivity analyses

---
## Reproducing results

### NBA analysis (`nba/`)

Reproduce the NBA analysis and plotting:

1. Run the main analysis (saves logs under `runout/`):

```bash
Rscript nba/NBAmain.R > runout/output-nba.log 2> runout/error-nba.log
```

2. Create plots:

```bash
Rscript nba/nbaPlot.R
Rscript nba/PanelPlot.R
Rscript nba/noPanelPlot.R
```

Files of interest:
- `NBAmain.R` — main script for NBA analysis
- `nbaPlot.R` — heatmap plotting
- `PanelPlot.R` — score estimation with penalty
- `noPanelPlot.R` — score estimation without penalty
- `PFunc.R`, `Func.R`, `FuncContinue.R` — supporting functions

### Grouping simulation (`simulation--grouping/`)

Reproduce the grouping simulation results:

1. Run the simulations (logs to `runout/`):

```bash
Rscript simulation--grouping/simu1.R > runout/output-simu1.log 2> runout/error-simu1.log
Rscript simulation--grouping/simu2.R > runout/output-simu2.log 2> runout/error-simu2.log
```

2. Generate plots:

```bash
Rscript simulation--grouping/simuPlot.R
```

Files of interest:
- `simu1.R`, `simu2.R` — simulation for different parameter sets
- `simuFunc.R` — helper functions used by the simulations
- `simuPlot.R` — script to produce plots

### Structure detection simulation (`simulation--structure/`)

Reproduce the structure-detection simulation and sensitivity analyses:

1. Run the main parameter experiments (logs to `runout/`):

```bash
Rscript simulation--structure/Para1.R > runout/output-Para1.log 2> runout/error-Para1.log
Rscript simulation--structure/Para2.R > runout/output-Para2.log 2> runout/error-Para2.log
```

2. Generate plots:

```bash
Rscript simulation--structure/simuPlot.R
```

3.  Sensitivity analysis:

```bash
Rscript simulation--structure/parameter_sensitivity.R > runout/output-sen1.log 2> runout/error-sen1.log
Rscript simulation--structure/parameter_sensitivity_para2.R > runout/output-sen2.log 2> runout/error-sen2.log
```


Files of interest:
- `Para1.R`, `Para2.R` — main parameter-setting scripts
- `parameter_sensitivity.R`, `parameter_sensitivity_para2.R` — sensitivity analyses
- `dpFunc.R`, `dpFuncContinue.R` — core dynamic programming functions
- `gpFunc.R` — grouping-related functions
- `simuPlot.R` — plotting utilities
