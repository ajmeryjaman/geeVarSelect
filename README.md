# geeVarSelect
#### This repository contains code to perform variable selection when analyzing correlated continuous outcomes with generalized estimating equations (GEE) using the boosting algorithm or the penalization approach. The goal is to reproduce one key result of the paper ["EEBoost: A General Method for Prediction and Variable Selection Based on Estimating Equations"](https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10098) by Wolfson (2011, JASA).

The repository contains the following folders:

Folder | Description
--- | ---
src | This folder includes all ".R" files containing the source codes.
data | This folder includes "seeds.Rdata" file which is used in simulations.
docs | This folder includes all ".Rd" format documentation files.

The src folder contains the following files:

File | Description
--- | ---
[geeFUNCTIONS.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/geeFUNCTIONS.R) | This file contains all the required functions to perform penalized GEE.
[penalizedGEE.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/penalizedGEE.R) | This file contains the main function to perform penalized GEE.
[geeBOOST.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/geeBOOST.R) | This file contains the functions to perform GEE boosting algorithm.
[simulation.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/simulation.R) | This file contains the function that generates a correlated dataset, performs estimation using the two algorithms, and evaluates the prediction errors on different test datasets.
[runSim_all.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/runSim_all.R) | This file runs independent simulations under different setting and saves the results.
[table2.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/table2.R) | This file loads the saved results, calculates the performance metrics and creates the output table.

To reproduce Table 2 of the JASA paper, first we need to run the file [runSim_all.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/runSim_all.R); this will create four ".Rdata" format files containing the results. Execution of [runSim_all.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/runSim_all.R) requires time. To create the final output, i.e. Table  2, we need to run [table2.R](https://github.com/ajmeryjaman/geeVarSelect/blob/main/src/table2.R).



