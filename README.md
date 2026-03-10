# Simulation Study: Dynamic Borrowing and Prognostic Covariate Adjustment for RCTs With External Controls

## Overview
This repository contains the code used to conduct a simulation study evaluating statistical methods for incorporating historical control data in randomized clinical trials. 

The simulations investigate the operating characteristics of several estimators under varying levels of compatibility between historical and current trial data. Compatibility is explored through scenarios involving mean shifts, covariate distribution shifts, and changes in the outcome-generating mechanism.


## Repository Structure
├── data_generating/
│ Functions defining the data-generating mechanisms used in the simulations.
│
├── analysis/
│ Implementation of the statistical methods evaluated in the study.
│ Includes both frequentist and Bayesian approaches.
│
├── simulation/
│ Scripts used to run the simulation scenarios and simulation replicates.
│ simulation.R contains a self-contained simulation script where the
│ configuration settings are specified directly within the file.
│ simulation_0.R provides a minimal simulation runner where the
│ configuration is supplied externally when running from the terminal, e.g.
│
│ Rscript simulation/simulation_0.R config/input_baseline.R
│
├── plots/
│ Code used to generate figures and visual summaries of the simulation results.
│
├── results/
│ Saved outputs from simulation runs.
│
├── utils/
│ Helper functions used across scripts and loading packages. 
│
└── README.md

