# Hierarchical Orthogonal Subsampling (HOSS)

## Overview

This repository contains the simulation code, real data analysis and toy example for the HOSS method. 

## Structure

- **Simulation_Code/** - Simulation studies
  - `Case1/` - Continuous multivariate setting
  - `Case2/` - Continuous uniform setting  
  - `Case3/` - Mixed multivariate setting
  - `Case4/` - Mixed uniform setting
  - `Additional_Experiments/` - Predictive performance (MSPE) and variance components, encoding sensitivity, model misspecification, sequential vs joint design, sample-size scaling, runtime vs feature dimension, and the balanced subsampling comparison. Shared logic sits in `common.R`, and each experiment has its own runner script.

- **Real_Data_Analysis/** - Real data application (London Airbnb data)

- **Toy_Example/** - Simple illustrative example

## Requirements

- R 3.6 or later
- MASS
- AlgDesign
- Rcpp
- RcppArmadillo
- A C++11-compatible compiler

## Usage

Run the R scripts in each folder to reproduce the results. The C++ code (`myoss.cpp`) needs to be compiled using `Rcpp::sourceCpp()`. The scripts in `Additional_Experiments/` do this automatically.

## License

MIT
