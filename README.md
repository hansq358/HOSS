# Hierarchical Orthogonal Subsampling (HOSS)

## Overview

This repository contains the simulation code, real data analysis and toy example for the HOSS method. 

## Structure

- **Simulation_Code/** - Simulation studies
  - `Case1/` - Continuous multivariate setting
  - `Case2/` - Continuous uniform setting  
  - `Case3/` - Mixed multivariate setting
  - `Case4/` - Mixed uniform setting
  - `Computational efficiency/` - Time complexity analysis

- **Real_Data_Analysis/** - Real data application (London Airbnb data)

- **Toy_Example/** - Simple illustrative example

## Requirements

- R (version 3.6+)
- Rcpp
- Other dependencies specified in individual scripts

## Usage

Run the R scripts in each folder to reproduce the results. The C++ code (`myoss.cpp`, `ortho_helper`) needs to be compiled using `Rcpp::sourceCpp()`.

## License

MIT

