# OPLSR
Orthogonal Partial Least Squares Regression (OPLSR)

## Our paper
“Feature Selection Using Distribution of Orthogonal PLS Regression Vectors in Spectral Data”
This paper introduces a useful combined approach of applying orthogonal signal correction (OSC) and permutation tests to PLS for the purpose of feature selection.

## Usage
  *There are three methods for feauture selection: OPLSR, FDR, and Lasso.
  *Algorithms used in experiments can be impledmented by running test_FDR_PCR, test_linear, and testPCRselection_Data.
  *The matters related to the creation and experimentation of simulation data are in ./Simulation.
  *All the datasets, Metabolomics and NIR spectra, are in ./data
  *To produce output, each methods require some useful functions for calculation, call and etc. The functions are in ./utils
