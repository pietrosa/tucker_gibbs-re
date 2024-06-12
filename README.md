# tucker_gibbs-re
Gibbs sampler for a Bayesian orthogonal Tucker decomposition, accommodating posterior updates based on compressed tensor data

This repository contains code for the article "A Bayesian Hierarchical Model for Orthogonal Tucker Decomposition with Oblivious Tensor Compression" by [Anonymous].

## Data and omissions

The analysies presented in the main text use the [ORL Database of Faces](https://cam-orl.co.uk/facedatabase.html) and a hippocampus dataset obtained from the [Alzheimer's Disease Neuroimaging Initiative](https://adni.loni.usc.edu/). The former is publicly available, but we provide it (after some reshaping) in this repo. The latter is restricted by data privacy protocols: as a result, no content related to this analysis (aside from training code) can be made available here.

The nature of the proposed analysis involves a high degree of randomness and a large number of posterior draws for numerous scalar/matrix/tensor-variate parameters. We do not provide our exact draws (aside from those in the face analysis used to generate figures and tables) as they would likely not be useful to anyone. Comparable results can be obtained by rerunning the analysis.

## Dependencies

- Training code runs in MATLAB and requires the [Tensor Toolbox package](https://www.tensortoolbox.org/) (v3.5 used here)
- Evaluation code runs in MATLAB (and R for some figures/summary statistics) and again requires the Tensor Toolbox and (for R) the scales package (v1.3.0 used here)

