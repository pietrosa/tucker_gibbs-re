# tucker_gibbs-re
Gibbs sampler for a Bayesian orthogonal Tucker decomposition, accommodating posterior updates based on compressed tensor data

This repository contains code for the article "A Bayesian Hierarchical Model for Orthogonal Tucker Decomposition with Oblivious Tensor Compression" by [Anonymous].

## Data and omissions

The analysies presented in the main text use the [ORL Database of Faces](https://cam-orl.co.uk/facedatabase.html) and a hippocampus dataset obtained from the [Alzheimer's Disease Neuroimaging Initiative](https://adni.loni.usc.edu/). The former is publicly available, but we provide it (after some reshaping, in `/analysis/face_dat.mat`) in this repo. The latter is restricted by data privacy protocols: as a result, no content related to this analysis (aside from training code) can be made available here.

The nature of the proposed analysis involves a high degree of randomness and a large number of posterior draws for numerous scalar/matrix/tensor-variate parameters. We do not provide our exact draws (aside from those in the face analysis used to generate figures and tables) as they would likely not be useful to anyone. Comparable results can be obtained by rerunning the analysis.

## Dependencies

- Training code runs in MATLAB and requires the [Tensor Toolbox package](https://www.tensortoolbox.org/) (v3.5 used here)
- Evaluation code runs in MATLAB (and R for some figures/summary statistics) and again requires the Tensor Toolbox and (for R) the scales package (v1.3.0 used here)

## Training

Training code is in `analysis/analysis_face.m` and `analysis/analysis_hippo.m`. Both procedures are similar. As mentioned above, hipppocampus data is not provided in this repo.

Both files provide the code used in the main article. For a quick demo, consider the face analysis with:
```
R_list = {[5 5 5]};
m_list = {size(X), [50 60 240]};
ndraws = 22;
thin = 2;
burnin = 2;
nsave = floor((ndraws-burnin)/thin);
nchains = 2;
```
This yields the following (random) table of saved draws for Frobenius reconstruction error, posterior draw time, and tau2. Column names appear `nchains` times, corresponding to each parallel chain. The prefixes `m#` and `R#` refer to indices of `R_list` and `m_list`, respectively, e.g., `m2` refers to `m=[50 60 240]`.

|R1_m1_err|R1_m1_err|R1_m1_time|R1_m1_time|R1_m1_tau2|R1_m1_tau2|R1_m2_err|R1_m2_err|R1_m2_time|R1_m2_time|R1_m2_tau2|R1_m2_tau2|
|---------|---------|----------|----------|----------|----------|---------|---------|----------|----------|----------|----------|
|249.3189 |239.859  |0.099951  |0.099225  |66.23076  |71.69608  |245.0925 |244.0667 |0.105617  |0.104867  |74.80457  |72.13396  |
|239.6087 |237.8957 |0.123427  |0.121331  |71.76725  |72.85041  |269.8845 |245.6298 |0.127229  |0.120984  |62.82951  |71.15046  |
|237.7802 |237.5639 |0.11339   |0.114968  |72.8204   |73.01592  |241.3461 |239.99   |0.063112  |0.071301  |80.43164  |77.9691   |
|237.5609 |237.622  |0.095804  |0.113694  |73.01824  |73.10432  |243.4717 |240.4419 |0.064766  |0.068775  |67.92567  |67.36817  |
|237.59   |237.5966 |0.112479  |0.091596  |72.92893  |72.91888  |239.5319 |239.6682 |0.0691    |0.070683  |74.64206  |71.72057  |

The `draws` output can be used to further explore draws for any model parameter or for posterior inference.

## Evaluation

Running `analysis/face_reconstruction.m` will fit new models and generate (random) face reconstructions. The value of `kk` can be tweaked to select the saved draw to use in the reconstruction. Taking `kk=nsave` (i.e., using the last draw) yields the following reconstructions (of 16 random faces):

![untitled](https://github.com/pietrosa/tucker_gibbs-re/assets/40504922/59bcfd55-14c4-4fb1-8cce-a1bddfae01b6)

For quantitative results, see `analysis/face_results.R`, which uses the results in `analysis/faces_chains.csv`. For use with other results, the "Load Results" block should be modified. The code yields the following chain plot, convergence checks (i.e., that $\hat R\in [0.9934088, 1.0042507]$ and $n_\text{eff}\in[47.12596,133.33333]$), the following relative error/time plot, and the below tabular summary.

![Rplot](https://github.com/pietrosa/tucker_gibbs-re/assets/40504922/ef04a74c-5ce5-443f-a5ea-b0c8a5e0e30b)

![Rplot02](https://github.com/pietrosa/tucker_gibbs-re/assets/40504922/da1620b4-90f7-42b2-8348-8afb3ffc073c)

|R  |DR  |quant|Rhat             |n_eff           |err          |time       |
|---|----|-----|-----------------|----------------|-------------|-----------|
|R5 |DR1 |err  |0.998115377470114|133.333333333333|237.6 (0.9)  |0.13 (0.02)|
|R5 |DR1 |tau2 |0.998587389674479|132.428226241168|237.6 (0.9)  |0.13 (0.02)|
|R5 |DR.8|err  |1.00314456501515 |58.7483181305308|238.2 (2.1)  |0.16 (0.01)|
|R5 |DR.8|tau2 |0.996451726772937|133.333333333333|238.2 (2.1)  |0.16 (0.01)|
|R5 |DR.6|err  |1.00185486820677 |133.333333333333|239.4 (2.2)  |0.10 (0.01)|
|R5 |DR.6|tau2 |0.998017583118948|129.23992459734 |239.4 (2.2)  |0.10 (0.01)|
|R5 |DR.4|err  |1.00191280395579 |133.333333333333|241.5 (2.1)  |0.07 (0.01)|
|R5 |DR.4|tau2 |0.995817749666305|133.333333333333|241.5 (2.1)  |0.07 (0.01)|
|R5 |DR.2|err  |0.996463124993707|133.333333333333|251.7 (5.2)  |0.06 (0.01)|
|R5 |DR.2|tau2 |0.993446461722014|133.333333333333|251.7 (5.2)  |0.06 (0.01)|
|R15|DR1 |err  |0.998425166863652|108.20340482605 |190.1 (17.3) |0.32 (0.02)|
|R15|DR1 |tau2 |1.00069517208229 |115.898582903817|190.1 (17.3) |0.32 (0.02)|
|R15|DR.8|err  |0.999165793042363|133.333333333333|192.4 (23.8) |0.32 (0.01)|
|R15|DR.8|tau2 |0.998740791891143|133.333333333333|192.4 (23.8) |0.32 (0.01)|
|R15|DR.6|err  |0.998384576064762|133.333333333333|195.5 (18.9) |0.30 (0.02)|
|R15|DR.6|tau2 |0.998963523742757|133.333333333333|195.5 (18.9) |0.30 (0.02)|
|R15|DR.4|err  |1.00073360488472 |133.333333333333|204.4 (21.8) |0.27 (0.02)|
|R15|DR.4|tau2 |1.00211201556454 |133.333333333333|204.4 (21.8) |0.27 (0.02)|
|R15|DR.2|err  |0.99943109108091 |47.1259636167842|367.4 (75.0) |0.25 (0.01)|
|R15|DR.2|tau2 |1.00089833498435 |133.333333333333|367.4 (75.0) |0.25 (0.01)|
|R30|DR1 |err  |0.993408757838817|133.333333333333|179.9 (54.8) |0.93 (0.03)|
|R30|DR1 |tau2 |0.994853697823012|133.333333333333|179.9 (54.8) |0.93 (0.03)|
|R30|DR.8|err  |1.00235300054161 |92.1230961388292|192.6 (71.8) |0.92 (0.02)|
|R30|DR.8|tau2 |1.00013448010199 |133.333333333333|192.6 (71.8) |0.92 (0.02)|
|R30|DR.6|err  |0.996651088346188|133.333333333333|205.9 (75.9) |0.87 (0.02)|
|R30|DR.6|tau2 |1.00425068667868 |133.333333333333|205.9 (75.9) |0.87 (0.02)|
|R30|DR.4|err  |1.00346414317305 |133.333333333333|296.9 (93.8) |0.83 (0.03)|
|R30|DR.4|tau2 |0.996114471399055|133.333333333333|296.9 (93.8) |0.83 (0.03)|
|R30|DR.2|err  |1.00047286113395 |121.180822207341|1073.5 (84.0)|0.78 (0.01)|
|R30|DR.2|tau2 |0.996184071293635|126.229125232   |1073.5 (84.0)|0.78 (0.01)|
