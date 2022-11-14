## Physically-motivated emulation of GMST variability

This repository provides code in supplement to "Separating internal and externally-forced contributions to global temperature variability using a Bayesian stochastic energy balance framework" (M. Schillinger et al. 2022) accepted in Chaos: An Interdisciplinary Journal of Nonlinear Science. A preprint of the manuscript is available at https://arxiv.org/abs/2206.14573 . 

**Authors**: M. Schillinger, B. Ellerhoff, R. Scheichl, K. Rehfeld
**Responsible for this repository:**  Beatrice Ellerhoff ([@bellerhoff](https://github.com/bellerhoff)),  Maybritt Schillinger ([@m-schilinger](https://github.com/m-schillinger))

---

## Organisation of this repository

| directories    | description             |
| -------------- | ----------------------- |
| `./output/` | contains pre-processed data from stochastic multibox EBM fit to data |
| `./output-spectra/` | contains pre-processed spectra of stochastic multibox EBM fit |
| `./plots/` | empty directory to store created figures |

| scripts       | description     |
| ------------- | --------------- |
| `init.R` | load metadata, required libraries and plotting settings |
| `T2_params.R` | Summary of estimated parameters of 2-box fit, creates Table 2 |
| `F2_hadcrut_demo.R` | Application of workflow to observations, creates Figure 2 |
| `F3-F4_spectra_variance.R` | Spectral analysis and emulation of timescale-depenedent variance, creates Figure 3 and 4 |
| `FS5-FS6-TS3_hadcrut_supp.R` | Supplementary validity check of choice of parameters and sampling of internal variance, creates Figues 5 and 6, as well as Table 3 (in Appendix) |
| `F7_plot_fits.R` | Forced response from 2-box fit for all considered runs, Figure 7 in Appendix  |


| additional files             | description                                                  |
| ---------------------------- | ------------------------------------------------------------ |
| `.gitignore`                 | Information for GIT version control to not add several file extensions to version control (e.g. `*.png`, `*.pdf`) |
| `license.md`/ `license.html` | Licensing information                                        |
| `README.md`                  | General README         

---

## Prerequisites

Running the code in this repository requires the following [R](https://www.r-project.org/) packages:

- `ClimBayes` from https://github.com/paleovar/ClimBayes (latest release v.0.1.1)
- `dplyr`
- `ggplot2` 
- `tibble` 
- `PaleoSpec` 
- `RColorBrewer`
- `zoo`
- `ggpubr`
- `purrr`
- `stringr`
- `latex2exp`
- `tidyr`

---

## Pre-processing

The directory `./output/` contains the pre-processed data, which we obtained from fitting the stochastic two-box energy balance model (EBM) to global mean surface temperature (GMST) data using the [*ClimBayes*](https://github.com/paleovar/ClimBayes) package in R. The target data (Table 1 of submitted manuscript) is available from the data holdings of the Climate Research Programme’s Working, from Schmidt et al., Eby et al. (https://climate.uvic.ca/EMICAR5/participants.html), and Morice et al.. The *ClimBayes* package provides detailed information on how to prepare and fit the target data (see vignette... ). To reproduce our runs, use the `ebm_fit_config.yml` in the respective `./output/` directory and the `ebm_fit()` function from *ClimBayes*. We used `n_boxes=2` and `detrending=2`. 

`ebm_fit(temp_data, 
         forc_data,
         start_year,
         end_year,
         n_boxes,
         config,
         detrending,
         config_file)`

---

## References and Acknowledgments

Please see the data availability and acknowledgment statement of the submitted manuscript (https://arxiv.org/abs/2206.14573) and the ClimBayes package. 

We acknowledge the [R Core team](https://www.R-project.org/) and all package developers of packages used in this study. We thank them for their time and dedication to provide R and the packages to the public. Please see `citation()` for details on the R Core Team and `citation("pkgname")` for details on the developers of individual packages.

Please report bugs to the authors (beatrice-marie.ellerhoff(at)uni-tuebingen.de, maybritt.schillinger(at)stat.math.ethz.ch).

*Beatrice Ellerhoff and Maybritt Schillinger, June 2022*
