# Type 2 diabetes, body mass index and cataract across ancestries: a two-sample Mendelian randomization study

Khalid Alsaif and Zhongxue Chen, College of Health Solutions, Arizona State University.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23383154.svg)](https://doi.org/10.5281/zenodo.23383154)

This repository holds the analysis code, the harmonized SNP-level data, every method-specific estimate and the R session information for the manuscript "Does the effect of type 2 diabetes and body mass index on cataract differ between East Asians and Europeans? A cross-ancestry Mendelian randomization study" (submitted to the International Journal of Epidemiology, 2026).

## What the study does

Earlier Mendelian randomization (MR) studies in Biobank Japan concluded that type 2 diabetes (T2D) and body mass index (BMI) raise cataract risk in East Asians but not in Europeans. The European check in those studies relied on a UK Biobank self-reported cataract dataset with 5,045 cases, or on a single-variant analysis. This study re-tests the claim with European-derived instruments in Biobank Japan (24,622 cases) and in three European outcome datasets of increasing size (UK Biobank 5,045 cases; FinnGen R5 26,758 cases; a UK Biobank plus FinnGen meta-analysis with 39,519 cases), runs every two-sample MR estimator in TwoSampleMR::mr_method_list(), tests the difference between ancestries formally, measures how well the European instruments transfer to East Asians, and fits a multivariable model with T2D and BMI together.

## Data

All data are publicly available GWAS summary statistics downloaded from the IEU OpenGWAS database (https://gwas.mrcieu.ac.uk) between 15 and 23 September 2026. No individual-level data were used.

| Role | Trait | Source (OpenGWAS id) | Ancestry | Cases / controls or n |
|---|---|---|---|---|
| Exposure | T2D | Xue 2018 (ebi-a-GCST006867) | European | 61,714 / 593,952 |
| Exposure | BMI | Yengo 2018 (ieu-b-40) | European | n 681,275 |
| Outcome | Cataract | Biobank Japan, Ishigaki 2020 (bbj-a-94) | East Asian | 24,622 / 187,831 |
| Outcome | Self-reported cataract | UK Biobank, Neale lab (ukb-d-20002_1278) | European | 5,045 / 356,096 |
| Outcome | Senile cataract | FinnGen R5 (finn-b-H7_CATARACTSENILE) | European | 26,758 / 189,604 |
| Outcome | Cataract | Sakaue 2021, UK Biobank + FinnGen (ebi-a-GCST90018814) | European | 39,519 / 452,358 |
| Transferability | T2D | Biobank Japan (bbj-a-153) | East Asian | 40,250 / 170,615 |
| Transferability | T2D | Spracklen 2020 (ebi-a-GCST010118) | East Asian | 77,418 / 356,122 |
| Transferability | BMI | Biobank Japan (bbj-a-1) | East Asian | n 158,284 |
| East Asian instruments | T2D | Biobank Japan, Suzuki 2019 (bbj-a-77) | East Asian | 36,614 / 155,150 |

Scripts 01 and 07 also process secondary exposures (systolic and diastolic blood pressure, alcohol, smoking initiation); those results are in the supplementary data of the manuscript.

## Folder map

```
scripts/     R scripts 01 to 11, run in order (see below)
data/        exposure instruments and the per-pair outcome extracts as downloaded (RDS)
results/     everything the scripts produce
  harmonised_all.csv / .rds        harmonized SNP-level data for every exposure-outcome pair
  mr_all_methods.csv / .rds        every estimator in mr_method_list() for every pair
  heterogeneity_Q.csv              Cochran's Q
  egger_intercept.csv              MR-Egger intercept test
  cross_ancestry_test.csv          z-test, Cochran's Q and I2 between ancestries
  power_calculations.csv           power for each outcome (Burgess 2014 approximation)
  ukb_logodds_conversion.csv       UK Biobank linear-scale estimates converted to log-odds
  sessionInfo.txt                  R and package versions used
  BST693_results_workbook_v2.xlsx  Excel summary compiled from the CSV files above
  sensitivity/                     MR-PRESSO, leave-one-out, single-SNP, Steiger, plots
  transferability/                 European instruments in East Asian exposure GWAS
  eas_arm_sensitivity/             like-for-like and East Asian-scale analyses
  eas_instruments/                 East Asian-discovered T2D instruments
  mvmr/                            multivariable MR (T2D and BMI together) and conditional F
  figures/                         manuscript figures (PDF and PNG)
```

## Scripts

| Script | What it does |
|---|---|
| 01_pipeline_all_methods.R | Checks every dataset id against the OpenGWAS `ao` table, downloads instruments and outcome data, harmonizes, runs every method in mr_method_list() plus Cochran's Q and the Egger intercept, saves harmonized data and results |
| 02_sensitivity_main_pairs.R | MR-PRESSO, leave-one-out, single-SNP estimates, scatter, funnel and leave-one-out plots for the eight main pairs |
| 03_steiger_direction.R | Steiger direction test, repeated over a grid of prevalence assumptions |
| 04_cross_ancestry_test.R | Formal East Asian versus European comparison (z-test, Cochran's Q, I2) for every method with a usable standard error |
| 05_transferability.R | IVW slope of the European instruments in the East Asian exposure GWAS |
| 06_power.R | Power for each outcome and conversion of the UK Biobank linear-scale estimates to log-odds |
| 07_eas_arm_sensitivity.R | Checks specific to European instruments in an East Asian outcome: palindromic and proxy SNPs, 1 Mb pruning, restriction to instruments with East Asian effects, East Asian-scale estimate |
| 08_mvmr_t2d_bmi.R | Multivariable MR with T2D and BMI together, on each outcome |
| 09_mvmr_strength.R | Conditional F-statistics and heterogeneity for the multivariable model (MVMR package) |
| 10_eas_instruments.R | T2D on Biobank Japan cataract with East Asian-discovered instruments (Biobank Japan and Spracklen 2020), all methods |
| 11_figures.R | Manuscript figures, drawn from the saved results only |

## How to reproduce

1. Install R (the analyses used R 4.5.1) and the packages TwoSampleMR 0.7.9, ieugwasr 1.1.0, mr.raps 0.4.3, MRPRESSO, MVMR and ggplot2. Versions are recorded in `results/sessionInfo.txt` (written by script 01) and `results/sessionInfo_additional_packages.txt` (MRPRESSO and MVMR, loaded by scripts 02 and 09). MRPRESSO and MVMR are installed from GitHub (`remotes::install_github("rondolab/MR-PRESSO")`, `remotes::install_github("WSpiller/MVMR")`).
2. Place this repository at `~/Desktop/BST693` or change `base_dir` at the top of each script.
3. To re-run from the saved downloads in `data/`, run scripts 01 to 11 in order; each script skips any download that is already saved, so no OpenGWAS access is needed.
4. To re-download everything from OpenGWAS instead, move the `data/` and `results/` folders aside (the scripts re-create them), store your OpenGWAS access token as the variable `OPENGWAS_JWT` in `~/.Renviron` (see https://mrcieu.github.io/ieugwasr) and run scripts 01 to 11 in order. OpenGWAS data may have been updated since September 2026, so numbers can differ slightly from the manuscript.

Running scripts 01 to 10 takes a few minutes apart from the downloads and MR-PRESSO, which runs 10,000 simulations per pair.

## Reporting

The manuscript follows the STROBE-MR reporting guideline; the completed checklist is in its supplementary data.

## Software assistance

Claude (Anthropic) was used to assist with the writing of the manuscript and with the development of the R code in this repository. The analyses were performed by the first author, who reviewed and verified the outputs and takes responsibility for their accuracy.

## Licence

The code is released under the MIT licence (see LICENSE). The result tables and harmonized SNP-level data are derived from publicly available GWAS summary statistics and are shared for reproducibility; please cite the original studies listed in the manuscript's Table 1 and the OpenGWAS database (Elsworth et al. 2020, bioRxiv, doi:10.1101/2020.08.10.244293) when re-using them.

## How to cite

Alsaif K, Chen Z. Analysis code and data for "Does the effect of type 2 diabetes and body mass index on cataract differ between East Asians and Europeans? A cross-ancestry Mendelian randomization study". Zenodo, 2026. https://doi.org/10.5281/zenodo.23383154
