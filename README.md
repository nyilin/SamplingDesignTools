SamplingDesignTools: Tools for Dealing with Complex Sampling Designs
================

## Installation

Install the SamplingDesignTools package from GitHub (package devtools
needed):

``` r
install.packages("devtools")
devtools::install_github("nyilin/SamplingDesignTools")
```
## Usage

The SamplingDesignTools package provides tools for working with study designs in
epidemiological research, specifically the following study designs:

**Case-control study:**

- `draw_mcc()` to draw samples with or without matching on confounders.

**Nested case-control (NCC) study:** 

- `compute_km_weights()` to compute Kaplan-Meier (KM) type weights for NCC
samples with or without the full cohort.
- `compute_km_weights_controls()` to compute KM type weights for newly collected
controls when reusing cases in a previously collected NCC sample.

**Counter-matched nested case-control study:**

- `draw_ncc_cm()` to draw counter-matched case-control samples with or without
matching on confounders.

**(More) extreme case-control (MECC) study:**

- `draw_mecc()` to draw MECC samples with or without matching on confounders.
- `analyse_mecc_cond()` to estimate exposure effects from MECC samples.

Help documents for main functions above include minimal examples to illustrate
function usage. Please refer to the 
[online guidebook](https://nyilin.github.io/SamplingDesignTools/) 
for detailed demonstration and discussion.
