<!-- README.md is generated from README.Rmd. Please edit that file -->

# MMSpower

<!-- badges: start -->

[![R-CMD-check](https://github.com/mrc-ide/MMSpower/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mrc-ide/MMSpower/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mrc-ide/MMSpower/graph/badge.svg)](https://app.codecov.io/gh/mrc-ide/MMSpower)
<!-- badges: end -->

**Sample-size, power, and analysis tools for malaria molecular
surveillance (MMS) studies.**

MMSpower is the R companion to the [MMS-SD Study Design
Workshop](https://mrc-ide.github.io/MMS-SD_workshop/). It turns the
workshop’s methods into functions you can run. Each study question comes
as a matched pair: one function to **design** the study (how many
samples, sites, or rounds you need) and one to **analyse** the data once
it is collected.

Every function accounts for:

- **imperfect diagnostics** – Rogan-Gladen correction for test
  sensitivity and specificity;
- **clustered sampling** – a Kish design effect from an intra-cluster
  correlation (ICC), either estimated from the data or supplied;
- **finite populations** – an optional finite-population correction.

## Installation

MMSpower is not on CRAN. Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("mrc-ide/MMSpower")
```

## Functions

| Study question | Design | Analysis |
|----|----|----|
| What is the prevalence, and how precise is the estimate? | `design_precision()` | `estimate_prevalence()` |
| Is prevalence above (or below) a decision threshold? | `design_threshold()` | `test_threshold()` |
| Does prevalence differ between two groups or sites? | `design_difference()` | `test_difference()` |
| Is prevalence changing over time? | `design_trend()` | `test_trend()` |
| Will the study detect a rare variant at all? | `design_detection()` | – |
| What will the study cost? | `estimate_cost()` | – |

## Example

Estimate prevalence from a five-site survey run with an imperfect assay,
adjusting for variation between sites:

``` r
library(MMSpower)

res <- estimate_prevalence(
  x = c(2, 14, 4, 18, 6),      # positives per site
  n = c(60, 60, 55, 65, 50),   # samples per site
  sensitivity = 0.90,
  specificity = 0.98
)

res[c("prevalence", "ci_lower", "ci_upper", "deff", "n_eff")]
#> $prevalence
#> [1] 0.1496865
#>
#> $ci_lower
#> [1] 0.04623187
#>
#> $ci_upper
#> [1] 0.2531412
#>
#> $deff
#> [1] 4.861557
#>
#> $n_eff
#> [1] 59.65167
```

The sites disagree more than sampling alone would explain, so the design
effect (`deff`) is well above 1 and the 290 samples carry only about 60
independent observations’ worth of information (`n_eff`).

## Learn more

- `vignette("methods", package = "MMSpower")` – full derivation of the
  prevalence and precision equations.
- `vignette("methods-and-sources", package = "MMSpower")` – every equation
  mapped to its place in the MMS-SD workshop lecture material, and the
  parts that extend beyond it.
