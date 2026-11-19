# MMSpower

MMSpower provides simple sample size, power and analysis calculations
for malaria molecular surveillance (MMS) studies. It allows for
imperfect diagnostic tests (sensitivity and specificity) and for
clustered sampling across sites (intra-cluster correlation and the
design effect).

## What it covers

Each common study question has a **design** function, which tells you
how many samples you need, and a matching **analysis** function, which
tells you what your data show.

| Question | Design | Analysis |
|----|----|----|
| How precisely can I estimate prevalence? | [`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md) | [`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md) |
| Is prevalence above (or below) a threshold? | [`design_threshold()`](https://mrc-ide.github.io/MMSpower/reference/design_threshold.md) | [`test_threshold()`](https://mrc-ide.github.io/MMSpower/reference/test_threshold.md) |
| Does prevalence differ between two groups? | [`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md) | [`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md) |
| What will the study cost? | [`estimate_cost()`](https://mrc-ide.github.io/MMSpower/reference/estimate_cost.md) |  |

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mrc-ide/MMSpower")
```

## Where next

- **[Get
  started](https://mrc-ide.github.io/MMSpower/articles/MMSpower.md)**: a
  short tour of the package.
- **Tutorials**: worked examples for each type of study (coming soon).
- **[Reference](https://mrc-ide.github.io/MMSpower/reference/index.md)**:
  full documentation for every function, including the methods and
  formulas behind each calculation.

## Development

MMSpower is developed by Sequoia Leuba and Bob Verity at Imperial
College London. The latest in-progress work is on the `develop` branch
of the [GitHub repository](https://github.com/mrc-ide/MMSpower), where
you can also [report
issues](https://github.com/mrc-ide/MMSpower/issues).
