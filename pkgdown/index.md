# MMSpower

MMSpower provides simple sample size, power and analysis calculations for
malaria molecular surveillance (MMS) studies. It allows for imperfect
diagnostic tests (sensitivity and specificity) and for clustered sampling
across sites (intra-cluster correlation and the design effect).

## What it covers

Each common study question has a **design** function, which tells you how
many samples you need, and a matching **analysis** function, which tells you
what your data show.

| Question | Design | Analysis |
|---|---|---|
| How precisely can I estimate prevalence? | [`design_precision()`](reference/design_precision.html) | [`estimate_prevalence()`](reference/estimate_prevalence.html) |
| Is prevalence above (or below) a threshold? | [`design_threshold()`](reference/design_threshold.html) | [`test_threshold()`](reference/test_threshold.html) |
| Is a variant present at all? | [`design_detection()`](reference/design_detection.html) | |
| Does prevalence differ between two groups? | [`design_difference()`](reference/design_difference.html) | [`test_difference()`](reference/test_difference.html) |
| Is prevalence changing over time? | [`design_trend()`](reference/design_trend.html) | [`test_trend()`](reference/test_trend.html) |
| What will the study cost? | [`estimate_cost()`](reference/estimate_cost.html) | |

## Installation

```r
# install.packages("remotes")
remotes::install_github("mrc-ide/MMSpower")
```

## Where next

- **[Get started](articles/MMSpower.html)**: a short tour of the package.
- **Tutorials**: worked examples for each type of study (coming soon).
- **[Reference](reference/index.html)**: full documentation for every
  function, including the methods and formulas behind each calculation.

## Development

MMSpower is developed by Sequoia Leuba and Bob Verity at Imperial College
London. The latest in-progress work is on the `develop` branch of the
[GitHub repository](https://github.com/mrc-ide/MMSpower), where you can also
[report issues](https://github.com/mrc-ide/MMSpower/issues).
