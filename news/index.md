# Changelog

## MMSpower 1.0.0

First release.

- Design functions for four common questions in malaria molecular
  surveillance:
  [`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
  (estimating a prevalence),
  [`design_threshold()`](https://mrc-ide.github.io/MMSpower/reference/design_threshold.md)
  (testing against a threshold),
  [`design_detection()`](https://mrc-ide.github.io/MMSpower/reference/design_detection.md)
  (detecting a rare variant) and
  [`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md)
  (comparing two groups). Each returns the sample size needed for a
  target, or, given a full design, the margin of error, power or
  detection probability that design achieves.
- Matching analysis functions:
  [`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md),
  [`test_threshold()`](https://mrc-ide.github.io/MMSpower/reference/test_threshold.md)
  and
  [`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md).
- [`estimate_cost()`](https://mrc-ide.github.io/MMSpower/reference/estimate_cost.md)
  costs a design by stage (screened, collected, sequenced) plus a
  per-site cost, and checks it against a budget.
- All functions allow for an imperfect diagnostic test and for clustered
  sampling across sites, through the intra-cluster correlation (ICC).
