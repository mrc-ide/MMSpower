# MMSpower 1.0.0

First release.

* Design functions for four common questions in malaria molecular
  surveillance: `design_precision()` (estimating a prevalence),
  `design_threshold()` (testing against a threshold), `design_detection()`
  (detecting a rare variant) and `design_difference()` (comparing two
  groups). Each returns the sample size needed for a target, or, given a
  full design, the margin of error, power or detection probability that
  design achieves.
* Matching analysis functions: `estimate_prevalence()`, `test_threshold()`
  and `test_difference()`.
* `estimate_cost()` costs a design by stage (screened, collected,
  sequenced) plus a per-site cost, and checks it against a budget.
* All functions allow for an imperfect diagnostic test and for clustered
  sampling across sites, through the intra-cluster correlation (ICC).
