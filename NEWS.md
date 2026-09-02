# MMSpower 0.1.0.9000 (development version)

First working version of the package. Each study question has a matched
**design** function (plan the study) and, where applicable, an **analysis**
function (interpret the data). All functions share the same variance model:
Rogan-Gladen correction for imperfect diagnostics, a Kish design effect for
clustered sampling, and an optional finite-population correction.

## Prevalence and precision

* `estimate_prevalence()` -- prevalence point estimate and confidence
  interval from per-cluster counts (Wald, Clopper-Pearson, or
  Agresti-Coull), with the design effect and ICC estimated from the data.
* `design_precision()` -- total sample size for a target margin of error,
  with three design modes (SRS, fixed cluster size, fixed number of
  clusters).

## Threshold decisions

* `test_threshold()` -- one-sided or two-sided test of whether prevalence
  is above/below a decision threshold, using the null-variance z-statistic.
* `design_threshold()` -- sample size to run `test_threshold()` at a target
  power.

## Group comparisons

* `test_difference()` -- two-group prevalence comparison: z-test plus a
  Wald or Newcombe confidence interval on the difference, with per-group
  design effects.
* `design_difference()` -- per-group sample size to detect a given
  prevalence difference at a target power.

## Time trends

* `test_trend()` -- weighted-least-squares test for a linear trend in
  prevalence over time, with a uniform design-effect inflation.
* `design_trend()` -- per-timepoint sample size (or, in reverse, the power
  a given size achieves) to detect a linear trend.

## Detection and cost

* `design_detection()` -- sample size for a target probability of observing
  at least one positive, or the probability a given size achieves. Note:
  `specificity` is currently accepted but not used -- see the README and
  `?design_detection`.
* `estimate_cost()` -- total study cost from per-sample and per-site costs,
  with a `print()` method.

## Documentation

* `vignette("methods")` -- full derivation of the prevalence and precision
  equations.
* `vignette("methods-and-sources")` -- every equation mapped to its place
  in the MMS-SD workshop lecture material, and the parts that extend
  beyond it.
