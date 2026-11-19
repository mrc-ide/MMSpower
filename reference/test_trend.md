# Test whether prevalence is changing linearly over time

Given the number of positives and the number tested at several
timepoints (optionally at several sites per timepoint), fits a straight
line to the measured prevalences and tests whether its slope differs
from zero with the Cochran-Armitage trend test. Returns the estimated
change in true prevalence per unit of time with its confidence interval,
the p-value and whether the trend is significant. It is the analysis
twin of
[`design_trend()`](https://mrc-ide.github.io/MMSpower/reference/design_trend.md).

**Hypotheses** (\\\beta\\ is the slope of true prevalence):

- `"two.sided"` (default):

  H0: \\\beta = 0\\ against H1: \\\beta \ne 0\\.

- `"greater"`:

  H0: \\\beta \le 0\\ against H1: \\\beta \> 0\\ (an increasing trend).

- `"less"`:

  H0: \\\beta \ge 0\\ against H1: \\\beta \< 0\\ (a decreasing trend).

## Usage

``` r
test_trend(
  x,
  n,
  time = NULL,
  alternative = c("two.sided", "greater", "less"),
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  icc = NULL
)
```

## Arguments

- x:

  Vector of whole-number counts of people who tested positive: one per
  timepoint, or one per site and timepoint (give the timepoint of each
  in `time`).

- n:

  Vector of whole-number totals tested, the same length as `x`.

- time:

  Numeric vector giving the timepoint of each value of `x` (for example
  calendar years or survey rounds). If `NULL` (default), `1, 2, ...` is
  used, which only makes sense with one value per timepoint.

- alternative:

  Direction of the test: `"two.sided"` (default), `"greater"`
  (increasing) or `"less"` (decreasing).

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1 (the test never misses a
  case).

- specificity:

  Diagnostic specificity in (0, 1\]; default 1 (the test never gives a
  false positive). `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level for the interval, in (0, 1); default 0.95. The
  significance level of the test is 1 - `conf_level`.

- icc:

  Optional. Intra-cluster correlation, in \[0, 1\]. If `NULL` (default),
  it is estimated from the timepoints with two or more sites. Set to `0`
  for no clustering.

## Value

A list with:

- slope:

  Estimated change in true prevalence per unit of time.

- ci_lower, ci_upper:

  Confidence interval for the slope; `ci_upper` is `Inf` for `"greater"`
  and `ci_lower` is `-Inf` for `"less"`.

- p_value:

  p-value for the chosen `alternative`.

- reject:

  `TRUE` if `p_value` is below 1 - `conf_level`.

- n_eff:

  Effective sample size, `sum(n) / deff`.

- icc:

  ICC used: as supplied, or estimated from the data.

- deff:

  Design effect used (1 with no clustering).

## Details

**Line.** The counts are added up within each timepoint, giving \\X_j\\
positives of \\N_j\\ tested at time \\t_j\\, and measured prevalence
\\\hat p_j = X_j / N_j\\. Each timepoint is weighted by its effective
number tested, \\w_j = N_j / D\_{eff}\\. The slope is
\$\$\hat\beta\_{app} = \frac{\sum_j w_j (t_j - \bar t)\\\hat p_j}
{\sum_j w_j (t_j - \bar t)^2},\$\$ where \\\bar t\\ is the weighted mean
time.

**Test.** With no trend every timepoint has the weighted mean prevalence
\\\bar p\\, and the slope has variance \\\bar p (1 - \bar p) / \sum_j
w_j (t_j - \bar t)^2\\. The z-statistic is the slope divided by its
standard error (the Cochran-Armitage trend test). Without sites it gives
the same p-value as R's
[`prop.trend.test()`](https://rdrr.io/r/stats/prop.trend.test.html).

**Interval.** For a possible slope \\\beta\\, take the line with that
slope through \\\bar p\\ at \\\bar t\\, cut to \[0, 1\]; call its value
at each timepoint \\f_j\\. If that line were true, the slope would have
variance \\V(\beta) = \sum_j w_j (t_j - \bar t)^2 f_j (1 - f_j) /
(\sum_j w_j (t_j - \bar t)^2)^2\\. The interval holds the slopes with
\\(\hat\beta\_{app} - \beta)^2 \le z^2 V(\beta)\\. At \\\beta = 0\\ this
is the test, so `reject` is `TRUE` exactly when 0 is outside the
interval. It is matched to `alternative`: one-sided for `"greater"` and
`"less"`. The slope and its interval are divided by \\Se + Sp - 1\\ to
give the trend in true prevalence.

**Design effect.** It is \\1 + (n_w - 1) ICC\\, where \\n_w = \sum n^2 /
\sum n\\ is the mean number of people per site, weighted by site size,
and is applied to every timepoint. With `icc = NULL`, an ICC is
estimated within each timepoint that has two or more sites (as in
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)),
and the estimates are averaged; timepoints with one person per site, or
with all negative or all positive, are skipped. With one site per
timepoint there is nothing to estimate from, and a supplied `icc` is
ignored with a warning.

## References

Armitage, P. (1955) Tests for linear trends in proportions and
frequencies. *Biometrics* **11**(3), 375-386.
[doi:10.2307/3001775](https://doi.org/10.2307/3001775)

Cochran, W. G. (1954) Some methods for strengthening the common
chi-square tests. *Biometrics* **10**(4), 417-451.
[doi:10.2307/3001616](https://doi.org/10.2307/3001616)

Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)

Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
of a screening test. *American Journal of Epidemiology* **107**(1),
71-76.
[doi:10.1093/oxfordjournals.aje.a112510](https://doi.org/10.1093/oxfordjournals.aje.a112510)

MMS-SD Study Design Workshop.
<https://mrc-ide.github.io/MMS-SD_workshop/>

## See also

[`design_trend()`](https://mrc-ide.github.io/MMSpower/reference/design_trend.md)

## Examples

``` r
# Five annual rounds, one survey per round
test_trend(x = c(5, 8, 9, 14, 17), n = rep(100, 5), time = 2019:2023)
#> $slope
#> [1] 0.03
#> 
#> $ci_lower
#> [1] 0.01096155
#> 
#> $ci_upper
#> [1] 0.0482649
#> 
#> $p_value
#> [1] 0.00205777
#> 
#> $reject
#> [1] TRUE
#> 
#> $n_eff
#> [1] 500
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 

# Two sites per round, ICC estimated from the data
test_trend(
  x    = c(2, 3, 4, 5, 6, 8, 7, 10),
  n    = c(50, 45, 55, 50, 48, 52, 47, 51),
  time = c(1, 1, 2, 2, 3, 3, 4, 4)
)
#> $slope
#> [1] 0.04175961
#> 
#> $ci_lower
#> [1] 0.01363743
#> 
#> $ci_upper
#> [1] 0.06858299
#> 
#> $p_value
#> [1] 0.003665485
#> 
#> $reject
#> [1] TRUE
#> 
#> $n_eff
#> [1] 398
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 

# Imperfect diagnostic test, one-sided (a decline)
test_trend(x = c(20, 15, 11, 7), n = rep(120, 4), time = 1:4,
           alternative = "less", sensitivity = 0.9, specificity = 0.97)
#> $slope
#> [1] -0.04118774
#> 
#> $ci_lower
#> [1] -Inf
#> 
#> $ci_upper
#> [1] -0.01705297
#> 
#> $p_value
#> [1] 0.002550328
#> 
#> $reject
#> [1] TRUE
#> 
#> $n_eff
#> [1] 480
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 
```
