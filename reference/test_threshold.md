# Test whether prevalence is above or below a decision threshold

Analysis function: given the number of positives and the number tested
at one or more sites, tests whether the true prevalence is above (or
below) a decision threshold, with an exact binomial test (default) or a
Wald z-test. Returns the p-value, whether the result is significant, and
the confidence interval for the true prevalence that matches the test:
Clopper-Pearson for the exact test, Wald for the Wald test. It is the
analysis twin of
[`design_threshold()`](https://mrc-ide.github.io/MMSpower/reference/design_threshold.md).

With a perfect diagnostic test, the share of people who test positive is
the prevalence, and it is compared with the threshold directly. With an
imperfect test it is not: some people with the condition test negative,
and some without it test positive. So the threshold is converted to the
share who would test positive if the true prevalence were exactly at the
threshold, \\\theta \cdot Se + (1 - \theta)(1 - Sp)\\, where \\\theta\\
is the threshold, \\Se\\ the sensitivity and \\Sp\\ the specificity. For
example, with a threshold of 5 percent, sensitivity 0.9 and specificity
0.95: of the 5 percent who have the condition, 0.9 test positive
(0.045), and of the 95 percent who do not, 0.05 wrongly test positive
(0.0475), so 9.25 percent would test positive. The test then asks
whether the share testing positive is far enough above 9.25 percent that
chance alone is unlikely to explain it: with 100 people tested, 12
positives (12 percent) is not enough (p = 0.21), but 15 (15 percent) is
(p = 0.04). This works because, when \\Se + Sp\\ is above 1, a higher
true prevalence always gives a higher share testing positive. The
prevalence and its interval are then converted back to the true scale
(Rogan-Gladen correction). Clustering is allowed for with a design
effect, as in
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md).

**Hypotheses** (\\p\\ is the true prevalence, \\\theta\\ the threshold):

- `"greater"` (default):

  H0: \\p \le \theta\\ against H1: \\p \> \theta\\: is prevalence above
  the threshold?

- `"less"`:

  H0: \\p \ge \theta\\ against H1: \\p \< \theta\\: is prevalence below
  the threshold?

- `"two.sided"`:

  H0: \\p = \theta\\ against H1: \\p \ne \theta\\. A threshold question
  is usually one-way (above or below), so `"greater"` or `"less"` is
  usually the right choice.

## Usage

``` r
test_threshold(
  x,
  n,
  threshold,
  alternative = c("greater", "less", "two.sided"),
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  icc = NULL,
  method = c("exact", "wald")
)
```

## Arguments

- x:

  Vector of whole-number counts of people who tested positive, one per
  site (a single number for one site).

- n:

  Vector of whole-number totals tested, one per site (same length as
  `x`).

- threshold:

  Decision threshold for the true prevalence, strictly between 0 and 1
  (for example 0.05 for the WHO pfhrp2/3 deletion threshold).

- alternative:

  Direction of the test: `"greater"` (default, prevalence above the
  threshold), `"less"` (below it) or `"two.sided"` (different from it).

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

  Optional. Intra-cluster correlation (ICC), in \[0, 1\]. If `NULL`
  (default), it is estimated from the sites. Set to `0` for no
  clustering. A supplied value is also used with only one site (see
  Details).

- method:

  `"exact"` (default): exact binomial test with the Clopper-Pearson
  interval. `"wald"`: Wald z-test with the Wald interval. Each interval
  is its own test turned around, so `reject` is `TRUE` exactly when
  `threshold` is outside the interval, whichever method is used.

## Value

A list with:

- p_value:

  p-value for the chosen `alternative`.

- reject:

  `TRUE` if `p_value` is below 1 - `conf_level`.

- prevalence:

  Estimated true prevalence (Rogan-Gladen corrected, between 0 and 1).

- ci_lower, ci_upper:

  Confidence interval for the true prevalence; `ci_upper` is 1 for
  `"greater"` and `ci_lower` is 0 for `"less"`.

- n_eff:

  Effective sample size, `sum(n) / deff`.

- icc:

  ICC used: as supplied, or estimated from the data.

- deff:

  Design effect used (1 with no clustering).

- df:

  Degrees of freedom of the estimated design effect, sites - 1; `Inf`
  when `icc` is supplied or not estimated.

## Details

**Design effect** (returned as `deff`). It is \\1 + (n_w - 1) ICC\\,
where \\n_w = \sum n^2 / \sum n\\ is the mean number of people per site,
weighted by site size. With `icc` supplied (for example from earlier
surveys), that value is used. Supply `icc` with one site only when the
site stands for a wider population, such as a district; if the site
itself is the population of interest, leave `icc = NULL`. With
`icc = NULL`, the ICC is estimated from the data as in
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md);
this needs at least two sites, and with one site the design effect is 1.

**Few sites.** An ICC estimated from \\k\\ sites rests on \\k - 1\\
degrees of freedom (returned as `df`), so the design effect is itself
uncertain. The tests allow for this, as
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)
does: the exact test and its interval use the effective sample size
shrunk to \\m (z / t\_{k-1})^2\\ (Korn and Graubard, 1998), with \\z\\
and \\t\_{k-1}\\ the normal and t quantiles for the tail area of the
test (\\\alpha\\, or \\\alpha / 2\\ for two-sided); the Wald test uses a
t distribution with \\k - 1\\ degrees of freedom in place of the normal.
With `icc` supplied, or one site, nothing changes. Because the shrink
depends on \\\alpha\\, the exact p-value then depends a little on
`conf_level`; `reject` and the interval always agree.

**Effective sample size.** \\m\\ (returned as `n_eff`) is the number
tested divided by the design effect. The number of positives is scaled
the same way, to \\\hat{p} m\\, where \\\hat{p} = \sum x / \sum n\\ is
the measured prevalence. Neither need be a whole number: the binomial
tail probabilities are worked out through the beta distribution, which
equals [`pbinom()`](https://rdrr.io/r/stats/Binomial.html) for whole
numbers.

**Exact test** (`method = "exact"`). Suppose the true prevalence is
exactly at the threshold. For `"greater"`, the p-value is the chance of
at least as many positives as observed; for `"less"`, the chance of at
most as many. The threshold is first converted to the measured scale,
\\\theta\_{app} = \theta \cdot Se + (1 - \theta)(1 - Sp)\\. With one
site and a design effect of 1, the one-sided p-values are exactly those
of R's [`binom.test()`](https://rdrr.io/r/stats/binom.test.html).

**Wald test** (`method = "wald"`). The normal approximation,

\$\$z = \frac{\hat{p} - \theta\_{app}}{\sqrt{\hat{p}(1 - \hat{p}) /
m}}\$\$

with the standard error worked out at the measured prevalence
\\\hat{p}\\. The one-sided p-values are the normal tail areas above and
below \\z\\ (t tail areas when the ICC is estimated; see Few sites). It
can be too optimistic for small samples or a prevalence near 0; with no
positives (or all positive) the standard error is 0, and any threshold
on the other side counts as rejected.

**Two-sided p-value.** For both methods, twice the smaller of the two
one-sided p-values, at most 1. For the exact test this is the "central"
two-sided exact test (Fay 2010), the one that matches the two-sided
Clopper-Pearson interval. It is cautious: the false alarm rate is at
most 1 - `conf_level`, and often well below it. It can differ from
[`binom.test()`](https://rdrr.io/r/stats/binom.test.html), which adds up
all counts that are less likely than the one observed.

**Interval.** The interval is the test turned around: it holds the
thresholds the test would not reject, so `reject` is `TRUE` exactly when
`threshold` is outside it. It is matched to `alternative`: \\\[L, 1\]\\
for `"greater"`, \\\[0, U\]\\ for `"less"`, and the two-sided interval
(\\\alpha/2\\ in each tail) for `"two.sided"`. For the exact test it is
the Clopper-Pearson interval,

\$\$L = B^{-1}(\alpha;\\ x\_{eff},\\ m - x\_{eff} + 1), \qquad U =
B^{-1}(1 - \alpha;\\ x\_{eff} + 1,\\ m - x\_{eff})\$\$

with \\x\_{eff} = \hat{p} m\\ and \\B^{-1}\\ the beta quantile function;
for the Wald test it is \\\hat{p} \pm z\_{1-\alpha} \sqrt{\hat{p}(1 -
\hat{p}) / m}\\ (\\\alpha/2\\ in place of \\\alpha\\ for two-sided).
With an estimated ICC, \\m\\ is the shrunk size in the Clopper-Pearson
interval and the t quantile replaces \\z\\ in the Wald interval. It is
built on the measured scale and then corrected for an imperfect test,
with ends kept within \[0, 1\]; if the whole corrected interval lies
outside \[0, 1\], the function warns that the data do not fit the
assumed sensitivity and specificity.

**When to use DRpower instead.** With the ICC supplied, this function
keeps the chance of a false alarm at or below the significance level.
With `icc = NULL` the degrees-of-freedom adjustment keeps it there too:
in simulation, with the true prevalence at a 5 percent threshold, an ICC
of 0.05 and 2 to 10 sites (equal or very different sizes), false alarms
at a 5 percent level were at most 4.6 percent. The price is low power
with few sites, as an ICC estimated from them says little. For pfhrp2/3
deletion surveys with several sites and an unknown ICC, where the
question is whether prevalence is above the threshold, consider DRpower
(<https://mrc-ide.github.io/DRpower/>). It estimates prevalence and ICC
together and has more power, especially with few sites, but more false
alarms when prevalence is close to the threshold.

## References

MMS-SD Study Design Workshop.
<https://mrc-ide.github.io/MMS-SD_workshop/>

Clopper, C. J. & Pearson, E. S. (1934) The use of confidence or fiducial
limits illustrated in the case of the binomial. *Biometrika* **26**(4),
404-413.
[doi:10.1093/biomet/26.4.404](https://doi.org/10.1093/biomet/26.4.404)

Fay, M. P. (2010) Two-sided exact tests and matching confidence
intervals for discrete data. *The R Journal* **2**(1), 53-58.
[doi:10.32614/RJ-2010-008](https://doi.org/10.32614/RJ-2010-008)

Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)

Korn, E. L. & Graubard, B. I. (1998) Confidence intervals for
proportions with small expected number of positive counts estimated from
survey data. *Survey Methodology* **24**(2), 193-201.
(Degrees-of-freedom adjustment.)

Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
of a screening test. *American Journal of Epidemiology* **107**(1),
71-76.
[doi:10.1093/oxfordjournals.aje.a112510](https://doi.org/10.1093/oxfordjournals.aje.a112510)

## See also

[`design_threshold`](https://mrc-ide.github.io/MMSpower/reference/design_threshold.md)

## Examples

``` r
# One site: is prevalence above 5 percent?
test_threshold(x = 8, n = 100, threshold = 0.05)
#> $p_value
#> [1] 0.1279605
#> 
#> $reject
#> [1] FALSE
#> 
#> $prevalence
#> [1] 0.08
#> 
#> $ci_lower
#> [1] 0.04042887
#> 
#> $ci_upper
#> [1] 1
#> 
#> $n_eff
#> [1] 100
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 
#> $df
#> [1] Inf
#> 

# Five sites, ICC estimated from the data
test_threshold(
  x = c(0, 6, 1, 11, 4),
  n = c(50, 60, 40, 80, 70),
  threshold = 0.05
)
#> $p_value
#> [1] 0.2631809
#> 
#> $reject
#> [1] FALSE
#> 
#> $prevalence
#> [1] 0.07333333
#> 
#> $ci_lower
#> [1] 0.0284146
#> 
#> $ci_upper
#> [1] 1
#> 
#> $n_eff
#> [1] 108.5637
#> 
#> $icc
#> [1] 0.02828913
#> 
#> $deff
#> [1] 2.763356
#> 
#> $df
#> [1] 4
#> 

# Two-sided, imperfect diagnostic test
test_threshold(x = 12, n = 100, threshold = 0.10,
               alternative = "two.sided",
               sensitivity = 0.9, specificity = 0.95)
#> $p_value
#> [1] 0.795509
#> 
#> $reject
#> [1] FALSE
#> 
#> $prevalence
#> [1] 0.08235294
#> 
#> $ci_lower
#> [1] 0.01596341
#> 
#> $ci_upper
#> [1] 0.1767479
#> 
#> $n_eff
#> [1] 100
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 
#> $df
#> [1] Inf
#> 

# One site, a clear result: the interval lies above the threshold
test_threshold(x = 15, n = 100, threshold = 0.05)
#> $p_value
#> [1] 0.0001358542
#> 
#> $reject
#> [1] TRUE
#> 
#> $prevalence
#> [1] 0.15
#> 
#> $ci_lower
#> [1] 0.09479401
#> 
#> $ci_upper
#> [1] 1
#> 
#> $n_eff
#> [1] 100
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 
#> $df
#> [1] Inf
#> 

# Wald test and interval instead of the exact test
test_threshold(x = 15, n = 100, threshold = 0.05, method = "wald")
#> $p_value
#> [1] 0.0025507
#> 
#> $reject
#> [1] TRUE
#> 
#> $prevalence
#> [1] 0.15
#> 
#> $ci_lower
#> [1] 0.09126698
#> 
#> $ci_upper
#> [1] 1
#> 
#> $n_eff
#> [1] 100
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 
#> $df
#> [1] Inf
#> 
```
