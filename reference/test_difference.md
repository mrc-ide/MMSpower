# Test whether prevalence differs between two groups

Analysis function: given the number of positives and the number tested
in two groups, tests whether their true prevalence differs. Returns the
estimated difference with the confidence interval that matches the test,
the p-value and whether the difference is significant. It is the
analysis twin of
[`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md).

The test compares the measured (apparent) prevalences, the share who
test positive. Because sensitivity and specificity are assumed the same
in both groups, the measured prevalences differ exactly when the true
ones do, in the same direction, so the p-value does not depend on
sensitivity or specificity. Only the reported difference and its
interval are corrected for an imperfect test (Rogan-Gladen correction).
Clustering is handled with a design effect for each group.

**Hypotheses** (\\p_1\\, \\p_2\\ are the true prevalences):

- `"two.sided"` (default):

  H0: \\p_1 = p_2\\ against H1: \\p_1 \ne p_2\\.

- `"greater"`:

  H0: \\p_1 \le p_2\\ against H1: \\p_1 \> p_2\\.

- `"less"`:

  H0: \\p_1 \ge p_2\\ against H1: \\p_1 \< p_2\\.

## Usage

``` r
test_difference(
  x1,
  n1,
  x2,
  n2,
  alternative = c("two.sided", "greater", "less"),
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  icc = NULL,
  method = c("score", "wald")
)
```

## Arguments

- x1:

  Vector of whole-number counts of people who tested positive in group
  1, one per site (a single number if group 1 is one site).

- n1:

  Vector of whole-number totals tested in group 1, one per site (same
  length as `x1`).

- x2:

  Vector of whole-number counts of people who tested positive in group
  2, one per site (a single number if group 2 is one site).

- n2:

  Vector of whole-number totals tested in group 2, one per site (same
  length as `x2`).

- alternative:

  Direction of the test: `"two.sided"` (default), `"greater"` (group 1
  higher) or `"less"` (group 1 lower).

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1. Assumed the same in both
  groups.

- specificity:

  Diagnostic specificity in (0, 1\]; default 1. Assumed the same in both
  groups. `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level for the interval on the difference, in (0, 1);
  default 0.95. The significance level of the test is 1 - `conf_level`.

- icc:

  Optional. Intra-cluster correlation (ICC), in \[0, 1\]. If `NULL`
  (default), it is estimated separately for each group from its own
  sites. Otherwise give one number used for both groups, or two,
  `c(group 1, group 2)`. A supplied value is also used for a group with
  only one site (see Details). Set to `0` for no clustering.

- method:

  `"score"` (default): pooled z-test with the score interval. `"wald"`:
  unpooled (Wald) z-test with the Wald interval. Each interval is its
  own test turned around, so the interval leaves out 0 exactly when
  `reject` is `TRUE`, whichever method is used.

## Value

A list with:

- difference:

  Estimated true difference, group 1 minus group 2: the measured
  difference divided by \\Se + Sp - 1\\, between -1 and 1.

- ci_lower, ci_upper:

  Confidence interval for the true difference; `ci_upper` is 1 for
  `"greater"` and `ci_lower` is -1 for `"less"`.

- p_value:

  p-value for the chosen `alternative`.

- reject:

  `TRUE` if `p_value` is below 1 - `conf_level`.

- prevalence1, prevalence2:

  Estimated true prevalence in each group (Rogan-Gladen corrected,
  between 0 and 1).

- n_eff1, n_eff2:

  Effective sample size of each group, `sum(n) / deff`.

- icc1, icc2:

  ICC used for each group: as supplied, or estimated from that group's
  data.

- deff1, deff2:

  Design effect used for each group (1 with no clustering).

- df:

  Degrees of freedom of the test (see Few sites); `Inf` when neither ICC
  is estimated.

## Details

**Design effect** (returned as `deff1` and `deff2`). For each group it
is \\1 + (n_w - 1) ICC\\, where \\n_w = \sum n^2 / \sum n\\ is the mean
number of people per site, weighted by site size. With `icc` supplied
(for example from earlier surveys), that value is used. Supply `icc` for
a group with one site only when the site stands for a wider population,
such as a district; if the site itself is the population of interest,
leave `icc = NULL`. With `icc = NULL`, the ICC of each group is
estimated from its own data as in
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md);
this needs at least two sites, and a group with one site gets a design
effect of 1.

**Few sites.** An ICC estimated from \\k\\ sites rests on only \\k - 1\\
degrees of freedom, so the design effect is itself uncertain. When
either group's ICC is estimated, the test statistic is compared with a t
distribution in place of the normal, and the interval uses the matching
t quantile. Its degrees of freedom (returned as `df`) combine the two
groups by Satterthwaite's (1946) formula,

\$\$df = \frac{(v_1 + v_2)^2}{v_1^2 / df_1 + v_2^2 / df_2}, \qquad v_g =
\hat{p}\_g (1 - \hat{p}\_g) / m_g\$\$

with \\df_g = k_g - 1\\ for a group whose ICC is estimated, and a group
whose ICC is supplied (or that has one site) left out of the
denominator. With both ICCs supplied, `df` is `Inf` and the test is the
normal one. In simulations with a true ICC of 0.1, 50 people per site
and prevalences of 0.2 to 0.5, a 95 percent interval contained the true
difference 96 to 97 percent of the time with 3 sites per group and 95 to
96 percent with 5. With 2 sites per group it is very wide (99 percent),
an honest picture of how little 2 sites say.

**Effective sample size.** For group \\g\\, \\m_g\\ (returned as
`n_eff1` and `n_eff2`) is the number tested divided by the design
effect.

**Score test** (`method = "score"`). With \\\hat{p}\_g\\ the measured
prevalence in group \\g\\,

\$\$z = \frac{\hat{p}\_1 - \hat{p}\_2}
{\sqrt{\bar{p}(1-\bar{p})\\(1/m_1 + 1/m_2)}}, \qquad \bar{p} =
\frac{\hat{p}\_1 m_1 + \hat{p}\_2 m_2}{m_1 + m_2}.\$\$

With no clustering this is the same test as R's
`prop.test(correct = FALSE)`, and the one
[`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md)
sizes for.

**Wald test** (`method = "wald"`). The same, but with the standard error
worked out from each group's own prevalence,

\$\$z = \frac{\hat{p}\_1 - \hat{p}\_2}
{\sqrt{\hat{p}\_1(1-\hat{p}\_1)/m_1 +
\hat{p}\_2(1-\hat{p}\_2)/m_2}}.\$\$

With small groups it can be too optimistic: with 20 people per group and
true prevalences of 0.05 and 0.10, a 95 percent Wald interval contains
the true difference only about 90 percent of the time (score 97). When
each group's measured prevalence is 0 or 1 the standard error is 0.

For both, the one-sided p-values are the normal tail areas above and
below \\z\\ (t tail areas with `df` degrees of freedom when an ICC is
estimated), and the two-sided p-value is twice the smaller of them.

**Interval.** The interval is the test turned around: it holds the
differences \\\delta\\ the test would not reject, so it leaves out 0
exactly when `reject` is `TRUE`. It matches `alternative`: two-sided for
`"two.sided"`, \\\[L, 1\]\\ for `"greater"` and \\\[-1, U\]\\ for
`"less"`.

For `"score"`, the two prevalences are estimated again by maximum
likelihood under the condition \\p_1 - p_2 = \delta\\ (the closed form
of Farrington and Manning, 1990), giving \\\tilde{p}\_1, \tilde{p}\_2\\,
and

\$\$z(\delta) = \frac{\hat{p}\_1 - \hat{p}\_2 - \delta}
{\sqrt{\tilde{p}\_1(1-\tilde{p}\_1)/m_1 +
\tilde{p}\_2(1-\tilde{p}\_2)/m_2}}\$\$

is set equal to the z value for `conf_level` and solved for \\\delta\\
(Mee, 1984). At \\\delta = 0\\, \\\tilde{p}\_1 = \tilde{p}\_2 =
\bar{p}\\ and this is exactly the score test above. Miettinen and
Nurminen (1985) multiply the variance by \\(m_1 + m_2)/(m_1 + m_2 -
1)\\; that factor is left out here so that the interval and the test
always agree.

For `"wald"`, the interval is \\\hat{p}\_1 - \hat{p}\_2 \pm z
\sqrt{\hat{p}\_1(1 - \hat{p}\_1)/m_1 + \hat{p}\_2(1 -
\hat{p}\_2)/m_2}\\.

The interval is built on the measured scale, then divided by \\Se + Sp -
1\\ to give an interval for the true difference (the rest of the
Rogan-Gladen correction cancels in a difference). Ends outside -1 to 1
are set to -1 or 1.

## References

MMS-SD Study Design Workshop.
<https://mrc-ide.github.io/MMS-SD_workshop/>

Fleiss, J. L., Levin, B. & Paik, M. C. (2003) *Statistical Methods for
Rates and Proportions*, 3rd ed. Wiley. (Two-group z-test.)

Mee, R. W. (1984) Confidence bounds for the difference between two
probabilities. *Biometrics* **40**(4), 1175-1176.

Miettinen, O. & Nurminen, M. (1985) Comparative analysis of two rates.
*Statistics in Medicine* **4**(2), 213-226.

Farrington, C. P. & Manning, G. (1990) Test statistics and sample size
formulae for comparative binomial trials with null hypothesis of
non-zero risk difference or non-unity relative risk. *Statistics in
Medicine* **9**(12), 1447-1454. (Prevalences estimated under a given
difference.)

Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)

Satterthwaite, F. E. (1946) An approximate distribution of estimates of
variance components. *Biometrics Bulletin* **2**(6), 110-114.
[doi:10.2307/3002019](https://doi.org/10.2307/3002019)

Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
of a screening test. *American Journal of Epidemiology* **107**(1),
71-76.
[doi:10.1093/oxfordjournals.aje.a112510](https://doi.org/10.1093/oxfordjournals.aje.a112510)

## See also

[`design_difference`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md)

## Examples

``` r
# One count per group: 8 of 100 versus 20 of 120
test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)
#> $difference
#> [1] -0.08666667
#> 
#> $ci_lower
#> [1] -0.1746785
#> 
#> $ci_upper
#> [1] 0.001913274
#> 
#> $p_value
#> [1] 0.05478913
#> 
#> $reject
#> [1] FALSE
#> 
#> $prevalence1
#> [1] 0.08
#> 
#> $prevalence2
#> [1] 0.1666667
#> 
#> $n_eff1
#> [1] 100
#> 
#> $n_eff2
#> [1] 120
#> 
#> $icc1
#> [1] 0
#> 
#> $icc2
#> [1] 0
#> 
#> $deff1
#> [1] 1
#> 
#> $deff2
#> [1] 1
#> 
#> $df
#> [1] Inf
#> 

# Three sites per group, one-sided, imperfect diagnostic test. The sites
# differ, so the estimated design effects are above 1 (about 2.7 and 3.4)
# and p is 0.12; treating everyone as independent would give 0.007.
test_difference(
  x1 = c(1, 9, 3), n1 = c(50, 60, 40),
  x2 = c(4, 18, 9), n2 = c(55, 70, 45),
  alternative = "less",
  sensitivity = 0.9, specificity = 0.97
)
#> $difference
#> [1] -0.1099842
#> 
#> $ci_lower
#> [1] -1
#> 
#> $ci_upper
#> [1] 0.06249794
#> 
#> $p_value
#> [1] 0.1151506
#> 
#> $reject
#> [1] FALSE
#> 
#> $prevalence1
#> [1] 0.0651341
#> 
#> $prevalence2
#> [1] 0.1751183
#> 
#> $n_eff1
#> [1] 55.38594
#> 
#> $n_eff2
#> [1] 49.50051
#> 
#> $icc1
#> [1] 0.03393911
#> 
#> $icc2
#> [1] 0.04231414
#> 
#> $deff1
#> [1] 2.708268
#> 
#> $deff2
#> [1] 3.434308
#> 
#> $df
#> [1] 3.549145
#> 

# Wald test and interval instead of the default score versions
test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, method = "wald")
#> $difference
#> [1] -0.09583333
#> 
#> $ci_lower
#> [1] -0.1774752
#> 
#> $ci_upper
#> [1] -0.01419142
#> 
#> $p_value
#> [1] 0.02141113
#> 
#> $reject
#> [1] TRUE
#> 
#> $prevalence1
#> [1] 0.0375
#> 
#> $prevalence2
#> [1] 0.1333333
#> 
#> $n_eff1
#> [1] 80
#> 
#> $n_eff2
#> [1] 90
#> 
#> $icc1
#> [1] 0
#> 
#> $icc2
#> [1] 0
#> 
#> $deff1
#> [1] 1
#> 
#> $deff2
#> [1] 1
#> 
#> $df
#> [1] Inf
#> 

# A different ICC for each group
test_difference(
  x1 = c(2, 5, 3), n1 = c(50, 60, 40),
  x2 = c(8, 11, 9), n2 = c(55, 70, 45),
  icc = c(0.02, 0.08)
)
#> $difference
#> [1] -0.09803922
#> 
#> $ci_lower
#> [1] -0.2720197
#> 
#> $ci_upper
#> [1] 0.02312369
#> 
#> $p_value
#> [1] 0.1202936
#> 
#> $reject
#> [1] FALSE
#> 
#> $prevalence1
#> [1] 0.06666667
#> 
#> $prevalence2
#> [1] 0.1647059
#> 
#> $n_eff1
#> [1] 74.75083
#> 
#> $n_eff2
#> [1] 30.34439
#> 
#> $icc1
#> [1] 0.02
#> 
#> $icc2
#> [1] 0.08
#> 
#> $deff1
#> [1] 2.006667
#> 
#> $deff2
#> [1] 5.602353
#> 
#> $df
#> [1] Inf
#> 
```
