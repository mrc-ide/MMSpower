# Estimate prevalence and its precision

Analysis function: given observed counts, returns a prevalence point
estimate and confidence interval. Three CI methods are available:

- `"wald"` (default):

  Symmetric unless an end is moved to 0 or 1. \\\hat{p} \pm z \cdot
  SE\\. Can be too narrow for small samples or prevalence near 0 or 1.

- `"clopper-pearson"`:

  Asymmetric. Exact binomial interval, via the beta distribution. For a
  simple random sample, its coverage is at least the stated confidence
  level, so it is usually wider than the Wald interval.

- `"agresti-coull"`:

  Asymmetric. A Wald interval calculated after adding \\z^2/2\\ positive
  and \\z^2/2\\ negative cases to the data (about 1.92 of each at 95
  percent). This moves the interval towards 0.5; the reported prevalence
  is unchanged.

All methods account for imperfect diagnostic tests (Rogan-Gladen
correction) and clustered sampling (Kish design effect).

If `icc` is not supplied, it is estimated from the data (the observed
variance between sites compared with the variance expected if every
person were independent). Setting `icc = 0` treats every person as
independent; if people at the same site tend to be similar, this gives
an interval that is too narrow.

## Usage

``` r
estimate_prevalence(
  x,
  n,
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  icc = NULL,
  ci_method = c("wald", "clopper-pearson", "agresti-coull")
)
```

## Arguments

- x:

  Vector of whole-number counts of people who tested positive, one per
  site (a single number for one site).

- n:

  Vector of whole-number totals tested, one per site (same length as
  `x`).

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1 (the test never misses a
  case).

- specificity:

  Diagnostic specificity in (0, 1\]; default 1 (the test never gives a
  false positive). `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level, in (0, 1); default 0.95.

- icc:

  Optional. Intra-cluster correlation, in \[0, 1\]. If `NULL` (default),
  ICC is estimated from the data. Set to `0` for no clustering. A
  supplied value is also used with a single site.

- ci_method:

  CI method: `"wald"` (default), `"clopper-pearson"`, or
  `"agresti-coull"`. See Description.

## Value

A list with:

- prevalence:

  Point estimate of true prevalence (Rogan-Gladen corrected).

- ci_lower, ci_upper:

  Confidence limits on the true-prevalence scale.

- moe:

  Half-width of the interval, `(ci_upper - ci_lower) / 2`.

- n_eff:

  Effective sample size, `sum(n) / deff`.

- icc:

  ICC used: as supplied, or estimated from the data.

- deff:

  Design effect applied (1 when there is no clustering).

## Details

The interval is built in five steps: apparent prevalence, a design
effect for clustering, an effective sample size, a confidence interval
by the chosen method, and finally the Rogan-Gladen correction for an
imperfect test.

**1. Apparent prevalence.** The pooled proportion of test positives
across all sites,

\$\$\hat{p} = \frac{\sum_i x_i}{\sum_i n_i}\$\$

It is corrected for test imperfection only at step 5, so the correction
applies identically to the point estimate and to both interval
endpoints.

**2. Design effect (Kish).** Observations within a site are positively
correlated, so a clustered sample carries less information than the same
number of independent people. The Kish (1965) design effect scales the
variance by

\$\$D\_{eff} = 1 + (n_w - 1)\\\rho\$\$

where \\n_w = \sum n^2 / \sum n\\ is the mean number of people per site,
weighted by site size (the plain mean when all sites are the same size),
and \\\rho\\ the intra-cluster correlation (ICC).

*ICC supplied*: \\\rho\\ is used directly, also for a single site.

*ICC estimated* (`icc = NULL`, the default): with two or more sites and
an average of more than one person per site, the design effect is first
estimated as the ratio of the observed variance of the site proportions
to the variance expected under simple random sampling,

\$\$\widehat{D} =
\frac{\mathrm{Var}(\hat{p}\_i)}{\\\overline{\hat{p}(1 - \hat{p}) /
n_i}\\}, \qquad \hat{p}\_i = \frac{x_i}{n_i}\$\$

and the ICC recovered by inverting the Kish formula, \\\hat{\rho} =
(\widehat{D} - 1) / (\bar{n} - 1)\\, kept between 0 and

1.  \\D\_{eff}\\ is then computed from \\\hat{\rho}\\ as above. With a
    single site, one person per site, or nobody (or everybody) positive,
    the ICC is 0. If everyone tested comes from one site, supply `icc`.

An ICC estimated from only a few sites is uncertain, and the interval
then tends to be too narrow. In simulations with a true ICC of 0.1, 50
people per site and a prevalence of 0.3, a 95 percent Wald interval
contained the true prevalence about 76 percent of the time with 2 sites,
87 percent with 5 and 90 percent with 10 (Clopper-Pearson: 80, 89 and 92
percent). With the true ICC supplied it was 92 to 95 percent. Supplying
`icc` from other studies avoids this.

**3. Effective sample size.**

\$\$n\_{eff} = \frac{n\_{total}}{D\_{eff}}\$\$

the number of independent observations carrying the same information as
the clustered sample.

**4. Confidence interval on apparent prevalence.** Let \\z =
\Phi^{-1}(1 - \alpha/2)\\, \\\alpha = 1 - \\ `conf_level`, and
\\x\_{eff} = \hat{p}\\n\_{eff}\\, the number of positives scaled to the
effective sample size. Any end below 0 is set to 0 and any end above 1
is set to 1.

*`"wald"`* – the normal-approximation interval,

\$\$\hat{p} \\\pm\\ z \sqrt{\frac{\hat{p}(1 - \hat{p})}{n\_{eff}}}\$\$

For small \\n\\ or prevalence near 0 or 1 it can be too narrow, and it
has zero width when \\\hat{p} = 0\\ or \\1\\. For example, with 30
people and a true prevalence of 5 percent, a 95 percent Wald interval
contains the true prevalence only about 78 percent of the time.

*`"clopper-pearson"`* – the exact binomial interval,

\$\$L = B^{-1}\\\left(\tfrac{\alpha}{2};\\ x\_{eff},\\ n\_{eff} -
x\_{eff} + 1\right)\$\$ \$\$U = B^{-1}\\\left(1 - \tfrac{\alpha}{2};\\
x\_{eff} + 1,\\ n\_{eff} - x\_{eff}\right)\$\$

where \\B^{-1}(q; a, b)\\ is the \\q\\-th quantile of a beta
distribution, and \\L = 0\\ when \\\hat{p} = 0\\ and \\U = 1\\ when
\\\hat{p} = 1\\. \\x\_{eff}\\ need not be a whole number, so the design
effect can be applied; coverage of at least the stated level then holds
only approximately.

*`"agresti-coull"`* – add \\z^2/2\\ positive and \\z^2/2\\ negative
cases, then take a Wald interval on the adjusted proportion,

\$\$\tilde{n} = n\_{eff} + z^2, \qquad \tilde{p} = \frac{x\_{eff} +
z^2/2}{\tilde{n}}\$\$ \$\$\tilde{p} \\\pm\\ z \sqrt{\frac{\tilde{p}(1 -
\tilde{p})}{\tilde{n}}}\$\$

**5. Rogan-Gladen correction.** An imperfect test inflates apparent
prevalence through false positives and deflates it through false
negatives. Inverting \\p\_{app} = p\\Se + (1 - p)(1 - Sp)\\ gives true
prevalence,

\$\$p = \frac{p\_{app} - (1 - Sp)}{Se + Sp - 1}\$\$

This same formula is applied to \\\hat{p}\\ and to both CI endpoints,
and the results are kept within \[0, 1\]. If the whole corrected
interval lies below 0 (fewer positives than false positives alone would
produce) or above 1, the function warns that the data do not fit the
assumed sensitivity and specificity.

**Margin of error.** `moe` is the half-width of the corrected interval,
`(ci_upper - ci_lower) / 2`. Clopper-Pearson and Agresti-Coull intervals
are asymmetric, and so is a Wald interval with an end moved to 0 or 1;
report `ci_lower` and `ci_upper` then.

## References

MMS-SD Study Design Workshop.
<https://mrc-ide.github.io/MMS-SD_workshop/>

Kish, L. (1965) *Survey Sampling*. Wiley.

Clopper, C. J. & Pearson, E. S. (1934) The use of confidence or fiducial
limits illustrated in the case of the binomial. *Biometrika* **26**(4),
404-413.
[doi:10.1093/biomet/26.4.404](https://doi.org/10.1093/biomet/26.4.404)

Agresti, A. & Coull, B. A. (1998) Approximate is better than "exact" for
interval estimation of binomial proportions. *The American Statistician*
**52**(2), 119-126.
[doi:10.1080/00031305.1998.10480550](https://doi.org/10.1080/00031305.1998.10480550)

Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
of a screening test. *American Journal of Epidemiology* **107**(1),
71-76.
[doi:10.1093/oxfordjournals.aje.a112510](https://doi.org/10.1093/oxfordjournals.aje.a112510)

## Examples

``` r
# Single site / simple random sample
estimate_prevalence(x = 8, n = 50)
#> $prevalence
#> [1] 0.16
#> 
#> $ci_lower
#> [1] 0.05838385
#> 
#> $ci_upper
#> [1] 0.2616161
#> 
#> $moe
#> [1] 0.1016161
#> 
#> $n_eff
#> [1] 50
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 

# Multi-site clustered data -- ICC estimated from the data
estimate_prevalence(
  x = c(0, 4, 0, 22, 25, 16, 12, 8),
  n = c(60, 80, 70, 100, 40, 60, 50, 90)
)
#> $prevalence
#> [1] 0.1581818
#> 
#> $ci_lower
#> [1] 0.01536317
#> 
#> $ci_upper
#> [1] 0.3010005
#> 
#> $moe
#> [1] 0.1428186
#> 
#> $n_eff
#> [1] 25.07848
#> 
#> $icc
#> [1] 0.2867281
#> 
#> $deff
#> [1] 21.93115
#> 

# Clopper-Pearson interval
estimate_prevalence(x = 3, n = 30, ci_method = "clopper-pearson")
#> $prevalence
#> [1] 0.1
#> 
#> $ci_lower
#> [1] 0.02111714
#> 
#> $ci_upper
#> [1] 0.2652885
#> 
#> $moe
#> [1] 0.1220857
#> 
#> $n_eff
#> [1] 30
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 

# Imperfect diagnostic test
estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95)
#> $prevalence
#> [1] 0.2941176
#> 
#> $ci_lower
#> [1] 0.1884508
#> 
#> $ci_upper
#> [1] 0.3997845
#> 
#> $moe
#> [1] 0.1056669
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

# Known ICC
estimate_prevalence(x = c(5, 8, 3), n = c(40, 40, 40), icc = 0.05)
#> $prevalence
#> [1] 0.1333333
#> 
#> $ci_lower
#> [1] 0.02887
#> 
#> $ci_upper
#> [1] 0.2377967
#> 
#> $moe
#> [1] 0.1044633
#> 
#> $n_eff
#> [1] 40.67797
#> 
#> $icc
#> [1] 0.05
#> 
#> $deff
#> [1] 2.95
#> 
```
