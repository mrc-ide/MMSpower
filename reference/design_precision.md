# Calculate sample size for a target margin of error when estimating prevalence

Given a target precision (margin of error), returns the minimum total
sample size required; or, given a full design, returns the margin of
error it gives. Handles imperfect diagnostic tests via the Rogan-Gladen
variance adjustment, and clustered sampling via the design effect (Kish
formula).

**What is solved for** – leave out the one quantity you want:

- `moe` only:

  Returns the number of people, with no sites.

- `moe` and `n_per_site`:

  Returns the number of sites.

- `moe` and `n_sites`:

  Returns the number of people per site.

- A full design, without `moe`:

  Returns the margin of error: give `n` (no sites), `n_sites` and
  `n_per_site`, or a vector `n_per_site` with one size per site (for
  sites of different sizes).

## Usage

``` r
design_precision(
  prevalence,
  moe = NULL,
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  ci_method = c("clopper-pearson", "wald", "agresti-coull"),
  n = NULL,
  n_sites = NULL,
  n_per_site = NULL,
  icc = 0
)
```

## Arguments

- prevalence:

  Expected true prevalence, in (0, 1).

- moe:

  Target margin of error (half-width of the confidence interval) on the
  true-prevalence scale, in (0, 0.5). A decimal fraction, e.g. 0.05 for
  +/-5 percentage points. Leave out with a full design.

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1 (perfect test).

- specificity:

  Diagnostic specificity in (0, 1\]; default 1.
  `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level, in (0, 1); default 0.95.

- ci_method:

  The interval the study will be analysed with, as in
  [`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md):
  `"clopper-pearson"` (default, the same default as
  [`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)),
  `"wald"` or `"agresti-coull"`. See Details.

- n:

  Optional whole number: the number of people in a full design with no
  sites. Cannot be used with `n_sites` or `n_per_site`.

- n_sites:

  Optional whole number. Fixes the number of sites; alone, the function
  solves for people per site.

- n_per_site:

  Optional. People per site: one whole number (alone, the function
  solves for the number of sites), or a vector with one size per site, a
  full design with sites of different sizes (then leave out `n_sites`).

- icc:

  Intra-cluster correlation, in \[0, 1\]; default 0 (no clustering). An
  `icc` above 0 needs `n_sites` or `n_per_site`.

## Value

A list with:

- n:

  Total sample size. With sites, `n_sites * n_per_site`.

- n_sites:

  Number of sites (`NULL` if unclustered).

- n_per_site:

  People per site (`NULL` if unclustered).

- deff:

  Design effect of the returned design, `1 + (n_per_site - 1) * icc`; 1
  if unclustered.

- n_eff:

  Unclustered sample size \\n_0\\ for the target precision, rounded up;
  for a full design, `n / deff`.

- moe:

  Margin of error of the returned design: the half-width of its expected
  interval.

## Details

**Margin of error and the interval method**: the margin of error is the
half-width of the interval
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)
would report if the measured prevalence came out exactly as expected, so
a design sized here gives back the target when analysed with the same
`ci_method`. For `"wald"` it is the usual closed form below. For
`"clopper-pearson"` and `"agresti-coull"` there is no formula: the
unclustered size \\n_0\\ is found by a numerical search, and the
clustering steps below then work as for Wald. Clopper-Pearson is the
widest of the three, so it needs the most people (341 against Wald's 323
for 30 percent +/- 5); the gap is largest at low prevalence.

**Rogan-Gladen variance adjustment**: the Rogan-Gladen transform is
linear in the apparent prevalence, so scaling its variance by the
(constant) RG slope is exact – not a delta-method approximation. The
unclustered Wald sample size is

\$\$n_0 = \frac{z^2 \\ p\_{app}(1-p\_{app})}{MOE^2 \cdot (Se + Sp -
1)^2}\$\$

where \\p\_{app} = p \cdot Se + (1-p)(1-Sp)\\ is the apparent prevalence
implied by the expected true prevalence and the test characteristics.
When \\Se = Sp = 1\\ this reduces to the standard MOE formula.

**Fixed site size**: the design effect \\D\_{eff} = 1 + (n\_{per
site} - 1) \cdot ICC\\ is known directly, so \\n = n_0 \cdot D\_{eff}\\
and the number of sites is `ceiling(n / n_per_site)`.

**Fixed number of sites**: Deff depends on the number of people per
site, \\n / n\_{sites}\\, which depends on n. Substituting and solving
gives the closed form:

\$\$n = \frac{n_0 \cdot n\_{sites} \cdot (1-ICC)}{n\_{sites} - n_0 \cdot
ICC}\$\$

A solution exists only when \\n\_{sites} \> n_0 \cdot ICC\\. Otherwise
the target MOE is unachievable: each site is worth at most \\1/ICC\\
independent people however many are sampled there, so the MOE floors at
the margin of error of \\n\_{sites} / ICC\\ independent people (for
Wald, \\z\sqrt{p\_{app}(1-p\_{app}) ICC / n\_{sites}} / (Se + Sp -
1)\\), and the function stops and reports this minimum.

**Rounding**: the number of sites (or people per site) is rounded up,
and the returned `n` is the number actually sampled,
`n_sites * n_per_site`, so it can be slightly above the continuous
solution. The returned `moe` is that of the rounded design, so it can be
slightly below the target.

**Full design**: the margin of error is that of \\n / D\_{eff}\\
independent people, for Wald \\z \sqrt{p\_{app}(1-p\_{app}) D\_{eff} /
n} / (Se + Sp - 1)\\. For sites of different sizes, \\D\_{eff} = 1 +
(n_w - 1) ICC\\ with \\n_w = \sum s^2 / \sum s\\ the size-weighted mean
site size, as in
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md).

## References

Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
screening test. American Journal of Epidemiology 107(1):71-76.

MMS-SD Study Design Workshop, Modules 1 (sampling), 2 (sample size from
margin of error) and 5 (ICC / design effect).
<https://mrc-ide.github.io/MMS-SD_workshop/>

## Examples

``` r
# Unclustered, perfect test
design_precision(prevalence = 0.3, moe = 0.05)
#> $n
#> [1] 341
#> 
#> $n_sites
#> NULL
#> 
#> $n_per_site
#> NULL
#> 
#> $deff
#> [1] 1
#> 
#> $n_eff
#> [1] 341
#> 
#> $moe
#> [1] 0.04993769
#> 

# Sized for a Wald interval instead of the default Clopper-Pearson
design_precision(prevalence = 0.3, moe = 0.05, ci_method = "wald")
#> $n
#> [1] 323
#> 
#> $n_sites
#> NULL
#> 
#> $n_per_site
#> NULL
#> 
#> $deff
#> [1] 1
#> 
#> $n_eff
#> [1] 323
#> 
#> $moe
#> [1] 0.04997542
#> 

# Imperfect test (sensitivity 90%, specificity 95%)
design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)
#> $n
#> [1] 472
#> 
#> $n_sites
#> NULL
#> 
#> $n_per_site
#> NULL
#> 
#> $deff
#> [1] 1
#> 
#> $n_eff
#> [1] 472
#> 
#> $moe
#> [1] 0.04999082
#> 

# Fixed site size: how many sites do I need?
design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
#> $n
#> [1] 500
#> 
#> $n_sites
#> [1] 50
#> 
#> $n_per_site
#> [1] 10
#> 
#> $deff
#> [1] 1.45
#> 
#> $n_eff
#> [1] 341
#> 
#> $moe
#> [1] 0.04965355
#> 

# Fixed number of sites: how many people per site?
design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)
#> $n
#> [1] 500
#> 
#> $n_sites
#> [1] 50
#> 
#> $n_per_site
#> [1] 10
#> 
#> $deff
#> [1] 1.45
#> 
#> $n_eff
#> [1] 341
#> 
#> $moe
#> [1] 0.04965355
#> 

# Full design: what margin of error do 6 sites of 40 give?
design_precision(0.2, n_sites = 6, n_per_site = 40, icc = 0.05)
#> $n
#> [1] 240
#> 
#> $n_sites
#> [1] 6
#> 
#> $n_per_site
#> [1] 40
#> 
#> $deff
#> [1] 2.95
#> 
#> $n_eff
#> [1] 81.35593
#> 
#> $moe
#> [1] 0.0919475
#> 

# Sites of different sizes
design_precision(0.2, n_per_site = c(250, 250, 60, 60, 40, 40), icc = 0.05)
#> $n
#> [1] 700
#> 
#> $n_sites
#> [1] 6
#> 
#> $n_per_site
#> [1] 250 250  60  60  40  40
#> 
#> $deff
#> [1] 10.62143
#> 
#> $n_eff
#> [1] 65.90451
#> 
#> $moe
#> [1] 0.1025958
#> 
```
