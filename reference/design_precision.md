# Calculate sample size for a target margin of error when estimating prevalence

Given a target precision (margin of error), returns the minimum total
sample size required. Handles imperfect diagnostic tests via the
Rogan-Gladen variance adjustment, and clustered sampling via the design
effect (Kish formula).

**Three design modes** – controlled by `n_sites` and `n_per_site`:

- Unclustered (`n_sites = NULL`, `n_per_site = NULL`):

  Treats all observations as independent. Returns total `n` only.

- Fixed site size (`n_per_site` supplied):

  Solves for the number of sites.

- Fixed number of sites (`n_sites` supplied):

  Solves for the number of people per site.

## Usage

``` r
design_precision(
  prevalence,
  moe,
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
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
  +/-5 percentage points.

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1 (perfect test).

- specificity:

  Diagnostic specificity in (0, 1\]; default 1.
  `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level, in (0, 1); default 0.95.

- n_sites:

  Optional whole number. Fixes the number of sites; the function solves
  for people per site. Cannot be used with `n_per_site`.

- n_per_site:

  Optional whole number. Fixes the number of people per site; the
  function solves for the number of sites. Cannot be used with
  `n_sites`.

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

  Unclustered sample size \\n_0\\ that gives the same precision, rounded
  up.

## Details

**Rogan-Gladen variance adjustment**: the Rogan-Gladen transform is
linear in the apparent prevalence, so scaling its variance by the
(constant) RG slope is exact – not a delta-method approximation. The
unclustered sample size is

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

\$\$MOE\_{min} = z\sqrt{\frac{p\_{app}(1-p\_{app}) \cdot ICC}
{n\_{sites} \cdot (Se + Sp - 1)^2}}\$\$

and the function stops and reports this minimum.

**Rounding**: the number of sites (or people per site) is rounded up,
and the returned `n` is the number actually sampled,
`n_sites * n_per_site`, so it can be slightly above the continuous
solution.

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

# Imperfect test (sensitivity 90%, specificity 95%)
design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)
#> $n
#> [1] 451
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
#> [1] 451
#> 

# Fixed site size: how many sites do I need?
design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
#> $n
#> [1] 470
#> 
#> $n_sites
#> [1] 47
#> 
#> $n_per_site
#> [1] 10
#> 
#> $deff
#> [1] 1.45
#> 
#> $n_eff
#> [1] 323
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
#> [1] 323
#> 
```
