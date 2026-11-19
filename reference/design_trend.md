# Calculate sample size per timepoint to detect a linear trend in prevalence

Given the prevalence at the first timepoint and either the prevalence at
the last timepoint or a slope, returns the number of people needed **per
timepoint** to detect a linear trend with the target power, using the
test in
[`test_trend()`](https://mrc-ide.github.io/MMSpower/reference/test_trend.md).
Supply `n_per_timepoint` instead to run it in reverse: the function then
reports the power that number of people per timepoint gives. Every
timepoint gets the same design.

**Three design modes** – controlled by `n_sites` and `n_per_site`, at
each timepoint:

- Unclustered (`n_sites = NULL`, `n_per_site = NULL`):

  Returns the number of people per timepoint.

- Fixed site size (`n_per_site` supplied):

  Solves for the number of sites per timepoint.

- Fixed number of sites (`n_sites` supplied):

  Solves for the number of people per site.

## Usage

``` r
design_trend(
  prevalence_start,
  prevalence_end = NULL,
  slope = NULL,
  n_timepoints = NULL,
  times = NULL,
  n_per_timepoint = NULL,
  power = 0.8,
  alternative = c("two.sided", "greater", "less"),
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  n_sites = NULL,
  n_per_site = NULL,
  icc = 0
)
```

## Arguments

- prevalence_start:

  Expected true prevalence at the first timepoint, in (0, 1).

- prevalence_end:

  Expected true prevalence at the last timepoint, in (0, 1). Supply
  exactly one of `prevalence_end` and `slope`.

- slope:

  Change in true prevalence per unit of time. Supply exactly one of
  `prevalence_end` and `slope`.

- n_timepoints:

  Whole number, 2 or more: equally spaced timepoints 0, 1, ...,
  `n_timepoints` - 1. Supply exactly one of `n_timepoints` and `times`.

- times:

  Numeric vector of at least two timepoints, for unequal spacing or real
  units (for example calendar years).

- n_per_timepoint:

  Optional whole number. If supplied, the function runs in reverse and
  reports the power this many people per timepoint give; `power` is then
  not used.

- power:

  Target power, in \[0.5, 1); default 0.80.

- alternative:

  Direction of the test, as used later in
  [`test_trend()`](https://mrc-ide.github.io/MMSpower/reference/test_trend.md):
  `"two.sided"` (default), `"greater"` (an increase) or `"less"` (a
  decrease).

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1.

- specificity:

  Diagnostic specificity in (0, 1\]; default 1.
  `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level, in (0, 1); default 0.95. The significance level of
  the test is 1 - `conf_level`.

- n_sites:

  Optional whole number. Fixes the number of sites per timepoint; the
  function solves for people per site. Cannot be used with `n_per_site`.

- n_per_site:

  Optional whole number. Fixes the number of people per site; the
  function solves for the number of sites per timepoint. Cannot be used
  with `n_sites`.

- icc:

  Intra-cluster correlation, in \[0, 1\]; default 0 (no clustering). An
  `icc` above 0 needs `n_sites` or `n_per_site`.

## Value

A list with:

- n_per_timepoint:

  People per timepoint (or as supplied, in reverse). With sites,
  `n_sites * n_per_site`.

- n:

  Total over all timepoints.

- n_sites:

  Sites per timepoint (`NULL` if unclustered).

- n_per_site:

  People per site (`NULL` if unclustered).

- deff:

  Design effect at each timepoint; 1 if unclustered.

- n_eff:

  Unclustered people per timepoint \\n_0\\, rounded up (in reverse,
  `n_per_timepoint / deff`).

- power:

  The target (forward), or the power achieved (reverse).

## Details

**Sample size per timepoint.** The test fits a straight line through the
measured prevalences at times \\t_1, \dots, t_T\\. With \\n\\ people per
timepoint, its slope has variance about \\\bar\sigma^2\_{app} / (n
S\_{xx})\\, where \\S\_{xx} = \sum_i (t_i - \bar t)^2\\ and
\\\bar\sigma^2\_{app} = \bar p\_{app}(1 - \bar p\_{app})\\ uses the
average of the first and last measured prevalences. Detecting a measured
slope \\\beta\_{app}\\ with the target power needs

\$\$n_0 = \frac{(z\_\alpha + z\_\beta)^2 \\ \bar\sigma^2\_{app}}
{S\_{xx} \\ \beta\_{app}^2},\$\$

with \\z\_\beta\\ the normal quantile at `power` and \\z\_\alpha\\ at
\\1 - \alpha/2\\ (two-sided) or \\1 - \alpha\\ (one-sided). The measured
slope is the true slope times \\Se + Sp - 1\\, so an imperfect test
needs more people (Dupont & Plummer 1998).

**Clustering**, as in
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md):
with `n_per_site` fixed, \\n = n_0 D\_{eff}\\ with \\D\_{eff} = 1 +
(n\_{per site} - 1) ICC\\; with `n_sites` fixed, \\n = n_0
n\_{sites}(1 - ICC) / (n\_{sites} - n_0 ICC)\\, which exists only when
\\n\_{sites} \> n_0 ICC\\; otherwise the function stops and reports the
highest power those sites can give. If fewer people are needed than
there are sites, fewer sites are used, with 1 person at each, and a
warning says so. With sites, `n_per_timepoint` is the number actually
sampled, `n_sites * n_per_site`.

**Reverse mode.** The same equation run backwards: the power of the
slope test with \\n / D\_{eff}\\ people per timepoint.

## References

Dupont, W. D. & Plummer, W. D. (1998) Power and sample size calculations
for studies involving linear regression. *Controlled Clinical Trials*
**19**(6), 589-601.

Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)

Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
of a screening test. *American Journal of Epidemiology* **107**(1),
71-76.
[doi:10.1093/oxfordjournals.aje.a112510](https://doi.org/10.1093/oxfordjournals.aje.a112510)

MMS-SD Study Design Workshop.
<https://mrc-ide.github.io/MMS-SD_workshop/>

## See also

[`test_trend()`](https://mrc-ide.github.io/MMSpower/reference/test_trend.md)

## Examples

``` r
# Prevalence rising from 10 to 20 percent over 5 annual rounds
design_trend(prevalence_start = 0.10, prevalence_end = 0.20, n_timepoints = 5)
#> $n_per_timepoint
#> [1] 161
#> 
#> $n
#> [1] 805
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
#> [1] 161
#> 
#> $power
#> [1] 0.8
#> 

# The same trend given as a slope
design_trend(prevalence_start = 0.10, slope = 0.025, n_timepoints = 5)
#> $n_per_timepoint
#> [1] 161
#> 
#> $n
#> [1] 805
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
#> [1] 161
#> 
#> $power
#> [1] 0.8
#> 

# Reverse: the power of 150 people per round
design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
             n_per_timepoint = 150)
#> $n_per_timepoint
#> [1] 150
#> 
#> $n
#> [1] 750
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
#> [1] 150
#> 
#> $power
#> [1] 0.7738757
#> 

# Fixed site size: how many sites per round?
design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
             n_per_site = 25, icc = 0.03)
#> $n_per_timepoint
#> [1] 300
#> 
#> $n
#> [1] 1500
#> 
#> $n_sites
#> [1] 12
#> 
#> $n_per_site
#> [1] 25
#> 
#> $deff
#> [1] 1.72
#> 
#> $n_eff
#> [1] 161
#> 
#> $power
#> [1] 0.8
#> 
```
