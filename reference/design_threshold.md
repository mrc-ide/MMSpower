# Calculate sample size to detect a prevalence above or below a threshold

Design function: given a decision threshold and the true prevalence
expected in the population, returns the number of people needed for
[`test_threshold()`](https://mrc-ide.github.io/MMSpower/reference/test_threshold.md)
(with its default exact binomial test) to show that the prevalence is
above (or below) the threshold with the target power.

Three designs are supported:

- No clustering (`n_sites` and `n_per_site` both `NULL`, default):

  Returns the number of people.

- Fixed number of people per site (`n_per_site`):

  Solves for the number of sites needed.

- Fixed number of sites (`n_sites`):

  Solves for the number of people needed per site.

## Usage

``` r
design_threshold(
  threshold,
  prevalence,
  power = 0.8,
  alternative = c("greater", "less", "two.sided"),
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  n_sites = NULL,
  n_per_site = NULL,
  icc = 0
)
```

## Arguments

- threshold:

  Decision threshold for the true prevalence, strictly between 0 and 1
  (for example 0.05 for the WHO pfhrp2/3 deletion threshold). This is
  the value
  [`test_threshold()`](https://mrc-ide.github.io/MMSpower/reference/test_threshold.md)
  tests against.

- prevalence:

  The true prevalence expected in the population, in (0, 1). It must be
  above `threshold` for `alternative = "greater"` and below it for
  `"less"`.

- power:

  Target power, in (0, 1); default 0.80.

- alternative:

  Which side of the threshold to look for, the same as in the
  [`test_threshold()`](https://mrc-ide.github.io/MMSpower/reference/test_threshold.md)
  call. `"greater"` (default) looks for a prevalence above the
  threshold, `"less"` for one below it, and `"two.sided"` for either.

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1 (the test never misses a
  case).

- specificity:

  Diagnostic specificity in (0, 1\]; default 1 (the test never gives a
  false positive). `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level, in (0, 1); default 0.95. The significance level of
  the test is 1 - `conf_level`.

- n_sites:

  Optional whole number. Fixes the number of sites; the function solves
  for the number of people per site. Cannot be used with `n_per_site`.

- n_per_site:

  Optional whole number. Fixes the number of people per site; the
  function solves for the number of sites. Cannot be used with
  `n_sites`.

- icc:

  Intra-cluster correlation, in \[0, 1\]; default 0 (no clustering). An
  `icc` above 0 needs `n_sites` or `n_per_site`.

## Value

A list with (every value `NA`, with a warning, if more than 100,000
people are needed without clustering):

- n:

  Total number of people needed. With sites, `n_sites * n_per_site`.

- n_sites:

  Number of sites (`NULL` with no clustering).

- n_per_site:

  People per site (`NULL` with no clustering).

- deff:

  Design effect of the returned design, `1 + (n_per_site - 1) * icc` (1
  with no clustering).

- n_eff:

  Number of independent people needed: the exact answer with no
  clustering.

## Details

**Exact power.** The power of a design is the chance that the exact
binomial test in
[`test_threshold()`](https://mrc-ide.github.io/MMSpower/reference/test_threshold.md)
rejects when the true prevalence is `prevalence`. For \\n\\ people and
`alternative = "greater"`, the test rejects when the number of positives
\\X\\ is at least the critical count \\k\\, the smallest count with
\\P(X \ge k \mid \theta\_{app}) \< \alpha\\; the power is then \\P(X \ge
k \mid p\_{app})\\. `"less"` is the mirror image, and `"two.sided"` adds
the two tails, each at \\\alpha / 2\\. This is worked out from the
binomial distribution itself, not from a normal approximation, which can
undersize the study badly at low prevalence.

**Saw-tooth power and the choice of n.** With an exact test, power does
not always go up when people are added: counts are whole numbers, so the
critical count jumps, and power goes up and down in a saw-tooth. For
example, with `threshold = 0.10`, `prevalence = 0.05`,
`alternative = "less"` and 95 percent confidence, 179 people give 81.3
percent power, 182 to 190 people all give less than 80 percent, and
every number from 203 on gives at least 80 percent. The function returns
the larger number (203): the smallest \\n\\ from which every larger
\\n\\ also reaches the target power. This is the safe choice, as a few
extra (or a few lost) people never drop the study below the target.
Sample sizes up to 100,000 are checked. A design needing more than that
without clustering (a prevalence very close to the threshold) gives a
warning that reports the power 100,000 people would reach, and returns
`NA` for every value.

**Imperfect test.** The test is worked out on the measured scale: a true
prevalence \\p\\ is measured as \\p\_{app} = p \cdot Se + (1-p)(1-Sp)\\,
and the same for the threshold. The measured gap is the true gap times
\\Se + Sp - 1\\, so an imperfect test needs more people.

**Clustering.** The exact search above gives \\n_0\\, the number of
independent people needed (returned as `n_eff`). Clustering is then
allowed for with the design effect, as in
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md):
with `n_per_site` fixed, \\n = n_0 \cdot D\_{eff}\\ with \\D\_{eff} =
1 + (n\_{per site} - 1) ICC\\; with `n_sites` fixed,

\$\$n = \frac{n_0 \cdot n\_{sites} \cdot (1-ICC)}{n\_{sites} - n_0 \cdot
ICC}\$\$

which exists only when \\n\_{sites} \> n_0 \cdot ICC\\ (each site is
worth at most \\1/ICC\\ independent people, however many are sampled
there); otherwise the function stops. The number of sites (or people per
site) is then rounded up, and `n` is `n_sites * n_per_site`. Scaling the
exact unclustered answer by the design effect is an approximation: the
effective sample size \\n / D\_{eff}\\ of the returned design is at
least \\n_0\\, but need not sit on the saw-tooth exactly where \\n_0\\
does. The power assumes the analysis uses the same ICC.

## References

MMS-SD Study Design Workshop.
<https://mrc-ide.github.io/MMS-SD_workshop/>

Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)

Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
of a screening test. *American Journal of Epidemiology* **107**(1),
71-76.
[doi:10.1093/oxfordjournals.aje.a112510](https://doi.org/10.1093/oxfordjournals.aje.a112510)

## See also

[`test_threshold`](https://mrc-ide.github.io/MMSpower/reference/test_threshold.md)

## Examples

``` r
# Show that prevalence is above a 5 percent threshold when it is really
# 10 percent, with 80 percent power
design_threshold(threshold = 0.05, prevalence = 0.10)
#> $n
#> [1] 191
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
#> [1] 191
#> 

# Higher power needs more people
design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.90)
#> $n
#> [1] 256
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
#> [1] 256
#> 

# Below the threshold instead
design_threshold(threshold = 0.10, prevalence = 0.05, alternative = "less")
#> $n
#> [1] 203
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
#> [1] 203
#> 

# Fixed number of people per site: how many sites?
design_threshold(threshold = 0.05, prevalence = 0.10,
                 n_per_site = 20, icc = 0.05)
#> $n
#> [1] 380
#> 
#> $n_sites
#> [1] 19
#> 
#> $n_per_site
#> [1] 20
#> 
#> $deff
#> [1] 1.95
#> 
#> $n_eff
#> [1] 191
#> 
```
