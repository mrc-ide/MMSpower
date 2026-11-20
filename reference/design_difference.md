# Calculate sample size to detect a difference in prevalence between two groups

Design function: given the expected true prevalence in two groups,
returns the sample size needed in each group for
[`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md)
(with its default score test) to detect that difference with the target
power; or, given a full design for both groups, returns the power it
gives.

Each group has its own arguments, `n1`, `n_sites1`, `n_per_site1` for
group 1 and `n2`, `n_sites2`, `n_per_site2` for group 2, as in
[`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md).
Each group-2 argument defaults to its group-1 value, so for two groups
with the same design only the group-1 arguments are needed.

**What is solved for** – leave out the one quantity you want, the same
in both groups:

- No sites:

  Returns the number of people per group.

- `n_per_site1` (and `n_per_site2`):

  Returns the number of sites per group, the same in both.

- `n_sites1` (and `n_sites2`):

  Returns the number of people per site, the same in both.

- A full design for both groups:

  Returns the power. For each group give `n` (no sites), `n_sites` and
  `n_per_site`, or a vector `n_per_site` with one size per site (for
  sites of different sizes). Leave `power` at its default.

When solving, the groups' sizes follow their fixed values: for example,
`n_sites1 = 10, n_sites2 = 20` gives group 2 twice as many people as
group 1.

## Usage

``` r
design_difference(
  prevalence1,
  prevalence2,
  power = 0.8,
  alternative = c("two.sided", "greater", "less"),
  sensitivity = 1,
  specificity = 1,
  conf_level = 0.95,
  n1 = NULL,
  n2 = n1,
  n_sites1 = NULL,
  n_sites2 = n_sites1,
  n_per_site1 = NULL,
  n_per_site2 = n_per_site1,
  icc = 0
)
```

## Arguments

- prevalence1:

  Expected true prevalence in group 1, in \[0, 1\].

- prevalence2:

  Expected true prevalence in group 2, in \[0, 1\]. The two prevalences
  must differ.

- power:

  Target power, in \[0.5, 1); default 0.80. Not used, and not to be
  given, with a full design.

- alternative:

  Which kind of difference to look for. `"two.sided"` (default) looks
  for any difference; `"greater"` looks for group 1 being higher
  (requires `prevalence1 > prevalence2`); `"less"` looks for group 1
  being lower (requires `prevalence1 < prevalence2`).

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1. Assumed the same in both
  groups.

- specificity:

  Diagnostic specificity in (0, 1\]; default 1. Assumed the same in both
  groups. `sensitivity + specificity` must exceed 1.

- conf_level:

  Confidence level, in \[0.5, 1); default 0.95. The significance level
  of the test is 1 - `conf_level`. Below 0.5 the one-sided z value turns
  negative and the sample size shrinks towards nonsense, so it is not
  allowed.

- n1, n2:

  Optional whole numbers: the number of people in each group, for a full
  design with no sites. `n2` defaults to `n1`.

- n_sites1, n_sites2:

  Optional whole numbers: the number of sites in each group. Alone, the
  function solves for the number of people per site. `n_sites2` defaults
  to `n_sites1`.

- n_per_site1, n_per_site2:

  Optional: people per site in each group. One whole number (alone, the
  function solves for the number of sites), or a vector with one size
  per site, a full design with sites of different sizes (then leave out
  that group's `n_sites`). `n_per_site2` defaults to `n_per_site1`.

- icc:

  Intra-cluster correlation, in \[0, 1\]: one number for both groups, or
  two, `c(group 1, group 2)`. Default 0 (no clustering). An `icc` above
  0 needs sites in that group.

## Value

A list with:

- n_per_group:

  People in each group, `c(group 1, group 2)`. With sites,
  `n_sites * n_per_site`.

- n:

  Total number of people, `sum(n_per_group)`.

- n_sites:

  Sites in each group (`NULL` with no clustering).

- n_per_site:

  People per site in each group (`NULL` with no clustering); a list of
  each group's site sizes if they differ within a group.

- deff:

  Design effect of each group (1 with no clustering).

- n_eff:

  People each group would need with no clustering, keeping the groups in
  the same ratio as the design; for a full design, `n_per_group / deff`.

- power:

  Power of the returned design.

## Details

**Sample size.** Write \\m_1\\, \\m_2\\ for the effective sizes of the
two groups (below). The score test in
[`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md)
reaches the target power when

\$\$z\_\alpha \sqrt{\bar{p}(1-\bar{p})\\(1/m_1 + 1/m_2)} + z\_\beta
\sqrt{p\_{1,app}(1-p\_{1,app})/m_1 + p\_{2,app}(1-p\_{2,app})/m_2} =
\|p\_{1,app} - p\_{2,app}\|\$\$

where \\p\_{g,app} = p_g \cdot Se + (1-p_g)(1-Sp)\\ is the prevalence
the test would measure in group \\g\\, \\\bar{p} = (p\_{1,app} m_1 +
p\_{2,app} m_2)/(m_1 + m_2)\\ is the pooled prevalence the test uses,
\\z\_\beta\\ = `qnorm(power)` (0.84 for 80 percent power), and
\\z\_\alpha\\ = `qnorm(1 - alpha/2)` for a two-sided test (1.96 at 95
percent) or `qnorm(1 - alpha)` for a one-sided test (1.64 at 95
percent), with \\\alpha = 1 -\\ `conf_level`. With no clustering and two
groups of \\n\\ people this is the formula of R's
[`power.prop.test()`](https://rdrr.io/r/stats/power.prop.test.html).

Because \\p\_{1,app} - p\_{2,app} = (p_1 - p_2)(Se + Sp - 1)\\, an
imperfect test shrinks the difference the test can see and increases the
sample size.

**Effective size.** For group \\g\\ with \\n_g\\ people in sites of
\\s_g\\ people, \\m_g = n_g / D_g\\ with design effect \\D_g = 1 +
(s_g - 1)\\ICC_g\\. For sites of different sizes, \\s_g\\ is the
size-weighted mean site size \\\sum s^2 / \sum s\\, as in
[`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md).
With one unknown \\x\\, each design gives \\m_g\\ directly:

- No clustering:

  \\m_g = x\\, people per group.

- Fixed `n_per_site1`, `n_per_site2`:

  \\m_g = x \\ s_g / D_g\\, with \\x\\ sites per group.

- Fixed `n_sites1`, `n_sites2`:

  \\m_g = k_g x / (1 + (x - 1) ICC_g)\\, with \\k_g\\ sites and \\x\\
  people per site.

Each \\m_g\\ grows with \\x\\, so the left side of the equation above
falls as \\x\\ grows, and the function finds the smallest \\x\\ that
meets it with a single root search
([`uniroot()`](https://rdrr.io/r/stats/uniroot.html)), then rounds up to
a whole number.

**Power.** For a full design, or the rounded design that was solved for,
the power is the equation above solved for \\z\_\beta\\:
\\\Phi((\|p\_{1,app} - p\_{2,app}\| - z\_\alpha \sqrt{\bar{p}(1-\bar{p})
(1/m_1 + 1/m_2)}) / \sqrt{p\_{1,app}(1-p\_{1,app})/m_1 +
p\_{2,app}(1-p\_{2,app})/m_2})\\.

**Highest power.** With `n_sites` fixed, adding people to each site
cannot raise the power past a limit: each site is worth at most
\\1/ICC\\ independent people, so group \\g\\ never counts as more than
\\k_g / ICC_g\\ independent people. For example, with ICC 0.05 each site
is worth at most 20 independent people, so 8 sites per group give at
most 160 independent people; for 10 versus 20 percent that gives 70.9
percent power, so 80 percent can never be reached. When the target is
out of reach, the function stops and reports the power at that limit.

**Rounding.** The number of sites, or people per site, is rounded up to
a whole number, the same in both groups. With `n_sites` fixed, every
site is used, with at least 1 person each, even if fewer people would do
(as in the other design functions); with `n_per_site` larger than
needed, 1 site per group is used. The people in each group are
`n_sites * n_per_site`, and the design effects are those of the rounded
design.

## References

MMS-SD Study Design Workshop.
<https://mrc-ide.github.io/MMS-SD_workshop/>

Fleiss, J. L., Levin, B. & Paik, M. C. (2003) *Statistical Methods for
Rates and Proportions*, 3rd ed. Wiley. (Two-group sample size.)

Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)

Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
of a screening test. *American Journal of Epidemiology* **107**(1),
71-76.
[doi:10.1093/oxfordjournals.aje.a112510](https://doi.org/10.1093/oxfordjournals.aje.a112510)

## See also

[`test_difference`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md)

## Examples

``` r
# Detect 10 versus 20 percent prevalence with 80 percent power
design_difference(prevalence1 = 0.10, prevalence2 = 0.20)
#> $n_per_group
#> [1] 199 199
#> 
#> $n
#> [1] 398
#> 
#> $n_sites
#> NULL
#> 
#> $n_per_site
#> NULL
#> 
#> $deff
#> [1] 1 1
#> 
#> $n_eff
#> [1] 199 199
#> 
#> $power
#> [1] 0.8000726
#> 

# Higher power needs more people
design_difference(0.10, 0.20, power = 0.90)
#> $n_per_group
#> [1] 266 266
#> 
#> $n
#> [1] 532
#> 
#> $n_sites
#> NULL
#> 
#> $n_per_site
#> NULL
#> 
#> $deff
#> [1] 1 1
#> 
#> $n_eff
#> [1] 266 266
#> 
#> $power
#> [1] 0.9001549
#> 

# Imperfect diagnostic test
design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95)
#> $n_per_group
#> [1] 317 317
#> 
#> $n
#> [1] 634
#> 
#> $n_sites
#> NULL
#> 
#> $n_per_site
#> NULL
#> 
#> $deff
#> [1] 1 1
#> 
#> $n_eff
#> [1] 317 317
#> 
#> $power
#> [1] 0.8012193
#> 

# Fixed number of people per site: how many sites per group?
design_difference(0.10, 0.20, n_per_site1 = 30, icc = 0.02)
#> $n_per_group
#> [1] 330 330
#> 
#> $n
#> [1] 660
#> 
#> $n_sites
#> [1] 11 11
#> 
#> $n_per_site
#> [1] 30 30
#> 
#> $deff
#> [1] 1.58 1.58
#> 
#> $n_eff
#> [1] 199 199
#> 
#> $power
#> [1] 0.8188356
#> 

# Fixed number of sites: how many people per site?
design_difference(0.10, 0.20, n_sites1 = 30, icc = 0.02)
#> $n_per_group
#> [1] 240 240
#> 
#> $n
#> [1] 480
#> 
#> $n_sites
#> [1] 30 30
#> 
#> $n_per_site
#> [1] 8 8
#> 
#> $deff
#> [1] 1.14 1.14
#> 
#> $n_eff
#> [1] 199 199
#> 
#> $power
#> [1] 0.82185
#> 

# Twice as many sites in group 2, and more clustering there
design_difference(0.10, 0.20, n_sites1 = 20, n_sites2 = 40,
                  icc = c(0.02, 0.05))
#> $n_per_group
#> [1] 200 400
#> 
#> $n
#> [1] 600
#> 
#> $n_sites
#> [1] 20 40
#> 
#> $n_per_site
#> [1] 10 10
#> 
#> $deff
#> [1] 1.18 1.45
#> 
#> $n_eff
#> [1] 155 309
#> 
#> $power
#> [1] 0.8121849
#> 

# Full design: the power of 200 people per group
design_difference(0.10, 0.20, n1 = 200)
#> $n_per_group
#> [1] 200 200
#> 
#> $n
#> [1] 400
#> 
#> $n_sites
#> NULL
#> 
#> $n_per_site
#> NULL
#> 
#> $deff
#> [1] 1 1
#> 
#> $n_eff
#> [1] 200 200
#> 
#> $power
#> [1] 0.8020484
#> 

# Full design with sites of different sizes in each group
design_difference(0.10, 0.20,
                  n_per_site1 = c(80, 60, 40, 40), n_per_site2 = c(100, 50, 50),
                  icc = 0.02)
#> $n_per_group
#> [1] 220 200
#> 
#> $n
#> [1] 420
#> 
#> $n_sites
#> [1] 4 3
#> 
#> $n_per_site
#> $n_per_site[[1]]
#> [1] 80 60 40 40
#> 
#> $n_per_site[[2]]
#> [1] 100  50  50
#> 
#> 
#> $deff
#> [1] 2.18 2.48
#> 
#> $n_eff
#> [1] 100.91743  80.64516
#> 
#> $power
#> [1] 0.4784116
#> 
```
