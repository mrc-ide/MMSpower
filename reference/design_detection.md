# Calculate sample size to detect at least one case

Given the expected true prevalence of a (possibly rare) variant or
condition, returns the number of people needed so that the chance of
finding at least one case reaches a target (`detection_prob`). Supply
`n` instead to run it in reverse: the function then reports the chance
of finding at least one case with that many people.

**Three design modes** – controlled by `n_sites` and `n_per_site`:

- Unclustered (`n_sites = NULL`, `n_per_site = NULL`):

  Returns the number of people.

- Fixed site size (`n_per_site` supplied):

  Solves for the number of sites.

- Fixed number of sites (`n_sites` supplied):

  Solves for the number of people per site.

## Usage

``` r
design_detection(
  prevalence,
  detection_prob = 0.8,
  n = NULL,
  sensitivity = 1,
  n_sites = NULL,
  n_per_site = NULL,
  icc = 0
)
```

## Arguments

- prevalence:

  Expected true prevalence, in (0, 1): the lowest prevalence the study
  should be able to detect.

- detection_prob:

  Target chance of finding at least one case, in (0, 1); default 0.80.
  Not used when `n` is supplied.

- n:

  Optional whole number. If supplied, the function runs in reverse and
  reports the chance of finding at least one case among `n` people.

- sensitivity:

  Diagnostic sensitivity in (0, 1\]; default 1 (the test never misses a
  case).

- n_sites:

  Optional whole number. Fixes the number of sites; the function solves
  for people per site. Cannot be used with `n_per_site`.

- n_per_site:

  Optional whole number. Fixes the number of people per site; the
  function solves for the number of sites. Cannot be used with
  `n_sites`.

- icc:

  Intra-cluster correlation, in \[0, 1); default 0 (no clustering). An
  `icc` above 0 needs `n_sites` or `n_per_site`.

## Value

A list with:

- n:

  Total sample size (or `n` as supplied, in reverse). With sites,
  `n_sites * n_per_site`.

- n_sites:

  Number of sites (`NULL` if unclustered).

- n_per_site:

  People per site (`NULL` if unclustered).

- deff:

  Design effect: `n` divided by the number of people with no clustering
  who give the same chance of finding a case; 1 if unclustered.

- n_eff:

  Number of people needed with no clustering, \\n_0\\ rounded up (in
  reverse, the number with no clustering who give the same chance of
  finding a case).

- detection_prob:

  The target (forward), or the chance achieved with `n` people
  (reverse).

## Details

**Sample size.** Each person tested has a chance \\q = prevalence \times
sensitivity\\ of being a true case whom the test detects. The chance
that no case is found among \\n\\ people is \\(1 - q)^n\\, so the number
of people needed is

\$\$n_0 = \frac{\log(1 - detection\\prob)}{\log(1 - q)}\$\$

rounded up (Cannon 2001, equation 4.6). Positives are assumed to be
confirmed, for example by a second test (Cannon 2001), so a false
positive is not counted as finding a case, and specificity is not used.

**Clustering.** Cases bunch together at sites, so the chance \\q\\ is
taken to vary between sites as a beta distribution with mean \\q\\ and
intra-cluster correlation ICC, as in DRpower. The chance that no case is
found among \\m\\ people at one site is

\$\$P_0(m) = \frac{B(\alpha, \beta + m)}{B(\alpha, \beta)}, \qquad
\alpha = q \left(\frac{1}{ICC} - 1\right), \quad \beta = (1 - q)
\left(\frac{1}{ICC} - 1\right),\$\$

and at \\k\\ sites it is \\P_0(m)^k\\. With `n_per_site` fixed, the
number of sites is \\\lceil \log(1 - detection\\prob) / \log P_0(n\_{per
site}) \rceil\\. With `n_sites` fixed there is no formula for the people
per site, so the smallest whole number that reaches `detection_prob` is
found by a numerical search. \\P_0(m)\\ always falls towards 0 as \\m\\
grows, but very slowly when the ICC is high and cases are rare; if even
\\10^{15}\\ people per site fall short, the function stops.

**Sites and rounding.** With `n_sites` fixed, every site is used, with
at least 1 person each, even if fewer people would do (as in the other
design functions). A site larger than needed gives 1 site of that size.
With sites, the returned `n` is the number actually sampled,
`n_sites * n_per_site`.

## References

Cannon, R. M. (2001) Sense and sensitivity: designing surveys based on
an imperfect test. *Preventive Veterinary Medicine* **49**(3-4),
141-163.
[doi:10.1016/S0167-5877(01)00184-2](https://doi.org/10.1016/S0167-5877%2801%2900184-2)
(Section 4.)

Verity, R. & Ruybal, S. DRpower: Study design and analysis for pfhrp2/3
deletion prevalence studies. R package version 1.0.4 (functions
`get_power_presence()` and `get_sample_size_presence()`).
<https://github.com/mrc-ide/DRpower>

MMS-SD Study Design Workshop, Module 7.
<https://mrc-ide.github.io/MMS-SD_workshop/>

## Examples

``` r
# 80 percent chance of finding a case at 5 percent prevalence
design_detection(prevalence = 0.05)
#> $n
#> [1] 32
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
#> [1] 32
#> 
#> $detection_prob
#> [1] 0.8
#> 

# Imperfect test
design_detection(prevalence = 0.05, sensitivity = 0.8)
#> $n
#> [1] 40
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
#> [1] 40
#> 
#> $detection_prob
#> [1] 0.8
#> 

# Reverse: the chance of finding a 1 percent variant with 139 people
design_detection(prevalence = 0.01, n = 139)
#> $n
#> [1] 139
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
#> [1] 139
#> 
#> $detection_prob
#> [1] 0.7526613
#> 

# Fixed site size: how many sites?
design_detection(prevalence = 0.05, n_per_site = 10, icc = 0.05)
#> $n
#> [1] 40
#> 
#> $n_sites
#> [1] 4
#> 
#> $n_per_site
#> [1] 10
#> 
#> $deff
#> [1] 1.223414
#> 
#> $n_eff
#> [1] 32
#> 
#> $detection_prob
#> [1] 0.8
#> 

# Fixed number of sites: how many people per site?
design_detection(prevalence = 0.01, n_sites = 5, icc = 0.1)
#> $n
#> [1] 1470
#> 
#> $n_sites
#> [1] 5
#> 
#> $n_per_site
#> [1] 294
#> 
#> $deff
#> [1] 9.17942
#> 
#> $n_eff
#> [1] 161
#> 
#> $detection_prob
#> [1] 0.8
#> 
```
