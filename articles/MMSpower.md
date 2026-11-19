# Get started with MMSpower

``` r

library(MMSpower)
```

MMSpower answers four common questions in malaria molecular surveillance
studies. For most, a **design** function says how many samples are
needed, and a matching **analysis** function says what the data show.
Every function allows for an imperfect diagnostic test (`sensitivity`,
and `specificity` where it applies) and for clustered sampling across
sites (`icc`, with `n_sites` or `n_per_site`).

This page is a short tour. Each function’s reference page explains its
methods and formulas in full.

## Estimating prevalence

How many people are needed to estimate a prevalence of around 30% to
within ±5 percentage points, with 10 people per site and an ICC of 0.05?

``` r

design_precision(prevalence = 0.3, moe = 0.05, n_per_site = 10, icc = 0.05)
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
```

Once the data are in, estimate the prevalence and its confidence
interval:

``` r

estimate_prevalence(x = rep(3, 47), n = rep(10, 47), icc = 0.05)
#> $prevalence
#> [1] 0.3
#> 
#> $ci_lower
#> [1] 0.2506
#> 
#> $ci_upper
#> [1] 0.353099
#> 
#> $moe
#> [1] 0.05124951
#> 
#> $n_eff
#> [1] 324.1379
#> 
#> $icc
#> [1] 0.05
#> 
#> $deff
#> [1] 1.45
#> 
#> $df
#> [1] Inf
```

## Testing against a threshold

How many people are needed to show that prevalence is above a 5%
threshold when it is really 10%, with 80% power?

``` r

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
```

Then test the data against the threshold:

``` r

test_threshold(x = 18, n = 191, threshold = 0.05)
#> $p_value
#> [1] 0.007726515
#> 
#> $reject
#> [1] TRUE
#> 
#> $prevalence
#> [1] 0.09424084
#> 
#> $ci_lower
#> [1] 0.06179912
#> 
#> $ci_upper
#> [1] 1
#> 
#> $n_eff
#> [1] 191
#> 
#> $icc
#> [1] 0
#> 
#> $deff
#> [1] 1
#> 
#> $df
#> [1] Inf
```

## Detecting a rare variant

How many people are needed for an 80% chance of finding at least one
case of a variant present in 1% of infections?

``` r

design_detection(prevalence = 0.01)
#> $n
#> [1] 161
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
#> $detection_prob
#> [1] 0.8
```

## Comparing two groups

How many people per group are needed to detect a difference between 10%
and 20% prevalence, with 30 sites per group and an ICC of 0.02?

``` r

design_difference(prevalence1 = 0.10, prevalence2 = 0.20, n_sites = 30, icc = 0.02)
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
```

Then test for a difference:

``` r

test_difference(x1 = 24, n1 = 240, x2 = 48, n2 = 240, icc = 0.02)
#> $difference
#> [1] -0.1
#> 
#> $ci_lower
#> [1] -0.2611286
#> 
#> $ci_upper
#> [1] 0.05779169
#> 
#> $p_value
#> [1] 0.2019338
#> 
#> $reject
#> [1] FALSE
#> 
#> $prevalence1
#> [1] 0.1
#> 
#> $prevalence2
#> [1] 0.2
#> 
#> $n_eff1
#> [1] 41.52249
#> 
#> $n_eff2
#> [1] 41.52249
#> 
#> $icc1
#> [1] 0.02
#> 
#> $icc2
#> [1] 0.02
#> 
#> $deff1
#> [1] 5.78
#> 
#> $deff2
#> [1] 5.78
#> 
#> $df
#> [1] Inf
```

## Costing a study

Finally, cost a design from its sample size and number of sites:

``` r

estimate_cost(n = 470, cost_per_sample = 50, n_sites = 47,
              fixed_cost_per_site = 500, transport_cost_per_site = 200)
#> $total_cost
#> [1] 56400
#> 
#> $variable_cost
#> [1] 23500
#> 
#> $site_cost
#> [1] 32900
#> 
#> $budget_remaining
#> NULL
```
