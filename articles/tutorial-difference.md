# Comparing two groups

``` r

library(MMSpower)
```

A common question in molecular surveillance is whether the prevalence of
a marker differs between two groups: two regions, or the same region
before and after an intervention. This article plans such a comparison
with
[`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md)
and then analyses the data with
[`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md).

Throughout, the scenario is a drug-resistance marker expected at around
10% in region A (group 1) and 20% in region B (group 2). We want 80%
power to detect that difference at the 5% significance level.

## How many people, with no clustering?

Start with the simplest case: everyone sampled is independent and the
diagnostic test is perfect.

``` r

des <- design_difference(prevalence1 = 0.10, prevalence2 = 0.20)
des
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
```

The answer is given per group (`n_per_group`, one value for each group)
and as a total (`n`). With no clustering this is the same calculation as
R’s [`power.prop.test()`](https://rdrr.io/r/stats/power.prop.test.html),
which gives a fractional answer that
[`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md)
rounds up:

``` r

power.prop.test(p1 = 0.10, p2 = 0.20, power = 0.80)$n
#> [1] 198.9634
```

## One-sided or two-sided?

The default, `alternative = "two.sided"`, looks for a difference in
either direction. If the only result that would change what you do is
region A being *lower* than region B (or, after an intervention,
prevalence going down), a one-sided test is enough, and it needs fewer
people:

``` r

design_difference(0.10, 0.20, alternative = "less")$n_per_group
#> [1] 157 157
```

`"less"` means group 1 lower than group 2; `"greater"` means group 1
higher. Decide the direction before seeing the data. A one-sided test
cannot detect a change in the other direction, however large.

## An imperfect test

If the marker assay misses some true positives (sensitivity below 1) or
produces some false positives (specificity below 1), the measured
prevalences move towards each other. The measured difference is the true
difference multiplied by `sensitivity + specificity - 1`, so a smaller
difference has to be detected and more people are needed:

``` r

design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95)$n_per_group
#> [1] 317 317
```

Sensitivity and specificity are assumed to be the same in both groups.
That is reasonable when the same laboratory and assay are used for both,
and worth checking if they are not.

## Clustered sampling

In practice people are sampled at sites (health facilities, villages or
schools), and people at the same site tend to be more alike than people
at different sites. The intra-cluster correlation (ICC) measures this.
The design effect, `1 + (people per site - 1) * ICC`, says how many
times more people a clustered design needs than an unclustered one.

There are two ways to plan a clustered study. If you know how many
people you can sample at each site, fix `n_per_site` and the function
solves for the number of sites per group:

``` r

design_difference(0.10, 0.20, n_per_site = 20, icc = 0.05)
#> $n_per_group
#> [1] 400 400
#> 
#> $n
#> [1] 800
#> 
#> $n_sites
#> [1] 20 20
#> 
#> $n_per_site
#> [1] 20 20
#> 
#> $deff
#> [1] 1.95 1.95
#> 
#> $n_eff
#> [1] 199 199
```

If the number of sites is fixed instead (for example, the health
facilities already in a surveillance network), fix `n_sites` and the
function solves for the number of people per site:

``` r

design_difference(0.10, 0.20, n_sites = 20, icc = 0.05)
#> $n_per_group
#> [1] 380 380
#> 
#> $n
#> [1] 760
#> 
#> $n_sites
#> [1] 20 20
#> 
#> $n_per_site
#> [1] 19 19
#> 
#> $deff
#> [1] 1.9 1.9
#> 
#> $n_eff
#> [1] 199 199
```

In both cases `n_eff` is what each group would need without clustering,
and `deff` is the design effect that turns one into the other.

## Groups that differ

The two groups need not share a design. Suppose region B has twice as
many sites in the network, and earlier surveys suggest more clustering
there. Give `n_sites` and `icc` as two values, `c(group 1, group 2)`:

``` r

diff_des <- design_difference(0.10, 0.20, n_sites = c(10, 20), icc = c(0.02, 0.05))
diff_des
#> $n_per_group
#> [1] 270 540
#> 
#> $n
#> [1] 810
#> 
#> $n_sites
#> [1] 10 20
#> 
#> $n_per_site
#> [1] 27 27
#> 
#> $deff
#> [1] 1.52 2.30
#> 
#> $n_eff
#> [1] 155 309
```

Every element of the result is now a pair of values, one per group. The
people per site is the same in both groups, so region B, with twice the
sites, has twice the people. With the same people per site, its larger
design effect comes from its higher ICC alone, and its `n_eff` is about
twice region A’s, in the same ratio as the design.

## When too few sites cannot reach the target

With the number of sites fixed, adding people to each site helps less
and less. Once a site has many people, another person there tells you
little that the others have not, because they share the same local
conditions. However many are sampled, each site is worth at most
`1 / ICC` independent people: with an ICC of 0.05, that is 20 people. So
8 sites per group can never give more than 160 independent people per
group, which is not enough for 80% power here. The function says so, and
reports the highest power that design can reach:

``` r

design_difference(0.10, 0.20, n_sites = 8, icc = 0.05)
#> Error:
#> ! Target power of 80% cannot be reached with these sites and ICC; the highest power they can give is 70.9%. Increase `n_sites` or lower the ICC.
```

Just above the limit, the answer is possible but not sensible:

``` r

design_difference(0.10, 0.20, n_sites = 10, icc = 0.05)$n_per_site
#> [1] 3647 3647
```

Ten sites can give at most 200 independent people per group, barely
above the 199 needed, so thousands of people per site are required. If
the number of people per site comes out implausibly large, the fix is
more sites, not more people.

## Analysing the data

Once the data are in,
[`test_difference()`](https://mrc-ide.github.io/MMSpower/reference/test_difference.md)
tests whether the prevalences differ. It takes the number of positives
(`x1`, `x2`) and the number tested (`n1`, `n2`) at each site in each
group.

The data below are **simulated**, not real: 20 sites per group with 19
people each (the `n_sites = 20` design above), true prevalences of 10%
and 20%, and an ICC of 0.05.

``` r

set.seed(2026)
sim_sites <- function(k, size, prev, icc) {
  # site prevalences vary around prev, with the given ICC (beta-binomial)
  a <- prev * (1 / icc - 1)
  b <- (1 - prev) * (1 / icc - 1)
  rbinom(k, size, rbeta(k, a, b))
}
x1 <- sim_sites(20, 19, 0.10, 0.05)
x2 <- sim_sites(20, 19, 0.20, 0.05)
n1 <- n2 <- rep(19, 20)
x1
#>  [1] 2 2 2 1 0 0 2 0 2 0 4 1 4 3 0 7 4 0 3 2
x2
#>  [1] 3 5 2 5 4 9 3 1 3 0 2 3 1 7 3 4 5 3 3 0
```

``` r

res <- test_difference(x1, n1, x2, n2)
res
#> $difference
#> [1] -0.07105263
#> 
#> $ci_lower
#> [1] -0.1408457
#> 
#> $ci_upper
#> [1] -0.001445529
#> 
#> $p_value
#> [1] 0.04579302
#> 
#> $reject
#> [1] TRUE
#> 
#> $prevalence1
#> [1] 0.1026316
#> 
#> $prevalence2
#> [1] 0.1736842
#> 
#> $n_eff1
#> [1] 194.5196
#> 
#> $n_eff2
#> [1] 213.5336
#> 
#> $icc1
#> [1] 0.05297391
#> 
#> $icc2
#> [1] 0.04330996
#> 
#> $deff1
#> [1] 1.95353
#> 
#> $deff2
#> [1] 1.779579
#> 
#> $df
#> [1] 36.89077
```

The estimated difference, group 1 minus group 2, is -7.1 percentage
points, with a 95% confidence interval from -14.1 to -0.1. Since no ICC
was given, it was estimated from each group’s own sites (`icc1`, `icc2`)
and used to work out the design effects. An estimated ICC is uncertain,
so the test allows for it by using a t distribution in place of the
normal, with `df` degrees of freedom (about 37 here, worked out from the
number of sites in each group). With few sites per group this makes the
test noticeably more cautious; if a trustworthy value from earlier
surveys is available, pass it with `icc`.

### Score or Wald?

The default, `method = "score"`, is the test that
[`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md)
sizes for. `method = "wald"` is the simpler textbook test, which works
out the standard error from each group’s own prevalence. In large
samples the two agree closely:

``` r

wald <- test_difference(x1, n1, x2, n2, method = "wald")
rbind(score = unlist(res[c("difference", "ci_lower", "ci_upper", "p_value")]),
      wald  = unlist(wald[c("difference", "ci_lower", "ci_upper", "p_value")]))
#>        difference   ci_lower     ci_upper    p_value
#> score -0.07105263 -0.1408457 -0.001445529 0.04579302
#> wald  -0.07105263 -0.1396387 -0.002466565 0.04269457
```

In small groups, or with prevalences near 0 or 1, the Wald interval can
be too narrow, so the score method is the safer default. With either
method, the interval leaves out 0 exactly when `reject` is `TRUE`.

### A one-sided test

If the study was planned as one-sided (region A lower), analyse it that
way. The interval is then one-sided too: it runs from -1 up to an upper
bound.

``` r

one <- test_difference(x1, n1, x2, n2, alternative = "less")
unlist(one[c("difference", "ci_lower", "ci_upper", "p_value")])
#>  difference    ci_lower    ci_upper     p_value 
#> -0.07105263 -1.00000000 -0.01338234  0.02289651
```

### An imperfect test

Correcting for sensitivity and specificity changes the estimated
prevalences, the difference and its interval, but not the p-value:

``` r

imp <- test_difference(x1, n1, x2, n2, sensitivity = 0.9, specificity = 0.95)
rbind(perfect   = unlist(res[c("difference", "ci_lower", "ci_upper", "p_value")]),
      imperfect = unlist(imp[c("difference", "ci_lower", "ci_upper", "p_value")]))
#>            difference   ci_lower     ci_upper    p_value
#> perfect   -0.07105263 -0.1408457 -0.001445529 0.04579302
#> imperfect -0.08359133 -0.1657008 -0.001700622 0.04579302
```

This is expected. With the same test in both groups, the measured
prevalences differ exactly when the true ones do, so the evidence for a
difference is the same either way. What changes is the size of the true
difference that the measured one implies: the measured difference is
divided by `sensitivity + specificity - 1`, which makes the difference
larger and its interval wider.
