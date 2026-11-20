# Estimating prevalence

``` r

library(MMSpower)
```

A national malaria programme wants to know how common a drug-resistance
marker (say, a validated *kelch13* mutation) is in one region. Earlier
reports suggest about 20% of infections carry it, and the programme
wants to estimate the prevalence to within ±5 percentage points, with
95% confidence.

Two functions do the work:
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
says how many samples to collect, and
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)
analyses them once they are in.

## Part 1: Designing the study

### How many samples, ignoring clustering?

Start with the simplest case: every sample is independent (as if drawn
by simple random sampling across the region) and the test is perfect.

``` r

d0 <- design_precision(prevalence = 0.2, moe = 0.05)
d0
#> $n
#> [1] 264
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
#> [1] 264
#> 
#> $moe
#> [1] 0.04996977
```

So 264 samples are needed. By default
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
sizes the study for a Clopper-Pearson interval, the same default as
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md),
so the analysis gives back the margin of error that was planned. There
is no formula for this, so the function finds the sample size by a
numerical search. The familiar $`z^2 p(1-p) / \text{MOE}^2`$ formula,
rounded up, is the Wald version, which you get with
`ci_method = "wald"`:

``` r

design_precision(prevalence = 0.2, moe = 0.05, ci_method = "wald")$n
#> [1] 246
```

This is the number most textbooks and the MMS-SD workshop give. The
Clopper-Pearson interval is a little wider, so it needs a few more
samples. The result also gives `moe`, the margin of error the rounded
design actually achieves, here 0.04997, just inside the target. Two
things drive it: the margin of error (halving it roughly quadruples the
sample size) and the expected prevalence (the sample size is largest at
50%). If the expected prevalence is uncertain, using a value closer to
50% is the cautious choice.

### An imperfect test

Molecular assays are good but not perfect. Suppose the assay has
sensitivity 90% and specificity 95%: it misses some true positives, and
some negatives are wrongly called positive.

``` r

d_test <- design_precision(prevalence = 0.2, moe = 0.05,
                           sensitivity = 0.9, specificity = 0.95)
d_test$n
#> [1] 387
```

The sample size rises from 264 to 387. The reason is that the analysis
will correct the measured prevalence back to the true scale (the
Rogan-Gladen correction), and this correction stretches the confidence
interval by a factor of $`1 / (Se + Sp - 1)`$. To keep the same margin
of error on the true scale, more samples are needed. Specificity matters
more than its size suggests: at low prevalence, false positives can make
up a large share of the positives.

### Clustering: samples come from sites

In practice samples are collected at a number of sites (health
facilities, villages), not drawn at random from the whole region. People
at the same site tend to be more alike than people at different sites,
so each extra sample from a site adds less information than a sample
from a new site. The intra-cluster correlation (ICC) measures how strong
this is, and the design effect
$`D_{eff} = 1 + (n_{per\,site} - 1)\,\text{ICC}`$ says how much the
sample size has to grow. For this tutorial we assume an ICC of 0.05, and
go back to a perfect test to keep the numbers simple.

There are two ways to fix the design, depending on what is constrained
in the field.

**Fixing people per site.** If each site can realistically collect 20
samples, solve for the number of sites:

``` r

d_per <- design_precision(prevalence = 0.2, moe = 0.05, n_per_site = 20, icc = 0.05)
d_per
#> $n
#> [1] 520
#> 
#> $n_sites
#> [1] 26
#> 
#> $n_per_site
#> [1] 20
#> 
#> $deff
#> [1] 1.95
#> 
#> $n_eff
#> [1] 264
#> 
#> $moe
#> [1] 0.04971169
```

The study needs 26 sites of 20 people, 520 samples in total.

**Fixing the number of sites.** If the programme already has 20 sentinel
sites, solve for the number of people per site:

``` r

d_sites <- design_precision(prevalence = 0.2, moe = 0.05, n_sites = 20, icc = 0.05)
d_sites
#> $n
#> [1] 740
#> 
#> $n_sites
#> [1] 20
#> 
#> $n_per_site
#> [1] 37
#> 
#> $deff
#> [1] 2.8
#> 
#> $n_eff
#> [1] 264
#> 
#> $moe
#> [1] 0.04994193
```

Each of the 20 sites now needs 37 people, 740 samples in total.

In both cases `n_eff` is the unclustered answer from the first section:
the clustered design carries the same information as 264 independent
samples. Fixing the sites at 20 needs more samples in total than fixing
20 people per site, because the people per site go up, and so does the
design effect.

**Too few sites.** With the number of sites fixed, there is a limit.
Each site is worth at most $`1 / \text{ICC}`$ independent samples,
however many people are tested there, so a small number of sites can
never reach the target precision:

``` r

design_precision(prevalence = 0.2, moe = 0.05, n_sites = 10, icc = 0.05)
#> Error:
#> ! Target MOE of 5.0% is unachievable with 10 sites and ICC = 0.050; the minimum achievable MOE is 5.8%.
```

The error reports the best margin of error that 10 sites can give. The
options are to add sites, accept a wider margin of error, or (if there
is evidence for it) assume a smaller ICC.

### More sites, fewer samples

The same precision can be reached with many different numbers of sites,
and spreading the samples over more sites reduces the total needed:

``` r

design_precision(prevalence = 0.2, moe = 0.05, n_sites = 15, icc = 0.05)$n
#> [1] 2070
design_precision(prevalence = 0.2, moe = 0.05, n_sites = 30, icc = 0.05)$n
#> [1] 450
design_precision(prevalence = 0.2, moe = 0.05, n_sites = 100, icc = 0.05)$n
#> [1] 300
```

With 15 sites the study needs 2070 samples; with 100 sites, only 300,
close to the 264 an unclustered survey would need. Each new site brings
in fresh information, whereas each extra sample at a site partly repeats
what that site has already shown. So, where the budget and logistics
allow, aim for more sites with fewer people at each. Site visits cost
money too, and the [Costing a
study](https://mrc-ide.github.io/MMSpower/articles/tutorial-cost.md)
tutorial shows how to find the cheapest design.

### The precision of a design you already have

Sometimes the design is set by circumstances: a network of sentinel
sites that can each collect a fixed number of samples, say. Give the
full design and leave out `moe`, and
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
returns the margin of error that design achieves:

``` r

d_have <- design_precision(prevalence = 0.2, n_sites = 20, n_per_site = 20, icc = 0.05)
d_have$moe
#> [1] 0.0569166
```

Twenty sites of 20 people give a margin of error of ±5.7 percentage
points, a little wider than the 5-point target.

Sites do not have to be the same size. Give `n_per_site` as one value
per site:

``` r

d_uneven <- design_precision(prevalence = 0.2, n_per_site = c(250, 250, 60, 60, 40, 40),
                             icc = 0.05)
d_uneven$moe
#> [1] 0.1025958

d_even <- design_precision(prevalence = 0.2, n_sites = 6, n_per_site = 117, icc = 0.05)
d_even$moe
#> [1] 0.08123823
```

These six sites hold 700 people and give a margin of error of ±10.3
points. Spreading about the same number evenly, 6 sites of 117, gives
±8.1. Large sites dominate the estimate, and many of their people repeat
what that site has already shown, so uneven site sizes cost precision.

### Allowing for dropout and suspected cases

The design functions give the number of samples needed **at analysis**.
More people than that have to be enrolled, for two reasons.

First, some samples are lost: consent is withdrawn, a sample is lost in
transit, or the assay fails. If about 10% are expected to be lost,
divide by $`1 - 0.1`$ and round up. Do this per site, not only for the
total, so that every site still delivers its share. Taking the design
with 20 samples at each of 26 sites:

``` r

per_site_analysis  <- d_per$n_per_site
per_site_confirmed <- ceiling(per_site_analysis / (1 - 0.1))
per_site_confirmed
#> [1] 23
```

Second, the marker is measured only in confirmed *P. falciparum*
infections, but people are enrolled as suspected cases, before the test
result is known. If about 40% of suspected cases are confirmed, divide
again by 0.4:

``` r

per_site_suspected <- ceiling(per_site_confirmed / 0.4)
per_site_suspected
#> [1] 58
```

So each site enrols 58 suspected cases, expecting 23 of them to be
confirmed, and 20 samples to reach the analysis. Across the 26 sites
that is 1508 suspected cases enrolled, 598 confirmed cases enrolled and
520 samples expected at analysis. Keep the three numbers apart when
planning and reporting: the design was calculated for the last one.

The confirmed fraction often differs between sites, for example with
transmission intensity. Use each site’s own value:

``` r

ceiling(per_site_confirmed / c(0.6, 0.4, 0.25))
#> [1] 39 58 92
```

A site where only a quarter of suspected cases are confirmed has to
enrol 92 suspected cases to deliver the same 20 samples. When budgeting,
cost the inflated numbers, since those are the people enrolled and
tested: in
[`estimate_cost()`](https://mrc-ide.github.io/MMSpower/reference/estimate_cost.md),
the suspected cases are `n_screened`, the confirmed cases `n_collected`
and the samples at analysis `n_sequenced` (see [Costing a
study](https://mrc-ide.github.io/MMSpower/articles/tutorial-cost.md)).

## Part 2: Analysing the data

Suppose the programme settles on 20 sites of about 30 people, and
records the number of samples carrying the marker (`x`) and the number
tested (`n`) at each site:

``` r

x <- c(4, 9, 6, 2, 11, 7, 5, 8, 3, 6, 10, 4, 7, 5, 12, 6, 3, 8, 5, 7)
n <- c(30, 32, 28, 31, 30, 29, 33, 30, 27, 31, 30, 32, 29, 30, 31, 28, 30, 30, 32, 29)
```

### The estimate

By default the ICC is estimated from the spread of the site prevalences,
and an exact (Clopper-Pearson) interval is returned:

``` r

est <- estimate_prevalence(x, n)
est
#> $prevalence
#> [1] 0.2126246
#> 
#> $ci_lower
#> [1] 0.172139
#> 
#> $ci_upper
#> [1] 0.2577502
#> 
#> $moe
#> [1] 0.04280562
#> 
#> $n_eff
#> [1] 424.3709
#> 
#> $icc
#> [1] 0.01434799
#> 
#> $deff
#> [1] 1.41857
#> 
#> $df
#> [1] 19
```

The estimated prevalence is 21.3%, with a 95% interval from 17.2% to
25.8%. `deff` is the design effect worked out from the estimated ICC,
and `n_eff` is the number of independent samples the clustered data are
worth. `df` is the number of sites minus one: the ICC was estimated from
20 sites, and the interval is widened a little to allow for the
uncertainty in that estimate (see below).

### Choosing an interval method

Three interval methods are available. They give the same point estimate
but different intervals:

``` r

cp <- estimate_prevalence(x, n, ci_method = "clopper-pearson")
c(cp$ci_lower, cp$ci_upper)
#> [1] 0.1721390 0.2577502

wald <- estimate_prevalence(x, n, ci_method = "wald")
c(wald$ci_lower, wald$ci_upper)
#> [1] 0.1710528 0.2541964

ac <- estimate_prevalence(x, n, ci_method = "agresti-coull")
c(ac$ci_lower, ac$ci_upper)
#> [1] 0.1739951 0.2571266
```

With this much data and a prevalence well away from 0 they are close.
The differences matter for small samples or a prevalence near 0 or 1.
There the Wald interval is too narrow, has zero width when no positives
are seen, and its limits can fall below 0 or above 1.
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)
sets such limits to 0 or 1, which keeps the interval sensible, but the
result is still less reliable than the other two methods in these
situations. Clopper-Pearson, the default, is the cautious choice, and
Agresti-Coull a good compromise. Use the same method at both stages:
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
takes the same `ci_method` argument, so a study sized for
Clopper-Pearson (the default in both) gives back the planned margin of
error when analysed with it. See
[`?estimate_prevalence`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)
for the details.

### Supplied or estimated ICC?

The interval depends on the ICC, and estimating an ICC needs data from
many sites: with only a handful, the estimate is unreliable.
[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)
allows for this by widening the interval according to the number of
sites the ICC was estimated from. With many sites this makes little
difference, but with a handful the interval becomes very wide. So with
few sites it is often better to supply a cautious value, for example
from earlier surveys in similar settings:

``` r

est_icc <- estimate_prevalence(x, n, icc = 0.05)
c(est_icc$ci_lower, est_icc$ci_upper)
#> [1] 0.1630936 0.2692575
```

With 20 sites the estimate is usable: it gave an ICC of 0.014 and an
interval from 17.2% to 25.8%, while supplying 0.05 gives 16.3% to 26.9%.
The [clustering
concepts](https://mrc-ide.github.io/MMSpower/articles/concept-clustering.md)
article shows how variable ICC estimates are with few sites.

### An imperfect test at analysis

If the assay is imperfect, give its sensitivity and specificity, and the
estimate and interval are corrected to the true-prevalence scale:

``` r

est_test <- estimate_prevalence(x, n, sensitivity = 0.9, specificity = 0.95)
c(est_test$prevalence, est_test$ci_lower, est_test$ci_upper)
#> [1] 0.1913230 0.1436929 0.2444121
```

The corrected prevalence is lower than the measured one, because at this
prevalence false positives outnumber missed cases, and the interval is
wider. The sensitivity and specificity given here should be the same
values assumed when the study was designed: the design only delivers the
planned margin of error if the analysis makes the same correction.
