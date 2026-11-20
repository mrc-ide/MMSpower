# Costing a study

``` r

library(MMSpower)
```

A design that meets its statistical target is only useful if it can be
paid for.
[`estimate_cost()`](https://mrc-ide.github.io/MMSpower/reference/estimate_cost.md)
puts a price on a design, and because it takes its inputs straight from
the design functions, it can be used to compare designs that are
statistically equivalent but cost very different amounts.

All costs in this article are **illustrative**, in an unnamed currency
unit. Replace them with figures from your own budget.

## The cost model

[`estimate_cost()`](https://mrc-ide.github.io/MMSpower/reference/estimate_cost.md)
splits a study’s cost into two parts.

**Stage costs.** A sample passes through up to three stages, each with
its own cost per person or sample:

- **screened** (`n_screened`, `cost_per_screened`): everyone screened at
  enrolment, covering consent, interview and a rapid test or microscopy;
- **collected** (`n_collected`, `cost_per_collected`): everyone a sample
  is collected from, covering a dried blood spot or blood draw, storage
  and shipping;
- **sequenced** (`n_sequenced`, `cost_per_sequenced`): every sample
  taken through the molecular laboratory, covering extraction, PCR or
  sequencing and data management.

The stages are nested, not added on top of each other. `n_collected` is
a subset of `n_screened`: for example, only those who test positive at
screening have a sample collected. `n_sequenced` is a subset of
`n_collected`: what is left after lost or failed samples. A person is
counted at every stage they reach, so `n_screened` is at least
`n_collected`, which is at least `n_sequenced`, and the function stops
with an error otherwise. Each stage’s cost is the cost of that stage
alone, not of everything up to it.

**Site costs.** Each site has a fixed cost (training, equipment,
administration) and a transport cost.

Give only the stages a study has. The design functions return the number
of samples needed at analysis, so their `n` is `n_sequenced`. With a
single all-in cost per sample, that one stage is enough:

``` r

estimate_cost(n_sequenced = 400, cost_per_sequenced = 40)
#> $total_cost
#> [1] 16000
#> 
#> $variable_cost
#> [1] 16000
#> 
#> $stage_cost
#>  screened collected sequenced 
#>         0         0     16000 
#> 
#> $site_cost
#> [1] 0
#> 
#> $budget_remaining
#> NULL
```

`stage_cost` gives the cost of each stage, and adds up to
`variable_cost`. With sites, give the number of sites and the cost per
site:

``` r

estimate_cost(n_sequenced = 400, cost_per_sequenced = 40, n_sites = 20,
              fixed_cost_per_site = 2000, transport_cost_per_site = 600)
#> $total_cost
#> [1] 68000
#> 
#> $variable_cost
#> [1] 16000
#> 
#> $stage_cost
#>  screened collected sequenced 
#>         0         0     16000 
#> 
#> $site_cost
#> [1] 52000
#> 
#> $budget_remaining
#> NULL
```

Site costs need not be the same everywhere. Remote sites usually cost
more to reach, so either per-site cost can be given as one value per
site. Adding a `budget` reports how much would be left over (negative if
the study is over budget):

``` r

transport <- c(200, 200, 300, 400, 600, 800, 1200, 1500)
estimate_cost(n_sequenced = 240, cost_per_sequenced = 40, n_sites = 8,
              fixed_cost_per_site = 2000, transport_cost_per_site = transport,
              budget = 35000)
#> $total_cost
#> [1] 30800
#> 
#> $variable_cost
#> [1] 9600
#> 
#> $stage_cost
#>  screened collected sequenced 
#>         0         0      9600 
#> 
#> $site_cost
#> [1] 21200
#> 
#> $budget_remaining
#> [1] 4200
```

## Costing the people actually enrolled

A design’s `n` is only the last stage. To reach it, a study has to
collect samples from more people, to allow for losses, and screen more
people still, if only those who test positive go on to have a sample
collected. “Allowing for dropout and suspected cases” in the [Estimating
prevalence](https://mrc-ide.github.io/MMSpower/articles/tutorial-prevalence.md)
tutorial works out these numbers; here we cost them.

That tutorial’s design has 20 people per site, with an ICC of 0.05:

``` r

d_prev <- design_precision(prevalence = 0.2, moe = 0.05, n_per_site = 20, icc = 0.05)
c(n_sites = d_prev$n_sites, n = d_prev$n)
#> n_sites       n 
#>      26     520
```

Work back from the analysis one stage at a time, per site so that every
site delivers its share. With 10% of samples expected to be lost, divide
by 0.9 to get the number collected; with 40% of those screened expected
to test positive, divide again by 0.4 to get the number screened:

``` r

per_site_sequenced <- d_prev$n_per_site
per_site_collected <- ceiling(per_site_sequenced / 0.9)
per_site_screened  <- ceiling(per_site_collected / 0.4)

n_sequenced <- d_prev$n
n_collected <- d_prev$n_sites * per_site_collected
n_screened  <- d_prev$n_sites * per_site_screened
c(screened = n_screened, collected = n_collected, sequenced = n_sequenced)
#>  screened collected sequenced 
#>      1508       598       520
```

Then cost all three stages, here 3 per person screened, 10 per sample
collected and 40 per sample sequenced:

``` r

staged <- estimate_cost(n_screened  = n_screened,  cost_per_screened  = 3,
                        n_collected = n_collected, cost_per_collected = 10,
                        n_sequenced = n_sequenced, cost_per_sequenced = 40,
                        n_sites = d_prev$n_sites, fixed_cost_per_site = 2000,
                        transport_cost_per_site = 600)
staged$stage_cost
#>  screened collected sequenced 
#>      4524      5980     20800
staged$total_cost
#> [1] 98904
```

Screening 1,508 people and collecting 598 samples add 10,504 to the
20,800 spent on sequencing. Costing only the 520 samples at analysis
would have missed both stages.

If samples are expected to fail in the laboratory rather than be lost
before it, the failed ones are still paid for: cost `n_sequenced` as the
number of samples sent to the laboratory, here `n_collected`, rather
than the number expected to succeed.

## Comparing designs

The main use of
[`estimate_cost()`](https://mrc-ide.github.io/MMSpower/reference/estimate_cost.md)
is choosing between designs. Suppose we want to estimate a prevalence of
around 30% to within ±5 percentage points, with an ICC of 0.05.
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
can meet that target with many combinations of sites and people per
site. Few people per site needs many sites; many people per site needs
fewer sites, but a larger design effect means more people in total.

Work out the design for each number of people per site from 5 to 60, and
cost each one:

``` r

cost_per_sample <- 40     # all-in cost per sample sequenced
site_cost       <- 2000   # fixed cost per site
transport_cost  <- 600    # transport cost per site

per_site <- 5:60
designs <- do.call(rbind, lapply(per_site, function(s) {
  des  <- design_precision(prevalence = 0.3, moe = 0.05, n_per_site = s, icc = 0.05)
  cost <- estimate_cost(n_sequenced = des$n, cost_per_sequenced = cost_per_sample,
                        n_sites = des$n_sites,
                        fixed_cost_per_site = site_cost,
                        transport_cost_per_site = transport_cost)
  data.frame(n_per_site = s, n_sites = des$n_sites, n = des$n,
             deff = des$deff, total_cost = cost$total_cost)
}))
head(designs)
#>   n_per_site n_sites   n deff total_cost
#> 1          5      82 410 1.20     229600
#> 2          6      71 426 1.25     201640
#> 3          7      64 448 1.30     184320
#> 4          8      58 464 1.35     169360
#> 5          9      53 477 1.40     156880
#> 6         10      50 500 1.45     150000
```

``` r

best <- designs[which.min(designs$total_cost), ]
best
#>    n_per_site n_sites   n deff total_cost
#> 32         36      26 936 2.75     105040
```

``` r

library(ggplot2)

ggplot(designs, aes(x = n_per_site, y = total_cost / 1000)) +
  geom_line(colour = "#1B2A41", linewidth = 0.8) +
  geom_point(colour = "#1B2A41", size = 1.2) +
  geom_point(data = best, colour = "#C0392B", size = 4) +
  geom_text(data = best, aes(label = sprintf("%d sites of %d", n_sites, n_per_site)),
            colour = "#6B7280", vjust = -1.4) +
  labs(x = "People per site", y = "Total cost (thousands)") +
  theme_classic(base_size = 13)
```

![Total study cost against people per site, falling steeply then
levelling off, with the cheapest design
marked.](tutorial-cost_files/figure-html/unnamed-chunk-10-1.png)

Every point on this curve gives the same precision. On the left,
sampling 5 people per site needs 82 sites, and the site costs dominate.
Moving right, each extra person per site saves sites, and the cost falls
quickly. Past a point the saving stops: the design effect keeps growing,
so the extra people needed to make up for it cost more than the sites
they save. The cheapest design here is 26 sites of 36 people (936 in
all, design effect 2.75), costing 105,040.

The curve is jagged because the number of sites is a whole number: a
design that needs 23.1 sites pays for 24. Near the bottom the curve is
also quite flat, so designs a little either side of the cheapest cost
much the same. That leaves room to choose on practical grounds, such as
how many people a field team can sample at a site in a day.

The shape depends on the ratio of site cost to sample cost. When sites
are expensive relative to samples, the cheapest design has more people
per site; when samples are expensive (for example, with whole-genome
sequencing), it has fewer. A higher ICC also pushes towards fewer people
per site, since each additional person at a site is worth less.

## Costing a comparison of two groups

The same works for any design function. For
[`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md),
note that `n` is already the total across both groups, but `n_sites` is
given per group, so the number of sites to cost is `sum(n_sites)`:

``` r

dd <- design_difference(prevalence1 = 0.10, prevalence2 = 0.20,
                        n_per_site1 = 20, icc = 0.05)
dd$n
#> [1] 800
dd$n_sites
#> [1] 20 20
estimate_cost(n_sequenced = dd$n, cost_per_sequenced = cost_per_sample,
              n_sites = sum(dd$n_sites),
              fixed_cost_per_site = site_cost,
              transport_cost_per_site = transport_cost)
#> $total_cost
#> [1] 136000
#> 
#> $variable_cost
#> [1] 32000
#> 
#> $stage_cost
#>  screened collected sequenced 
#>         0         0     32000 
#> 
#> $site_cost
#> [1] 104000
#> 
#> $budget_remaining
#> NULL
```

If the two regions have different costs (one is much further away, say),
call
[`estimate_cost()`](https://mrc-ide.github.io/MMSpower/reference/estimate_cost.md)
once per region and add the totals.

## Fitting a budget

To find which designs fit a fixed budget, keep the rows of the table
above that come in under it. With a budget of 110,000:

``` r

budget <- 110000
fits <- designs[designs$total_cost <= budget, ]
range(fits$n_per_site)
#> [1] 25 54
```

25 of the 56 designs meet the precision target within budget, with
between 25 and 54 people per site. Because the curve is jagged, not
every value in that range fits; check the table rather than assuming it
does. If no row fits, the target itself has to change: a wider margin of
error, or a cheaper laboratory method. For the chosen design, `budget`
shows what is left for contingencies:

``` r

estimate_cost(n_sequenced = best$n, cost_per_sequenced = cost_per_sample,
              n_sites = best$n_sites,
              fixed_cost_per_site = site_cost,
              transport_cost_per_site = transport_cost,
              budget = budget)$budget_remaining
#> [1] 4960
```
