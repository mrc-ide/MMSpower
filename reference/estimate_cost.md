# Estimate total study cost

Costs an MMS study from its design choices, as a cost for each stage a
sample passes through plus a per-site cost:

- **Stage costs**: for each of up to three stages, the number of people
  or samples reaching that stage times that stage's cost per person or
  sample. The stages are:

  - screened:

    people screened for infection at enrolment: consent, interview, a
    rapid test or microscopy.

  - collected:

    people a sample is collected from: a dried blood spot or blood draw,
    storage and shipping.

  - sequenced:

    samples taken through the molecular laboratory: extraction, PCR or
    sequencing, data management.

- **Site cost** = the sum over sites of (fixed cost per site + transport
  cost per site). The fixed cost (training, equipment, administration)
  and the transport cost may each be given once for all sites or as one
  value per site.

The stages are nested: everyone collected was first screened, and every
sample sequenced was first collected. So `n_collected` is a subset of
`n_screened`, not in addition to it (for example, only those who test
positive), and `n_sequenced` is a subset of `n_collected` (for example,
after lost or failed samples). A person is counted at every stage they
reach, and each stage's cost is the cost of that stage alone, not of
everything up to it. Give only the stages a study has; at least one is
needed. A study costed with a single all-in cost per sample can give
that as `n_sequenced` and `cost_per_sequenced`.

Pairs with the design functions: their `$n` is the number of samples
needed at analysis, so it is `n_sequenced`; scale it up for losses and
for suspected cases that turn out negative to get `n_collected` and
`n_screened`. Take `n_sites` from the design function's result too. To
cost several regions, call `estimate_cost()` once per region and add the
totals. All monetary inputs and outputs are in one currency unit.

## Usage

``` r
estimate_cost(
  n_screened = NULL,
  cost_per_screened = NULL,
  n_collected = NULL,
  cost_per_collected = NULL,
  n_sequenced = NULL,
  cost_per_sequenced = NULL,
  n_sites = NULL,
  fixed_cost_per_site = 0,
  transport_cost_per_site = 0,
  budget = NULL
)
```

## Arguments

- n_screened, n_collected, n_sequenced:

  Optional whole numbers. The number of people or samples reaching each
  stage (see Description). Each is no larger than the one before it. For
  `n_sequenced`, `$n` from any design function (for
  [`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md),
  the total across both groups).

- cost_per_screened, cost_per_collected, cost_per_sequenced:

  Non-negative numbers. The cost of each stage per person or sample
  reaching it. Give each one with its `n_` argument, or neither.

- n_sites:

  Optional whole number. Total number of sites to cost: `$n_sites` from
  a design function (for
  [`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md),
  `sum(n_sites)`). `NULL` (default) costs the stages only.

- fixed_cost_per_site:

  Non-negative fixed cost per site: one number for every site, or one
  value per site (length `n_sites`). Default 0.

- transport_cost_per_site:

  Non-negative transport cost per site: one number for every site, or
  one value per site. Default 0.

- budget:

  Optional positive number. Total money available.

## Value

A list with:

- total_cost:

  `variable_cost + site_cost`.

- variable_cost:

  The sum of `stage_cost`.

- stage_cost:

  Named vector of the cost of each stage (`screened`, `collected`,
  `sequenced`): its number times its cost, or 0 for a stage not given.

- site_cost:

  Sum over sites of fixed plus transport cost; 0 when `n_sites` is
  `NULL`.

- budget_remaining:

  `budget - total_cost` (negative if over budget), or `NULL` when no
  `budget` is given.

## References

MMS-SD Study Design Workshop, budget-officer activity; Modules 2 and 7.
<https://mrc-ide.github.io/MMS-SD_workshop/>

## Examples

``` r
# One all-in cost per sample
estimate_cost(n_sequenced = 323, cost_per_sequenced = 50)
#> $total_cost
#> [1] 16150
#> 
#> $variable_cost
#> [1] 16150
#> 
#> $stage_cost
#>  screened collected sequenced 
#>         0         0     16150 
#> 
#> $site_cost
#> [1] 0
#> 
#> $budget_remaining
#> NULL
#> 

# Three stages: 1,950 screened, the 780 positives sampled, 700 sequenced
# after losses; 12 sites at 5,000 fixed + 1,000 transport each
estimate_cost(n_screened  = 1950, cost_per_screened  = 2,
              n_collected = 780,  cost_per_collected = 15,
              n_sequenced = 700,  cost_per_sequenced = 40,
              n_sites = 12, fixed_cost_per_site = 5000,
              transport_cost_per_site = 1000)
#> $total_cost
#> [1] 115600
#> 
#> $variable_cost
#> [1] 43600
#> 
#> $stage_cost
#>  screened collected sequenced 
#>      3900     11700     28000 
#> 
#> $site_cost
#> [1] 72000
#> 
#> $budget_remaining
#> NULL
#> 

# A different fixed cost at each of 3 sites, checked against a budget
estimate_cost(n_sequenced = 50, cost_per_sequenced = 10, n_sites = 3,
              fixed_cost_per_site = c(100, 200, 300), budget = 1500)
#> $total_cost
#> [1] 1100
#> 
#> $variable_cost
#> [1] 500
#> 
#> $stage_cost
#>  screened collected sequenced 
#>         0         0       500 
#> 
#> $site_cost
#> [1] 600
#> 
#> $budget_remaining
#> [1] 400
#> 

# Cost a design, allowing 10% of samples to be lost and 40% of suspected
# cases to be positive
des <- design_precision(prevalence = 0.15, moe = 0.05, n_per_site = 30, icc = 0.05)
n_collected <- ceiling(des$n / 0.9)
estimate_cost(n_screened  = ceiling(n_collected / 0.4), cost_per_screened  = 2,
              n_collected = n_collected,                cost_per_collected = 15,
              n_sequenced = des$n,                      cost_per_sequenced = 40,
              n_sites = des$n_sites, fixed_cost_per_site = 5000,
              transport_cost_per_site = 1200)
#> $total_cost
#> [1] 145200
#> 
#> $variable_cost
#> [1] 33600
#> 
#> $stage_cost
#>  screened collected sequenced 
#>      3000      9000     21600 
#> 
#> $site_cost
#> [1] 111600
#> 
#> $budget_remaining
#> NULL
#> 

# Cost a two-group design: design_difference() gives one n_sites per
# group, so add them up
des2 <- design_difference(0.10, 0.20, n_per_site1 = 30, icc = 0.05)
estimate_cost(n_sequenced = des2$n, cost_per_sequenced = 50,
              n_sites = sum(des2$n_sites), fixed_cost_per_site = 5000)
#> $total_cost
#> [1] 221000
#> 
#> $variable_cost
#> [1] 51000
#> 
#> $stage_cost
#>  screened collected sequenced 
#>         0         0     51000 
#> 
#> $site_cost
#> [1] 170000
#> 
#> $budget_remaining
#> NULL
#> 
```
