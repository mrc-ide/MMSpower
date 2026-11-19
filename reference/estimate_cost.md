# Estimate total study cost

Costs an MMS study from its design choices, as a variable cost plus a
per-site cost:

- **Variable cost** = total samples enrolled times the cost per sample
  (collection, laboratory testing, consumables, data management).

- **Site cost** = the sum over sites of (fixed cost per site + transport
  cost per site). The fixed cost (training, equipment, administration)
  and the transport cost may each be given once for all sites or as one
  value per site.

Pairs with the design functions: take `n` and `n_sites` from a design
function's result. To explore trade-offs, change `n` or `n_sites` and
call it again. To cost several regions, call `estimate_cost()` once per
region and add the totals. All monetary inputs and outputs are in one
currency unit.

## Usage

``` r
estimate_cost(
  n,
  cost_per_sample,
  n_sites = NULL,
  fixed_cost_per_site = 0,
  transport_cost_per_site = 0,
  budget = NULL
)
```

## Arguments

- n:

  Whole number. Total number of samples to be enrolled: `$n` from any
  design function (for
  [`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md),
  the total across both groups).

- cost_per_sample:

  Non-negative number. Variable cost per enrolled sample.

- n_sites:

  Optional whole number. Total number of sites to cost: `$n_sites` from
  a design function (for
  [`design_difference()`](https://mrc-ide.github.io/MMSpower/reference/design_difference.md),
  `sum(n_sites)`). `NULL` (default) costs the samples only.

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

  `n * cost_per_sample`.

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
# Variable cost only
estimate_cost(n = 323, cost_per_sample = 50)
#> $total_cost
#> [1] 16150
#> 
#> $variable_cost
#> [1] 16150
#> 
#> $site_cost
#> [1] 0
#> 
#> $budget_remaining
#> NULL
#> 

# Same costs at every site (5,000 fixed + 1,000 transport per site)
estimate_cost(n = 600, cost_per_sample = 50, n_sites = 12,
              fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
#> $total_cost
#> [1] 102000
#> 
#> $variable_cost
#> [1] 30000
#> 
#> $site_cost
#> [1] 72000
#> 
#> $budget_remaining
#> NULL
#> 

# A different fixed cost at each of 3 sites, checked against a budget
estimate_cost(n = 50, cost_per_sample = 10, n_sites = 3,
              fixed_cost_per_site = c(100, 200, 300), budget = 1500)
#> $total_cost
#> [1] 1100
#> 
#> $variable_cost
#> [1] 500
#> 
#> $site_cost
#> [1] 600
#> 
#> $budget_remaining
#> [1] 400
#> 

# Cost a design
des <- design_precision(prevalence = 0.15, moe = 0.05, n_per_site = 30, icc = 0.05)
estimate_cost(n = des$n, cost_per_sample = 50, n_sites = des$n_sites,
              fixed_cost_per_site = 5000, transport_cost_per_site = 1200)
#> $total_cost
#> [1] 123200
#> 
#> $variable_cost
#> [1] 24000
#> 
#> $site_cost
#> [1] 99200
#> 
#> $budget_remaining
#> NULL
#> 

# Cost a two-group design: design_difference() gives one n_sites per
# group, so add them up
des2 <- design_difference(0.10, 0.20, n_per_site = 30, icc = 0.05)
estimate_cost(n = des2$n, cost_per_sample = 50, n_sites = sum(des2$n_sites),
              fixed_cost_per_site = 5000)
#> $total_cost
#> [1] 221000
#> 
#> $variable_cost
#> [1] 51000
#> 
#> $site_cost
#> [1] 170000
#> 
#> $budget_remaining
#> NULL
#> 
```
