# MMSpower

![MMSpower hex logo](logo.svg)

# MMSpower

Sample size, power and analysis for malaria molecular surveillance
studies, allowing for imperfect diagnostic tests and sampling clustered
across sites.

[Get started](https://mrc-ide.github.io/MMSpower/articles/MMSpower.md)
[Tutorials](https://mrc-ide.github.io/MMSpower/articles/index.md)
[Reference](https://mrc-ide.github.io/MMSpower/reference/index.md)

## Why MMSpower

Malaria programmes increasingly use molecular surveillance, for example
of drug-resistance markers or *pfhrp2/3* deletions, to make decisions.
These surveys are often clustered across sites, rely on imperfect
assays, and need to be tailored to local constraints. MMSpower provides
straightforward functions for designing a study and analysing its
results, with a focus on malaria molecular surveillance (MMS)
applications.

## A first example

How many samples are needed to estimate a prevalence of around 30% to
within ±5 percentage points, sampling 10 people per site with an ICC of
0.05?

``` r

library(MMSpower)

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

47 sites of 10 people. Without clustering, 323 people would be enough;
the design effect of 1.45 accounts for people at the same site
resembling each other.

The [tutorials](https://mrc-ide.github.io/MMSpower/articles/index.md)
take each type of study from design through to analysis, starting with
[estimating
prevalence](https://mrc-ide.github.io/MMSpower/articles/tutorial-prevalence.md).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mrc-ide/MMSpower")
```

## About

MMSpower is developed by Sequoia Leuba and Bob Verity at Imperial
College London. The source code is on
[GitHub](https://github.com/mrc-ide/MMSpower), where you can also
[report issues](https://github.com/mrc-ide/MMSpower/issues).
