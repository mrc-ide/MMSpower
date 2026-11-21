# MMSpower

<div class="mms-hero">
<img src="logo.svg" alt="MMSpower hex logo">
<div>
<h1>MMSpower</h1>
<p class="lead">Sample size, power and analysis for malaria molecular surveillance studies, allowing for imperfect diagnostic tests and sampling clustered across sites.</p>
<a class="btn btn-primary btn-lg" href="articles/MMSpower.html">Get started</a>
<a class="btn btn-outline-primary btn-lg" href="articles/index.html">Tutorials</a>
<a class="btn btn-outline-secondary btn-lg" href="reference/index.html">Reference</a>
</div>
</div>

## Why MMSpower

Malaria programmes increasingly use molecular surveillance, for example of
drug-resistance markers or *pfhrp2/3* deletions, to make decisions. These
surveys are often clustered across sites, rely on imperfect assays, and need
to be tailored to local constraints. MMSpower provides straightforward
functions for designing a study and analysing its results, with a focus on
malaria molecular surveillance (MMS) applications.

## A first example

How many samples are needed to estimate a prevalence of around 30% to within
±5 percentage points, sampling 10 people per site with an ICC of 0.05?

```r
library(MMSpower)

design_precision(prevalence = 0.3, moe = 0.05, n_per_site = 10, icc = 0.05)
#> $n
#> [1] 500
#>
#> $n_sites
#> [1] 50
#>
#> $n_per_site
#> [1] 10
#>
#> $deff
#> [1] 1.45
#>
#> $n_eff
#> [1] 341
#>
#> $moe
#> [1] 0.04965355
```

50 sites of 10 people, giving a margin of error just under the 5-point
target. Without clustering, 341 people would be enough; the design effect of
1.45 accounts for people at the same site resembling each other.

The [tutorials](articles/index.html) take each type of study from design
through to analysis, starting with
[estimating prevalence](articles/tutorial-prevalence.html).

## Installation

Install the latest release from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("mrc-ide/MMSpower")
```

or the development version with
`remotes::install_github("mrc-ide/MMSpower@develop")`.

## About

MMSpower is developed by Sequoia Leuba and Bob Verity at Imperial College
London. The source code is on
[GitHub](https://github.com/mrc-ide/MMSpower), where you can also
[report issues](https://github.com/mrc-ide/MMSpower/issues).
