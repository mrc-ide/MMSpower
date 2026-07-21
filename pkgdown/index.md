# MMSpower

**Status: early development.** This site is a working test build of the
`MMSpower` R package, currently used to trial the documentation and CI/CD
setup ahead of the package containing real statistical functionality.

## What this package will do

`MMSpower` is being developed to provide power calculations for MMS
(Malaria Molecular Surveillance) Study design. The scope of exactly which calculations are included is still
being finalised — see the [package overview](reference/MMSpower-package.html)
for the current agreed scope note.

## Current state

At this stage, the package contains a placeholder function
([`hello_mmspower()`](reference/hello_mmspower.html)) used to confirm that
the package builds, documents, tests, and deploys correctly end to end.
No statistical functions have been added yet.

## Getting started

```r
# install a specific released version
remotes::install_github("mrc-ide/MMSpower@v0.1.0")

# or install the latest development version
remotes::install_github("mrc-ide/MMSpower")
```

See the [Articles](articles/index.html) section for a walkthrough, and the
[Reference](reference/index.html) section for full function documentation.

## Development

This package is developed by Sequoia Leuba and Bob Verity, following a GitFlow
workflow (`main` / `develop` / feature branches) with continuous
integration on every pull request. See the repository's `develop` branch
for the latest in-progress work.
