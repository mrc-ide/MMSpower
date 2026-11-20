# Imperfect diagnostic tests

``` r

library(MMSpower)
```

Every MMSpower function takes a `sensitivity` and a `specificity`. This
article explains what they do: how an imperfect test changes the
prevalence you measure, how MMSpower corrects for it, and why it means
you need more samples.

## Sensitivity and specificity

A diagnostic test, whether a rapid diagnostic test, microscopy or a PCR
assay for a resistance marker, can go wrong in two ways.

- **Sensitivity** ($`Se`$) is the chance that the test is positive in
  someone who truly is positive. A test with $`Se = 0.9`$ misses 1 true
  positive in 10.
- **Specificity** ($`Sp`$) is the chance that the test is negative in
  someone who truly is negative. A test with $`Sp = 0.95`$ gives a false
  positive in 1 true negative in 20.

The default in every function is $`Se = Sp = 1`$, a perfect test
([`design_detection()`](https://mrc-ide.github.io/MMSpower/reference/design_detection.md)
takes only the sensitivity, as positives are assumed to be confirmed
before a variant is declared present). Values for a particular test
usually come from a validation study against a reference standard.

## Apparent prevalence

The proportion of positive results you actually see is the **apparent**
(or measured) prevalence. With true prevalence $`p`$, positives come
from true positives the test detects plus true negatives it gets wrong:

``` math
p_{app} = p \cdot Se + (1 - p)(1 - Sp).
```

The plot shows $`p_{app}`$ against $`p`$ for a few tests. The grey
diagonal is a perfect test.

``` r

p <- seq(0, 1, length.out = 201)
tests <- data.frame(se  = c(0.95, 0.90, 0.80),
                    sp  = c(0.98, 0.95, 0.90),
                    col = c("#1B2A41", "#C0392B", "#8A96A8"))

par(mar = c(4.5, 4.5, 1, 1), las = 1, bty = "l")
plot(p, p, type = "l", col = "#6B7280", lty = 2, lwd = 2,
     xlab = "True prevalence", ylab = "Apparent (measured) prevalence")
for (i in seq_len(nrow(tests))) {
  lines(p, p * tests$se[i] + (1 - p) * (1 - tests$sp[i]),
        col = tests$col[i], lwd = 2)
}
legend("topleft", bty = "n", lwd = 2, lty = c(2, 1, 1, 1),
       col = c("#6B7280", tests$col),
       legend = c("Perfect test",
                  sprintf("Se = %.2f, Sp = %.2f", tests$se, tests$sp)))
```

![](concept-diagnostic-tests_files/figure-html/apparent-plot-1.png)

Two things stand out. First, each line starts at $`1 - Sp`$, not at 0:
even where nobody is infected, false positives give a measured
prevalence of $`1 - Sp`$. Second, each line is flatter than the
diagonal. Its slope is $`Se + Sp - 1`$, so a given change in true
prevalence produces a smaller change in what is measured. At low
prevalence false positives dominate and the test overstates prevalence;
at high prevalence missed cases dominate and it understates it.

## The Rogan-Gladen correction

Because $`p_{app}`$ is a straight line in $`p`$, it can be turned round
to recover true prevalence from what was measured (Rogan and Gladen,
1978):

``` math
p = \frac{p_{app} - (1 - Sp)}{Se + Sp - 1}.
```

This only works if the denominator is positive, which is why every
function requires $`Se + Sp > 1`$. When $`Se + Sp = 1`$ the line in the
plot is flat: the test is positive with the same chance whether or not
the person is infected, so the result says nothing about prevalence.
When $`Se + Sp < 1`$ the test does worse than tossing a coin, and
swapping “positive” and “negative” would give a better one.

``` r

design_precision(prevalence = 0.3, moe = 0.05, sensitivity = 0.5, specificity = 0.5)
#> Error:
#> ! `sensitivity` + `specificity` must be greater than 1.
```

## Why an imperfect test needs more samples

The slope $`Se + Sp - 1`$ is the key quantity. The difference between
two true prevalences $`p_1`$ and $`p_2`$ appears in the measurements as

``` math
p_{1,app} - p_{2,app} = (p_1 - p_2)(Se + Sp - 1).
```

An imperfect test shrinks every difference you are trying to detect, and
also shrinks the gap between a prevalence and a threshold. To see a
smaller difference you need more samples. For estimation the same idea
applies in reverse: an error of $`\pm e`$ on the measured scale becomes
$`\pm e / (Se + Sp - 1)`$ after correction. In the Wald formula for the
sample size this appears as a division by $`(Se + Sp - 1)^2`$:

``` math
n_0 = \frac{z^2 \, p_{app}(1-p_{app})}{MOE^2 \cdot (Se + Sp - 1)^2}.
```

[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
uses this formula with `ci_method = "wald"`. Its default,
Clopper-Pearson, finds the sample size by a numerical search, but the
effect of an imperfect test is the same.

## Specificity matters most at low prevalence

The table below holds sensitivity at 0.95 and lowers specificity from 1
to 0.9. It uses a margin of error of 2 percentage points at a prevalence
of 5%, and 5 percentage points at a prevalence of 30%, and gives the
sample size from
[`design_precision()`](https://mrc-ide.github.io/MMSpower/reference/design_precision.md)
together with its ratio to the sample size at perfect specificity.

``` r

spec <- c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90)
n_low  <- sapply(spec, function(sp)
  design_precision(0.05, moe = 0.02, sensitivity = 0.95, specificity = sp)$n)
n_high <- sapply(spec, function(sp)
  design_precision(0.30, moe = 0.05, sensitivity = 0.95, specificity = sp)$n)

data.frame(specificity = spec,
           n_prev_5  = n_low,  ratio_5  = round(n_low / n_low[1], 2),
           n_prev_30 = n_high, ratio_30 = round(n_high / n_high[1], 2))
#>    specificity n_prev_5 ratio_5 n_prev_30 ratio_30
#> 1         1.00      536    1.00       366     1.00
#> 2         0.99      639    1.19       379     1.04
#> 3         0.98      744    1.39       392     1.07
#> 4         0.97      851    1.59       405     1.11
#> 5         0.96      962    1.79       419     1.14
#> 6         0.95     1075    2.01       433     1.18
#> 7         0.94     1190    2.22       447     1.22
#> 8         0.93     1309    2.44       462     1.26
#> 9         0.92     1430    2.67       477     1.30
#> 10        0.91     1554    2.90       493     1.35
#> 11        0.90     1682    3.14       508     1.39
```

``` r

par(mar = c(4.5, 4.5, 1, 1), las = 1, bty = "l")
plot(spec, n_low / n_low[1], type = "b", pch = 19, col = "#1B2A41", lwd = 2,
     xlim = c(1, 0.9), ylim = c(1, max(n_low / n_low[1])),
     xlab = "Specificity (sensitivity fixed at 0.95)",
     ylab = "Sample size relative to Sp = 1")
lines(spec, n_high / n_high[1], type = "b", pch = 19, col = "#C0392B", lwd = 2)
legend("topleft", bty = "n", pch = 19, lwd = 2, col = c("#1B2A41", "#C0392B"),
       legend = c("Prevalence 5%", "Prevalence 30%"))
```

![](concept-diagnostic-tests_files/figure-html/spec-plot-1.png)

Dropping specificity to 0.9 multiplies the sample size by 3.1 at 5%
prevalence but only by 1.4 at 30%. At 5% prevalence a specificity of 0.9
means false positives (about 10% of the 95% who are uninfected)
outnumber the true positives, so most of what is measured is noise.
Specificity matters most when the condition is rare, which is often the
case in molecular surveillance, for example when looking for *pfhrp2/3*
deletions or a newly emerging resistance marker near a 5% threshold.

## Correcting an estimate

[`estimate_prevalence()`](https://mrc-ide.github.io/MMSpower/reference/estimate_prevalence.md)
builds its confidence interval on the apparent prevalence first, then
applies the Rogan-Gladen formula to the point estimate and to both ends
of the interval. Because the correction is a straight line, this keeps
the interval’s coverage. With 30 positives out of 100 and a test with
sensitivity 0.9 and specificity 0.95:

``` r

est <- estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95)
est
#> $prevalence
#> [1] 0.2941176
#> 
#> $ci_lower
#> [1] 0.1910664
#> 
#> $ci_upper
#> [1] 0.4115467
#> 
#> $moe
#> [1] 0.1102402
#> 
#> $n_eff
#> [1] 100
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

The apparent prevalence is 0.30. The corrected prevalence is
$`(0.30 - 0.05) / 0.85 =`$ 0.294, and the interval is stretched by
$`1 / 0.85`$, so the margin of error is 0.11 rather than the 0.094 of a
perfect test.

Corrected values below 0 or above 1 are set to 0 or 1. If the **whole**
corrected interval lies outside $`[0, 1]`$, the function warns. This
happens when there are fewer positives than false positives alone would
produce. With a specificity of 0.95, about 10 of 200 truly negative
people would test positive, so 1 positive out of 200 does not fit:

``` r

estimate_prevalence(x = 1, n = 200, specificity = 0.95)
#> Warning: The corrected interval lies entirely outside [0, 1], so the data do
#> not fit the assumed sensitivity and specificity.
#> $prevalence
#> [1] 0
#> 
#> $ci_lower
#> [1] 0
#> 
#> $ci_upper
#> [1] 0
#> 
#> $moe
#> [1] 0
#> 
#> $n_eff
#> [1] 200
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

The returned interval is then $`[0, 0]`$, which should not be reported
as a precise estimate of zero. The warning means that the assumed
specificity is too low for these data, or that the test performed better
in this study than in its validation. Check where the sensitivity and
specificity values came from before reporting a result.

## Summary

- An imperfect test measures $`p_{app} = p \cdot Se + (1-p)(1-Sp)`$, not
  $`p`$.
- The Rogan-Gladen formula recovers $`p`$, and needs $`Se + Sp > 1`$.
- Differences and errors on the true scale are shrunk by $`Se + Sp - 1`$
  on the measured scale, so sample sizes grow by roughly
  $`1 / (Se + Sp - 1)^2`$ and more when false positives inflate
  $`p_{app}(1 - p_{app})`$.
- At low prevalence, small losses in specificity are costly.

Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
screening test. *American Journal of Epidemiology* 107(1):71-76.
