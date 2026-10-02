#' Estimate prevalence and its precision
#'
#' @description
#' Analysis function: given observed counts, returns a prevalence point
#' estimate and confidence interval. Three CI methods are available:
#'
#' \describe{
#'   \item{`"wald"` (default)}{Symmetric unless an end is moved to 0 or 1.
#'     \eqn{\hat{p} \pm z \cdot SE}. Can be too narrow for small samples or
#'     prevalence near 0 or 1.}
#'   \item{`"clopper-pearson"`}{Asymmetric. Exact binomial interval, via
#'     the beta distribution. For a simple random sample, its coverage is
#'     at least the stated confidence level, so it is usually wider than
#'     the Wald interval.}
#'   \item{`"agresti-coull"`}{Asymmetric. A Wald interval calculated after
#'     adding extra positive and negative cases to the data, the same number
#'     of each. The number depends on the confidence level: about 1.35 of
#'     each for 90 percent, 1.92 for 95 percent, and 3.32 for 99 percent.
#'     This moves the interval towards 0.5; the reported prevalence is
#'     unchanged.}
#' }
#'
#' All methods account for imperfect diagnostic tests (Rogan-Gladen
#' correction), clustered sampling (Kish design effect), and
#' finite-population corrections.
#'
#' If `icc` is not supplied, it is estimated from the data (the observed
#' variance between clusters compared with the variance expected if every
#' person were independent). Setting `icc = 0` treats every person as
#' independent; if people at the same site tend to be similar, this gives an
#' interval that is too narrow.
#'
#' @param x Vector of integer counts of people who tested positive, one per
#'   cluster/site.
#' @param n Vector of integer totals tested, one per cluster/site (same
#'   length as `x`).
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case). Set below 1 to correct for
#'   missed cases (Rogan-Gladen correction).
#' @param specificity Diagnostic specificity in (0, 1]: the share of people
#'   who truly do not have the condition that the test correctly reads as
#'   negative. Default 1 (the test never gives a false positive). Set below
#'   1 to correct for false positives (Rogan-Gladen correction).
#' @param conf_level Confidence level, in (0, 1); default 0.95.
#' @param icc Optional. Intra-cluster correlation, in \[0, 1\]. If `NULL`
#'   (default), ICC is estimated from the data. Set to `0` for no
#'   clustering. Only relevant when `x` and `n` have more than one
#'   element. Values below about 1.5e-8 are treated as 0.
#' @param fpc_N Optional positive integer. Total population size for a
#'   finite-population correction; must be larger than `sum(n)`. `NULL`
#'   (default) = no finite-population correction applied.
#' @param method CI method: `"wald"` (default), `"clopper-pearson"`, or
#'   `"agresti-coull"`. See Description. Clopper-Pearson and Agresti-Coull
#'   produce asymmetric intervals. A Wald interval becomes asymmetric when an
#'   end is moved to 0 or 1, and the function prints a message when that
#'   happens (see Details).
#'
#' @details
#' The interval is built in six steps: apparent prevalence, a design
#' effect for clustering, an effective sample size, an optional
#' finite-population correction, a confidence interval by the chosen
#' method, and finally the Rogan-Gladen correction for an imperfect test.
#' Each step's equation follows, with when and why it is used.
#'
#' \strong{1. Apparent prevalence.} The pooled proportion of test
#' positives across all clusters,
#'
#' \deqn{\hat{p} = \frac{\sum_i x_i}{\sum_i n_i}}
#'
#' Every CI method below operates on \eqn{\hat{p}}. It is corrected for
#' test imperfection only at step 6, so the correction applies identically
#' to the point estimate and to both interval endpoints.
#'
#' \strong{2. Design effect (Kish).} Observations within a cluster are
#' positively correlated, so a clustered sample carries less information
#' than the same number of independent people. The Kish (1965) design
#' effect scales the variance by
#'
#' \deqn{D_{eff} = 1 + (\bar{n} - 1)\,\rho}
#'
#' where \eqn{\bar{n}} is the mean cluster size and \eqn{\rho} the
#' intra-cluster correlation (ICC).
#'
#' \emph{ICC supplied} (\code{icc} set): \eqn{\rho} is used directly. A
#' value below about 1.5e-8 is treated as 0. The ICC is ignored (with a
#' warning) and \eqn{D_{eff} = 1} when there is a single cluster, or when
#' every cluster has only one person, so no two people share a site.
#'
#' \emph{ICC estimated} (\code{icc = NULL}, the default): with two or more
#' clusters and an average cluster size above 1, the design effect is
#' estimated as the ratio of the observed between-cluster variance of the
#' site proportions to the variance expected under simple random sampling,
#'
#' \deqn{\widehat{D_{eff}} = \frac{\mathrm{Var}(\hat{p}_i)}{\;\overline{\hat{p}(1 - \hat{p}) / n_i}\;}, \qquad \hat{p}_i = \frac{x_i}{n_i}}
#'
#' floored at 1. The implied ICC is recovered by inverting the Kish formula,
#' \eqn{\hat{\rho} = (\widehat{D_{eff}} - 1) / (\bar{n} - 1)},
#' kept between 0 and 1; \eqn{D_{eff}} is then recomputed from that
#' \eqn{\hat{\rho}} so the returned \code{deff} and \code{icc_used} stay
#' mutually consistent. If nobody tests positive, or everybody does, the
#' expected variance is 0 and \eqn{D_{eff} = 1}. \eqn{D_{eff} = 1} also
#' when there is a single cluster (no variation between sites to measure)
#' or every cluster has only one person (no two people share a site).
#'
#' \strong{3. Effective sample size.}
#'
#' \deqn{n_{eff} = \frac{n_{total}}{D_{eff}}}
#'
#' the number of independent observations carrying the same information as
#' the clustered sample. Returned as \code{n_eff}.
#'
#' \strong{4. Finite-population correction.} Used only when \code{fpc_N}
#' (the population size \eqn{N}) is given. Testing a large share of the
#' population reduces the variance by the factor
#'
#' \deqn{f = \frac{N - n_{total}}{N - 1}}
#'
#' where \eqn{n_{total}} is the total number of people tested. Reducing the
#' variance by \eqn{f} is the same as dividing the sample size by \eqn{f},
#' so the function applies it to the effective sample size,
#'
#' \deqn{n_{eff,adj} = \frac{n_{eff}}{f}}
#'
#' which all three CI methods then use. Without \code{fpc_N},
#' \eqn{n_{eff,adj} = n_{eff}}. If everyone in the population was tested
#' (\code{fpc_N} equal to \code{sum(n)}), there is no sampling uncertainty
#' left, so the function stops with an error. \code{fpc_N} counts people,
#' not sites. If a study tested at nearly all of the sites in the
#' population, the interval will be wider than necessary, because the
#' correction does not remove the variation between sites.
#'
#' \strong{5. Confidence interval on apparent prevalence.} Let
#' \eqn{z = \Phi^{-1}(1 - \alpha/2)}, \eqn{\alpha = 1 - } \code{conf_level},
#' and \eqn{x_{eff} = \hat{p}\,n_{eff,adj}}, the number of positives scaled
#' to the adjusted effective sample size.
#'
#' \emph{\code{"wald"}} -- the normal-approximation interval,
#'
#' \deqn{\hat{p} \;\pm\; z \sqrt{\frac{\hat{p}(1 - \hat{p})}{n_{eff,adj}}}}
#'
#' with any end below 0 set to 0 and any end above 1 set to 1. For small
#' \eqn{n} or prevalence near 0 or 1 it can be too narrow (it covers the
#' true prevalence less often than the stated confidence level), and it has
#' zero width when \eqn{\hat{p} = 0} or \eqn{1}.
#'
#' \emph{\code{"clopper-pearson"}} -- the exact binomial interval,
#'
#' \deqn{L = B^{-1}\!\left(\tfrac{\alpha}{2};\; x_{eff},\; n_{eff,adj} - x_{eff} + 1\right)}
#' \deqn{U = B^{-1}\!\left(1 - \tfrac{\alpha}{2};\; x_{eff} + 1,\; n_{eff,adj} - x_{eff}\right)}
#'
#' where \eqn{B^{-1}(q; a, b)} is the \eqn{q}-th quantile of a beta
#' distribution with parameters \eqn{a} and \eqn{b}, and \eqn{L = 0} when
#' \eqn{\hat{p} = 0} and \eqn{U = 1} when \eqn{\hat{p} = 1}. For a simple
#' random sample, it covers the true prevalence at least as often as the
#' stated confidence level, so it is usually wider than the Wald interval.
#' Here \eqn{x_{eff}} does not need to be an integer, so the clustering
#' and population adjustments can be applied. When they are, "at least as
#' often as the stated confidence level" holds only approximately.
#'
#' \emph{\code{"agresti-coull"}} -- add \eqn{z^2/2} positive cases and
#' \eqn{z^2/2} negative cases (about 1.92 of each at 95 percent), then take
#' a Wald interval on the adjusted proportion,
#'
#' \deqn{\tilde{n} = n_{eff,adj} + z^2, \qquad \tilde{p} = \frac{x_{eff} + z^2/2}{\tilde{n}}}
#' \deqn{\tilde{p} \;\pm\; z \sqrt{\frac{\tilde{p}(1 - \tilde{p})}{\tilde{n}}}}
#'
#' with any end below 0 set to 0 and any end above 1 set to 1. The added
#' cases pull the interval towards 0.5, so unlike the Wald interval it does
#' not have zero width at \eqn{\hat{p} = 0} or \eqn{1}.
#'
#' \strong{6. Rogan-Gladen correction.} An imperfect test inflates apparent
#' prevalence through false positives and deflates it through false
#' negatives. Inverting \eqn{p_{app} = p\,Se + (1 - p)(1 - Sp)} gives true
#' prevalence,
#'
#' \deqn{p = \frac{p_{app} - (1 - Sp)}{Se + Sp - 1}}
#'
#' This same formula is applied to \eqn{\hat{p}} and to both CI endpoints.
#' A corrected value below 0 is then set to 0, and one above 1 is set to 1.
#' When \eqn{Se = Sp = 1} it leaves the values unchanged.
#'
#' The denominator \eqn{Se + Sp - 1} must be positive, as the correction
#' divides by it. Thus, as it approaches 0, the corrected interval becomes
#' very wide, and the function warns when it is below 0.1.
#'
#' If the whole corrected interval lies below 0, both ends are set to 0 and
#' the interval has zero width (likewise above 1). This happens when there
#' are fewer positives than false positives alone would produce at the
#' given specificity (or more positives than the sensitivity allows). The
#' function warns that this zero width is not real precision.
#'
#' \strong{Margin of error.} From the corrected point estimate and
#' endpoints,
#'
#' \deqn{\mathrm{moe} = \frac{ci_{upper} - ci_{lower}}{2}, \quad
#'       \mathrm{moe\_lower} = p - ci_{lower}, \quad
#'       \mathrm{moe\_upper} = ci_{upper} - p}
#'
#' For \code{"wald"} the interval is usually symmetric and all three
#' coincide -- unless an end is moved to 0 or 1 (very low or high
#' prevalence), which makes it asymmetric. For \code{"clopper-pearson"}
#' and \code{"agresti-coull"} the interval is asymmetric by design:
#' \code{moe} is only the average of the two half-widths, so \code{moe_lower}
#' and \code{moe_upper} should be reported together. For \code{"wald"}, the
#' function prints a message when an end is moved to 0 or 1.
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley.
#'
#' Cochran, W. G. (1977) \emph{Sampling Techniques}, 3rd ed. Wiley.
#' (Finite-population correction.)
#'
#' Clopper, C. J. & Pearson, E. S. (1934) The use of confidence or fiducial
#' limits illustrated in the case of the binomial. \emph{Biometrika}
#' \strong{26}(4), 404-413. \doi{10.1093/biomet/26.4.404}
#'
#' Agresti, A. & Coull, B. A. (1998) Approximate is better than "exact" for
#' interval estimation of binomial proportions. \emph{The American
#' Statistician} \strong{52}(2), 119-126. \doi{10.1080/00031305.1998.10480550}
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. \emph{American Journal of Epidemiology}
#' \strong{107}(1), 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{prevalence}{Point estimate of true prevalence (Rogan-Gladen
#'     corrected)}
#'   \item{ci_lower}{Lower confidence limit on the true-prevalence scale}
#'   \item{ci_upper}{Upper confidence limit on the true-prevalence scale}
#'   \item{moe}{Half-width of the interval, \code{(ci_upper - ci_lower) / 2},
#'     the average of \code{moe_lower} and \code{moe_upper}. The two are
#'     equal only when the interval is symmetric; otherwise report them both
#'     (see Details).}
#'   \item{moe_lower}{\code{prevalence - ci_lower}: distance from point estimate
#'     to lower limit}
#'   \item{moe_upper}{\code{ci_upper - prevalence}: distance from point estimate
#'     to upper limit}
#'   \item{method}{CI method used (as supplied)}
#'   \item{n_total}{Total number of people tested, \code{sum(n)}}
#'   \item{n_eff}{Effective independent sample size before the
#'     finite-population correction: \code{n_total / deff}}
#'   \item{n_eff_adj}{Effective sample size the CI is actually built from:
#'     \code{n_eff} divided by \eqn{f}, the finite-population correction
#'     factor (see Details, step 4). Equals \code{n_eff} when \code{fpc_N}
#'     is \code{NULL}.}
#'   \item{conf_level}{Confidence level (as supplied)}
#'   \item{sensitivity}{Sensitivity (as supplied)}
#'   \item{specificity}{Specificity (as supplied)}
#'   \item{icc_used}{ICC used in the calculation. If \code{icc = NULL}, it is
#'     the value estimated from the data. If \code{icc} was supplied, it is
#'     that value, unless the supplied value was below about 1.5e-8 or was
#'     ignored (a single cluster, or every cluster has only one person), then
#'     it is 0.}
#'   \item{deff}{Design effect applied (1 when there is no clustering)}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}}
#'
#' @export
#'
#' @examples
#' # Single site / simple random sample
#' estimate_prevalence(x = 8, n = 50)
#'
#' # Multi-site clustered data -- ICC estimated from the data
#' estimate_prevalence(
#'   x = c(0, 4, 0, 22, 25, 16, 12, 8),
#'   n = c(60, 80, 70, 100, 40, 60, 50, 90)
#' )
#'
#' # Clopper-Pearson interval
#' estimate_prevalence(x = 3, n = 30, method = "clopper-pearson")
#'
#' # Imperfect diagnostic test
#' estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95)
#'
#' # Known ICC and population size
#' estimate_prevalence(x = c(5, 8, 3), n = c(40, 40, 40), icc = 0.05,
#'                     fpc_N = 2000)
estimate_prevalence <- function(x,
                                n,
                                sensitivity = 1,
                                specificity = 1,
                                conf_level  = 0.95,
                                icc         = NULL,
                                fpc_N       = NULL,
                                method      = "wald") {

  # ---- validate x and n ----
  if (is.logical(x) || is.logical(n))
    stop("`x` and `n` must be numeric, not logical (got class `",
         class(x)[1], "` for x, `", class(n)[1], "` for n). ",
         "Note: a plain `NA` is logical in R. Remove missing observations ",
         "before calling.")
  if (!is.numeric(x) || !is.numeric(n))
    stop("`x` and `n` must be numeric (got class `", class(x)[1], "` for x, ",
         "`", class(n)[1], "` for n).")
  if (length(x) == 0 || length(n) == 0)
    stop("`x` and `n` must each contain at least one value. ",
         "Supply at least one cluster's count and total.")
  # For each check below, `i` is the position of the first bad element.
  if (!all(is.finite(x))) {
    i <- which(!is.finite(x))[1]
    stop("`x` contains a missing or infinite value (found x[", i, "] = ",
         x[i], "). The number of positives at each cluster must be a finite ",
         "non-negative integer.")
  }
  if (!all(is.finite(n))) {
    i <- which(!is.finite(n))[1]
    stop("`n` contains a missing or infinite value (found n[", i, "] = ",
         n[i], "). The number tested at each cluster must be a finite positive ",
         "integer.")
  }
  if (length(x) != length(n))
    stop("`x` and `n` must have the same length ",
         "(got length(x) = ", length(x), ", length(n) = ", length(n), "). ",
         "Each value of `x` is the number of positives at one cluster, and ",
         "each value of `n` is the number tested at that cluster.")
  if (any(x < 0)) {
    i <- which(x < 0)[1]
    stop("`x` must be non-negative (found x[", i, "] = ", x[i], "). ",
         "Counts cannot be negative.")
  }
  if (any(n <= 0)) {
    i <- which(n <= 0)[1]
    stop("`n` must be positive for every cluster (found n[", i, "] = ", n[i], "). ",
         "Each cluster must have tested at least one person.")
  }
  if (any(x != floor(x))) {
    i <- which(x != floor(x))[1]
    stop("`x` must contain integers, not decimals ",
         "(found x[", i, "] = ", x[i], ").")
  }
  if (any(n != floor(n))) {
    i <- which(n != floor(n))[1]
    stop("`n` must contain integers, not decimals ",
         "(found n[", i, "] = ", n[i], ").")
  }
  if (any(x > n)) {
    i <- which(x > n)[1]
    stop("`x` cannot exceed `n` (found x[", i, "] = ", x[i],
         " > n[", i, "] = ", n[i], "). ",
         "The number of positives cannot be more than the number tested at ",
         "that cluster.")
  }

  # ---- validate single-value parameters ----
  # Each of the parameters below must be a single value.
  if (length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got class `",
         class(sensitivity)[1], "`, length ", length(sensitivity), "). ",
         "Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. ",
         "Pass 1 if the test never misses a case.")
  if (length(specificity) != 1 || !is.numeric(specificity))
    stop("`specificity` must be a single number in (0, 1] (got class `",
         class(specificity)[1], "`, length ", length(specificity), "). ",
         "Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. ",
         "Pass 1 if the test never gives a false positive.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number in (0, 1) (got class `",
         class(conf_level)[1], "`, length ", length(conf_level), ").")
  if (!is.null(icc) && (length(icc) != 1 || !is.numeric(icc)))
    stop("`icc` must be a single number in [0, 1] (got class `",
         class(icc)[1], "`, length ", length(icc), "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.character(method) || length(method) != 1)
    stop("`method` must be one of ",
         "'wald', 'clopper-pearson', or 'agresti-coull' (got class `",
         class(method)[1], "`, length ", length(method), ").")
  if (!method %in% c("wald", "clopper-pearson", "agresti-coull"))
    stop("`method` must be one of 'wald', 'clopper-pearson', or 'agresti-coull' ",
         "(got '", method, "').")

  # Check for NA, NaN and Inf before the range checks below.
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", sensitivity, ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a single finite number (got ", specificity, ").")
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")
  if (!is.null(icc) && !is.finite(icc))
    stop("`icc` must be a single finite number (got ", icc, "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")

  if (sensitivity < 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, "). ",
         "`sensitivity` is a diagnostic probability and cannot be ",
         if (sensitivity < 0) "negative." else "greater than 1.")
  if (sensitivity == 0)
    stop("`sensitivity` must be in (0, 1] (got 0). ",
         "A sensitivity of 0 means everyone who truly has the condition ",
         "tests negative.")
  if (specificity < 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, "). ",
         "`specificity` is a diagnostic probability and cannot be ",
         if (specificity < 0) "negative." else "greater than 1.")
  if (specificity == 0)
    stop("`specificity` must be in (0, 1] (got 0). ",
         "A specificity of 0 means everyone who truly does not have the ",
         "condition tests positive.")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("`sensitivity` + `specificity` must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "The correction divides by `sensitivity` + `specificity` - 1, so a sum ",
         "of 1 or less cannot be corrected.")
  if (correction < 0.1)
    warning("`sensitivity` + `specificity` = ", round(sensitivity + specificity, 4),
            ", which is very close to 1. Correcting for this much test error ",
            "makes the prevalence estimate and its confidence interval very ",
            "wide, and small changes in the assumed sensitivity or specificity ",
            "will change them a lot.")
  if (!is.null(icc) && (icc < 0 || icc > 1))
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "`icc` is a correlation and cannot be ",
         if (icc < 0) "negative. " else "greater than 1. ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be in (0, 1) (got ", conf_level, "). ",
         "For example, use 0.95 for a 95% confidence interval.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || length(fpc_N) != 1 || !is.finite(fpc_N) ||
      fpc_N < 1 || fpc_N != floor(fpc_N)))
    stop("`fpc_N` must be a single finite positive integer (got ",
         if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
         else if (length(fpc_N) != 1) paste0("length = ", length(fpc_N))
         else fpc_N,
         "). `fpc_N` is the total population size. ",
         "Set `fpc_N = NULL` to skip the finite-population correction.")

  # A supplied icc below this tiny threshold (e.g. 1e-12, a rounding
  # leftover) is treated as exactly 0. There is then no "ignored" warning,
  # and deff stays exactly 1.
  if (!is.null(icc) && icc < sqrt(.Machine$double.eps))
    icc <- 0

  n_clusters <- length(n)
  n_total    <- sum(n)

  if (!is.null(fpc_N) && fpc_N < n_total)
    stop("`fpc_N` (", fpc_N, ") is smaller than the total sample size (",
         n_total, "): you cannot test more people than there are in the ",
         "population.")
  if (!is.null(fpc_N) && fpc_N == n_total)
    stop("`fpc_N` (", fpc_N, ") equals the total sample size: the whole ",
         "population was tested (a census), so there is no sampling ",
         "uncertainty and no confidence interval to compute. ",
         "Set `fpc_N = NULL` if no FPC is needed.")

  p_hat <- sum(x) / n_total   # apparent prevalence

  # -----------------------------------------------------------------
  # Design effect / ICC (Kish formula, Module 5)
  # -----------------------------------------------------------------
  # TODO(review): cluster-size convention for the Kish design effect.
  # We use the arithmetic mean cluster size n_bar = mean(n), which matches
  # the workshop (Module 5, slide "Why is the ICC useful?": "n_bar =
  # average cluster size"). The classical Kish deff for UNEQUAL clusters
  # uses the size-weighted mean sum(n^2)/sum(n) instead, which is larger
  # and inflates deff more. With a supplied `icc` and very unequal
  # clusters the two diverge substantially. DECISION NEEDED: match the
  # lecture (mean) or match standard Kish (weighted)? Flagged for the
  # team's statistical review -- do not "fix" silently.
  n_bar <- mean(n)

  if (is.null(icc)) {
    if (n_clusters < 2 || n_bar == 1) {
      # Single cluster or all clusters of size 1 -- Kish denominator is 0.
      icc_used <- 0
      deff     <- 1
    } else {
      p_i     <- x / n
      var_obs <- stats::var(p_i)
      # TODO(review): which prevalence goes into Var_SRS? We use the pooled
      # p_hat = sum(x) / sum(n). The workshop's worked example (Module 5,
      # "The Design Effect - worked example", Deff = 17.73) uses the
      # unweighted mean of the site prevalences instead. On that example's
      # data the pooled version gives Deff = 20.43; the mean-of-sites version
      # reproduces the slide (17.73 from the slide's rounded site values,
      # 17.94 from the exact counts). DECISION NEEDED -- do not change
      # silently (see test EP-C-2).
      var_srs <- mean(p_hat * (1 - p_hat) / n)

      deff <- if (var_srs > 0) var_obs / var_srs else 1
      deff <- max(deff, 1)

      icc_used <- (deff - 1) / (n_bar - 1)
      icc_used <- min(max(icc_used, 0), 1)
      deff     <- 1 + (n_bar - 1) * icc_used  # keep pair mutually consistent
    }
  } else if (n_clusters < 2 || n_bar == 1) {
    # A supplied icc has no effect with a single cluster (or clusters all of
    # size 1): the Kish denominator is undefined, so fall back to SRS -- same
    # as the icc = NULL path above.
    if (icc > 0)
      warning("`icc` = ", icc, " was ignored: the design effect needs a ",
              "cluster structure (>= 2 clusters, mean size > 1). ",
              "`icc_used` is reported as 0.")
    icc_used <- 0
    deff     <- 1
  } else {
    icc_used <- icc
    deff     <- 1 + (n_bar - 1) * icc_used
  }

  n_eff <- n_total / deff

  # Finite-population correction. The FPC is a property of the sampling
  # fraction of the units actually drawn -- the collected count n_total,
  # NOT the design-effect-adjusted n_eff (Cochran 1977 sec. 2.8; Kish
  # 1965). It is a separate adjustment from the design effect: deff
  # measures the inefficiency of the design, the FPC measures how much of
  # the population was observed. So the factor uses n_total, matching
  # design_precision() / design_threshold(). (This treats the population
  # as finite in individuals -- `fpc_N` is a headcount. A study that
  # sampled nearly all *sites* would also shrink the between-cluster
  # variance, which this single-FPC shortcut does not separately model.)
  # FPC factor (1 when fpc_N is NULL); fpc_N > n_total is guaranteed above.
  fpc <- if (!is.null(fpc_N)) sqrt((fpc_N - n_total) / (fpc_N - 1)) else 1

  # -----------------------------------------------------------------
  # Confidence interval on apparent prevalence
  # All three methods use n_eff_adj = n_eff / fpc^2, which collapses to
  # n_eff when there is no FPC (fpc = 1). This is the variance-equivalent
  # simple-random-sample size: sqrt(p*(1-p)/n_eff_adj) == sqrt(p*(1-p)/n_eff)*fpc.
  # -----------------------------------------------------------------
  z         <- stats::qnorm(1 - (1 - conf_level) / 2)
  n_eff_adj <- n_eff / (fpc^2)   # incorporates both Deff and FPC
  alpha     <- 1 - conf_level

  if (method == "wald") {
    se        <- sqrt(p_hat * (1 - p_hat) / n_eff_adj)
    ci_lo_app <- max(p_hat - z * se, 0)
    ci_hi_app <- min(p_hat + z * se, 1)

  } else if (method == "clopper-pearson") {
    x_eff <- p_hat * n_eff_adj   # effective successes (continuous)

    ci_lo_app <- if (p_hat == 0) 0 else
      stats::qbeta(alpha / 2,     x_eff,     n_eff_adj - x_eff + 1)
    ci_hi_app <- if (p_hat == 1) 1 else
      stats::qbeta(1 - alpha / 2, x_eff + 1, n_eff_adj - x_eff)

    ci_lo_app <- max(ci_lo_app, 0)
    ci_hi_app <- min(ci_hi_app, 1)

  } else {   # agresti-coull
    x_eff     <- p_hat * n_eff_adj
    n_tilde   <- n_eff_adj + z^2
    p_tilde   <- (x_eff + z^2 / 2) / n_tilde
    se_tilde  <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)

    ci_lo_app <- max(p_tilde - z * se_tilde, 0)
    ci_hi_app <- min(p_tilde + z * se_tilde, 1)
  }

  # -----------------------------------------------------------------
  # Rogan-Gladen correction: apparent -> true prevalence
  # -----------------------------------------------------------------
  rg <- function(p) .rogan_gladen(p, sensitivity, specificity)

  prevalence <- max(0, min(1, rg(p_hat)))
  ci_lower   <- max(0, min(1, rg(ci_lo_app)))
  ci_upper   <- max(0, min(1, rg(ci_hi_app)))

  moe       <- (ci_upper - ci_lower) / 2
  moe_lower <- prevalence - ci_lower
  moe_upper <- ci_upper - prevalence

  # Rogan-Gladen overshoot: the apparent-scale CI has real width, but after
  # correcting for an imperfect test both endpoints map outside [0, 1] and
  # clamp to the same boundary, so the corrected interval collapses and
  # `moe` reads as 0 -- false precision, not a genuinely exact estimate.
  # (Distinct from the Wald interval legitimately being [0, 0] at x = 0
  # with sensitivity = specificity = 1, where the apparent CI is already
  # degenerate.)
  eps <- .Machine$double.eps^0.5
  if ((ci_hi_app - ci_lo_app) > eps &&
      (ci_upper - ci_lower) < eps &&
      (prevalence == 0 || prevalence == 1))
    warning("Rogan-Gladen overshoot: after correcting for test error, the ",
            "estimate and both ends of the confidence interval are all cut to ",
            prevalence, ", because the corrected interval lies entirely ",
            "outside 0 to 1 for this sensitivity and specificity. `moe` is ",
            "reported as 0, but that is not real precision -- use a more ",
            "accurate test or a larger sample.")

  # Only for "wald": its interval is normally symmetric, so one with an
  # endpoint cut at 0 or 1 is lopsided unexpectedly. Clopper-Pearson and
  # Agresti-Coull are asymmetric by design, and the help page says so.
  if (method == "wald" && moe > 0 && abs(moe_lower - moe_upper) > eps)
    message(method, " CI is asymmetric: moe_lower = ", round(moe_lower, 4),
            ", moe_upper = ", round(moe_upper, 4),
            ". moe = ", round(moe, 4), " is the average half-width; ",
            "report moe_lower and moe_upper separately.")

  list(
    prevalence  = prevalence,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    moe         = moe,
    moe_lower   = moe_lower,
    moe_upper   = moe_upper,
    method      = method,
    n_total     = n_total,
    n_eff       = n_eff,
    n_eff_adj   = n_eff_adj,
    conf_level  = conf_level,
    sensitivity = sensitivity,
    specificity = specificity,
    icc_used    = icc_used,
    deff        = deff,
    fpc_N       = fpc_N
  )
}
