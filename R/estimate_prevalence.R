#' Estimate prevalence and its precision
#'
#' @description
#' Analysis function: given observed counts, returns a prevalence point
#' estimate and confidence interval. Three CI methods are available:
#'
#' \describe{
#'   \item{`"wald"` (default)}{Symmetric interval: \eqn{\hat{p} \pm z \cdot SE}.
#'     Fast and familiar; can extend below 0 or above 1 near the boundaries
#'     (clamped to \[0, 1\]).}
#'   \item{`"clopper-pearson"`}{Exact binomial interval via the beta distribution.
#'     Asymmetric; conservative (guaranteed coverage). Preferred for small
#'     samples or extreme prevalences.}
#'   \item{`"agresti-coull"`}{Adjusted-proportion interval. Asymmetric; better
#'     coverage than Wald for moderate n, less conservative than
#'     Clopper-Pearson.}
#' }
#'
#' All methods account for imperfect diagnostic tests (Rogan-Gladen
#' correction), clustered sampling (Kish design effect), and
#' finite-population corrections.
#'
#' If `icc` is not supplied, it is estimated from the data as the ratio
#' of observed-to-expected variance across clusters, rather than assumed
#' to be 0 -- assuming independence would understate uncertainty for any
#' genuinely clustered survey.
#'
#' @param x Integer vector of positive counts per cluster/site.
#' @param n Integer vector of total samples per cluster/site (same length as x).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (perfect
#'   test). Set below 1 to activate the Rogan-Gladen correction.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Confidence level; default 0.95.
#' @param icc Optional. Intra-cluster correlation. If `NULL` (default), ICC
#'   is estimated from the data. Set to `0` to force the SRS (no-clustering)
#'   case. Only relevant when `x` and `n` have more than one element.
#' @param fpc_N Optional. Total population size for a finite-population
#'   correction. `NULL` (default) = no FPC applied.
#' @param method CI method: `"wald"` (default), `"clopper-pearson"`, or
#'   `"agresti-coull"`. See Description. Clopper-Pearson and Agresti-Coull
#'   produce asymmetric intervals; a `message()` is emitted when they differ
#'   noticeably from the symmetric summary `moe`.
#'
#' @section Inputs and outputs:
#' \strong{Inputs} (function arguments):
#' \itemize{
#'   \item \code{x}, \code{n} -- \emph{required}. Equal-length integer
#'     vectors: per-cluster positive counts and per-cluster totals. One
#'     element each is a simple random sample; several elements is a
#'     clustered design.
#'   \item \code{sensitivity}, \code{specificity} -- diagnostic test
#'     characteristics. Both \code{1} (the default) means a perfect test and
#'     no Rogan-Gladen correction.
#'   \item \code{conf_level} -- confidence level (default \code{0.95}).
#'   \item \code{icc} -- clustering strength. \code{NULL} (default) estimates
#'     it from the data; \code{0} forces the simple-random-sample case; a
#'     value in \[0, 1\] fixes it.
#'   \item \code{fpc_N} -- total population size for the finite-population
#'     correction. \code{NULL} (default) applies none.
#'   \item \code{method} -- CI method: \code{"wald"} (default),
#'     \code{"clopper-pearson"}, or \code{"agresti-coull"}.
#' }
#' \strong{Outputs} (named elements of the returned list):
#' \itemize{
#'   \item \emph{Estimate} -- \code{prevalence}, \code{ci_lower},
#'     \code{ci_upper}: point estimate and confidence limits, all on the
#'     true-prevalence scale (Rogan-Gladen corrected).
#'   \item \emph{Precision} -- \code{moe}, \code{moe_lower},
#'     \code{moe_upper}: interval half-width and the two one-sided distances
#'     (all equal for \code{"wald"}).
#'   \item \emph{Design quantities} -- \code{n_total}, \code{n_eff},
#'     \code{deff}, \code{icc_used}: what the clustering and FPC adjustments
#'     resolved to.
#'   \item \emph{Echoed inputs} -- \code{method}, \code{conf_level},
#'     \code{sensitivity}, \code{specificity}, \code{fpc_N}: returned
#'     unchanged so a result is self-describing.
#' }
#' Field-by-field definitions are under \strong{Value}, below.
#'
#' @details
#' The interval is built in six stages: apparent prevalence, a design
#' effect for clustering, an effective sample size, an optional
#' finite-population correction, a confidence interval by the chosen
#' method, and finally the Rogan-Gladen correction for an imperfect test.
#' Every equation used, and the reason it is used, follows. The
#' \code{methods} vignette (\code{vignette("methods", "MMSpower")}) gives
#' the same material with derivations and a worked example.
#'
#' \strong{1. Apparent prevalence.} The pooled proportion of test
#' positives across all clusters,
#'
#' \deqn{\hat{p} = \frac{\sum_i x_i}{\sum_i n_i}}
#'
#' Every CI method below operates on \eqn{\hat{p}}. It is corrected for
#' test imperfection only at stage 6, so the correction applies identically
#' to the point estimate and to both interval endpoints.
#'
#' \strong{2. Design effect (Kish).} Observations within a cluster are
#' positively correlated, so a clustered sample carries less information
#' than its nominal size. The Kish (1965) design effect scales the
#' variance by
#'
#' \deqn{D_{eff} = 1 + (\bar{n} - 1)\,\rho}
#'
#' where \eqn{\bar{n}} is the mean cluster size and \eqn{\rho} the
#' intra-cluster correlation (ICC). Taking \eqn{\rho = 0} would understate
#' uncertainty for any genuinely clustered survey.
#'
#' \emph{ICC supplied} (\code{icc} set): \eqn{\rho} is used directly.
#'
#' \emph{ICC estimated} (\code{icc = NULL}, the default): with two or more
#' clusters of size > 1, the design effect is estimated as the ratio of the
#' observed between-cluster variance of the site proportions to the
#' variance expected under simple random sampling,
#'
#' \deqn{\widehat{D_{eff}} = \frac{\mathrm{Var}(\hat{p}_i)}{\;\overline{\hat{p}(1 - \hat{p}) / n_i}\;}, \qquad \hat{p}_i = \frac{x_i}{n_i}}
#'
#' floored at 1. The implied ICC is recovered by inverting the Kish
#' formula, \eqn{\hat{\rho} = (\widehat{D_{eff}} - 1) / (\bar{n} - 1)},
#' clamped to \[0, 1\]; \eqn{D_{eff}} is then recomputed from the clamped
#' \eqn{\hat{\rho}} so the returned \code{deff} and \code{icc_used} stay
#' mutually consistent. Estimating \eqn{\rho} rather than assuming it avoids
#' silently analysing a clustered design as if it were independent. With a
#' single cluster, or clusters all of size 1, the Kish denominator is
#' undefined and the SRS case (\eqn{D_{eff} = 1}) is used.
#'
#' \strong{3. Effective sample size.}
#'
#' \deqn{n_{eff} = \frac{n_{total}}{D_{eff}}}
#'
#' the number of independent observations carrying the same information as
#' the clustered sample. Returned as \code{n_eff}.
#'
#' \strong{4. Finite-population correction.} When the sample is a
#' non-trivial fraction of a population of known size \eqn{N}
#' (\code{fpc_N}), sampling without replacement reduces the sampling
#' variance by the factor
#'
#' \deqn{f = \frac{N - n_{eff}}{N - 1}}
#'
#' Rather than multiply each method's variance by \eqn{f}, the correction
#' is folded into a single adjusted sample size
#'
#' \deqn{n_{eff,adj} = \frac{n_{eff}}{f}}
#'
#' shared by all three CI methods (with no FPC, \eqn{f = 1} and
#' \eqn{n_{eff,adj} = n_{eff}}). This makes the methods agree in the
#' large-sample limit while each keeps its own boundary behaviour.
#' \eqn{n_{eff,adj} \to \infty} as \eqn{n_{eff} \to N}; the function stops
#' before that point, where the sampling variance would be zero.
#'
#' \strong{5. Confidence interval on apparent prevalence.} Let
#' \eqn{z = \Phi^{-1}(1 - \alpha/2)}, \eqn{\alpha = 1 - } \code{conf_level},
#' and \eqn{x_{eff} = \hat{p}\,n_{eff,adj}}.
#'
#' \emph{\code{"wald"}} -- the normal-approximation interval,
#'
#' \deqn{\hat{p} \;\pm\; z \sqrt{\frac{\hat{p}(1 - \hat{p})}{n_{eff,adj}}}}
#'
#' clamped to \[0, 1\]. Fast and familiar, but under-covers for small
#' \eqn{n} or prevalence near 0 or 1, and gives a zero-width interval when
#' \eqn{\hat{p} = 0} or \eqn{1}.
#'
#' \emph{\code{"clopper-pearson"}} -- the exact binomial interval, via the
#' beta-quantile identity,
#'
#' \deqn{L = B^{-1}\!\left(\tfrac{\alpha}{2};\; x_{eff},\; n_{eff,adj} - x_{eff} + 1\right)}
#' \deqn{U = B^{-1}\!\left(1 - \tfrac{\alpha}{2};\; x_{eff} + 1,\; n_{eff,adj} - x_{eff}\right)}
#'
#' with \eqn{L = 0} when \eqn{\hat{p} = 0} and \eqn{U = 1} when
#' \eqn{\hat{p} = 1}. Guarantees at least nominal coverage, so it is the
#' safe choice for small samples or extreme prevalence, at the cost of
#' being conservative. Continuous \eqn{x_{eff}} replaces an integer count so
#' the clustering and FPC adjustments carry through.
#'
#' \emph{\code{"agresti-coull"}} -- add \eqn{z^2} pseudo-observations, then
#' take a Wald interval on the adjusted proportion,
#'
#' \deqn{\tilde{n} = n_{eff,adj} + z^2, \qquad \tilde{p} = \frac{x_{eff} + z^2/2}{\tilde{n}}}
#' \deqn{\tilde{p} \;\pm\; z \sqrt{\frac{\tilde{p}(1 - \tilde{p})}{\tilde{n}}}}
#'
#' clamped to \[0, 1\]. Recovers most of Clopper-Pearson's coverage gain
#' over Wald without the full conservatism; a reasonable default for
#' moderate \eqn{n}.
#'
#' \strong{6. Rogan-Gladen correction.} An imperfect test inflates apparent
#' prevalence through false positives and deflates it through false
#' negatives. Inverting \eqn{p_{app} = p\,Se + (1 - p)(1 - Sp)} gives true
#' prevalence,
#'
#' \deqn{p = \frac{p_{app} - (1 - Sp)}{Se + Sp - 1}}
#'
#' This affine map is applied identically to \eqn{\hat{p}} and to both CI
#' endpoints, each result then clamped to \[0, 1\]. When \eqn{Se = Sp = 1}
#' it is the identity. The denominator \eqn{Se + Sp - 1} must be positive (a
#' test better than chance); as it approaches 0 the correction inflates the
#' estimate and its interval without bound, and the function warns below
#' 0.1.
#'
#' \strong{Margin of error.} From the corrected point estimate and
#' endpoints,
#'
#' \deqn{\mathrm{moe} = \frac{ci_{upper} - ci_{lower}}{2}, \quad
#'       \mathrm{moe\_lower} = p - ci_{lower}, \quad
#'       \mathrm{moe\_upper} = ci_{upper} - p}
#'
#' For \code{"wald"} the interval is usually symmetric and all three
#' coincide -- unless an endpoint is clamped to \[0, 1\] (very low or high
#' prevalence), which makes it asymmetric too. For \code{"clopper-pearson"}
#' and \code{"agresti-coull"} the interval is asymmetric by construction:
#' \code{moe} is only the average of the two half-widths, so \code{moe_lower}
#' and \code{moe_upper} should be reported together. A \code{message()} is
#' emitted whenever the two half-widths differ by more than 10\% of
#' \code{moe}, whichever method produced it.
#'
#' @section Equations and sources:
#' Mostly direct workshop material (MMS-SD Study Design Workshop,
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}):
#' \itemize{
#'   \item \emph{Apparent prevalence} \eqn{\hat p = x/n} and the
#'     \emph{Wald interval} \eqn{\hat p \pm z_{1-\alpha/2}\sqrt{\hat p
#'     (1-\hat p)/n}} -- Module 1 "Sampling from a population", slides "The
#'     Wald confidence interval" / "Putting it all together" (lecture
#'     slides p. 11).
#'   \item \emph{Design effect}
#'     \eqn{D_{eff} = \mathrm{Var}_{obs}/\mathrm{Var}_{SRS}},
#'     \emph{effective sample size} \eqn{N_{eff} = N/D_{eff}},
#'     \eqn{D_{eff} = 1 + (\bar n - 1)\,r} and
#'     \eqn{r = (D_{eff}-1)/(\bar n - 1)} -- Module 5 "Dealing with
#'     over-dispersion in multi-cluster studies", slides "The Design
#'     Effect", "The effective sample size", "The intra-cluster
#'     correlation coefficient" (lecture slides pp. 4-5).
#'   \item \emph{Clustered (generalised) Wald interval}
#'     \eqn{\hat p \pm z_{1-\alpha/2}\sqrt{\hat p(1-\hat p)/N \cdot
#'     D_{eff}}} -- Module 5, slide "How can we design multi-cluster
#'     studies?" (p. 5).
#'   \item \emph{Over-dispersion check} (sites expected within
#'     \eqn{\hat p \pm \sqrt{\hat p(1-\hat p)/n_i}}; more than ~10\% of
#'     sites outside implies over-dispersion) -- Module 5, slide
#'     "Detecting over-dispersion" (p. 4).
#'   \item \emph{Rogan-Gladen correction}
#'     \eqn{\hat p_{true} = (\hat p_{app} - (1-Sp))/(Se + Sp - 1)} and its
#'     delta-method variance, \emph{Clopper-Pearson} and \emph{Agresti-Coull}
#'     intervals, and the \emph{finite-population correction} -- \strong{not}
#'     in the workshop; see the references below and Cochran (1977) for the
#'     FPC.
#' }
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 1 (sampling / Wald interval) and 5
#' (ICC / design effect). \url{https://mrc-ide.github.io/MMS-SD_workshop/}
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
#' @return A named list. The following fields are always present:
#'   \item{prevalence}{Point estimate of true prevalence (Rogan-Gladen corrected)}
#'   \item{ci_lower}{Lower confidence limit on the true-prevalence scale}
#'   \item{ci_upper}{Upper confidence limit on the true-prevalence scale}
#'   \item{moe}{Half-width of the interval: \code{(ci_upper - ci_lower) / 2}.
#'     For \code{"wald"} it equals both \code{moe_lower} and \code{moe_upper}
#'     and fully describes the precision -- unless an endpoint was clamped to
#'     \[0, 1\], which makes even the Wald interval asymmetric. For
#'     \code{"clopper-pearson"} and \code{"agresti-coull"} the interval is
#'     asymmetric by construction: \code{moe} is only the average of the two
#'     half-widths and does not describe either side. Whenever the interval is
#'     asymmetric, report \code{moe_lower} and \code{moe_upper} together.}
#'   \item{moe_lower}{\code{prevalence - ci_lower}: distance from point estimate
#'     to lower limit}
#'   \item{moe_upper}{\code{ci_upper - prevalence}: distance from point estimate
#'     to upper limit}
#'   \item{method}{CI method used (as supplied)}
#'   \item{n_total}{Total samples across all clusters}
#'   \item{n_eff}{Effective independent sample size before the FPC:
#'     \code{n_total / deff}}
#'   \item{n_eff_adj}{Effective sample size the CI is actually built from:
#'     \code{n_eff} divided by the squared FPC factor (equals \code{n_eff}
#'     when \code{fpc_N} is \code{NULL})}
#'   \item{conf_level}{Confidence level (as supplied)}
#'   \item{sensitivity}{Sensitivity (as supplied)}
#'   \item{specificity}{Specificity (as supplied)}
#'   \item{icc_used}{ICC applied: estimated from data if \code{icc = NULL},
#'     else as supplied}
#'   \item{deff}{Design effect applied (1 for SRS)}
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
#' # Exact binomial interval (better for small samples)
#' estimate_prevalence(x = 3, n = 30, method = "clopper-pearson")
#'
#' # Imperfect diagnostic test
#' estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95)
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
         "Note: plain `NA` is logical in R -- filter missing observations before calling, ",
         "or use NA_real_ if you need a typed NA placeholder.")
  if (!is.numeric(x) || !is.numeric(n))
    stop("`x` and `n` must be numeric vectors (got class `", class(x)[1], "` for x, ",
         "`", class(n)[1], "` for n). Supply integer or double counts.")
  if (length(x) == 0 || length(n) == 0)
    stop("`x` and `n` must be non-empty vectors. ",
         "Supply at least one cluster's count and total.")
  if (!all(is.finite(x)))
    stop("`x` contains NA, NaN, or infinite values. ",
         "All counts must be finite non-negative integers.")
  if (!all(is.finite(n)))
    stop("`n` contains NA, NaN, or infinite values. ",
         "All cluster totals must be finite positive integers.")
  if (length(x) != length(n))
    stop("`x` and `n` must have the same length ",
         "(got length(x) = ", length(x), ", length(n) = ", length(n), "). ",
         "Each element of `x` is the positive count for one cluster and each ",
         "element of `n` is that cluster's total.")
  if (any(x != floor(x)))
    stop("`x` must contain whole numbers -- counts cannot be fractional ",
         "(found x[", which(x != floor(x))[1], "] = ", x[which(x != floor(x))[1]], ").")
  if (any(n != floor(n)))
    stop("`n` must contain whole numbers -- sample sizes cannot be fractional ",
         "(found n[", which(n != floor(n))[1], "] = ", n[which(n != floor(n))[1]], ").")
  if (any(n <= 0))
    stop("`n` must be positive for every cluster ",
         "(found n[", which(n <= 0)[1], "] = ", n[which(n <= 0)[1]], "). ",
         "A cluster with zero or negative total is undefined.")
  if (any(x < 0))
    stop("`x` must be non-negative ",
         "(found x[", which(x < 0)[1], "] = ", x[which(x < 0)[1]], "). ",
         "Counts cannot be negative.")
  if (any(x > n))
    stop("`x` cannot exceed `n` ",
         "(found x[", which(x > n)[1], "] = ", x[which(x > n)[1]],
         " > n[", which(x > n)[1], "] = ", n[which(x > n)[1]], "). ",
         "Positive counts cannot exceed the total tested per cluster.")

  # ---- validate scalar parameters ----
  if (is.null(sensitivity) || length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got ",
         if (is.null(sensitivity)) "NULL"
         else if (!is.numeric(sensitivity)) paste0("class `", class(sensitivity)[1], "`")
         else paste0("length = ", length(sensitivity)), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (is.null(specificity) || length(specificity) != 1 || !is.numeric(specificity))
    stop("`specificity` must be a single number in (0, 1] (got ",
         if (is.null(specificity)) "NULL"
         else if (!is.numeric(specificity)) paste0("class `", class(specificity)[1], "`")
         else paste0("length = ", length(specificity)), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number (got ",
         if (!is.numeric(conf_level)) paste0("class `", class(conf_level)[1], "`")
         else paste0("length = ", length(conf_level)), "). Use, e.g., 0.95.")
  if (!is.null(icc) && (length(icc) != 1 || !is.numeric(icc)))
    stop("`icc` must be a single number or NULL (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`")
         else paste0("length = ", length(icc)), "). ",
         "Set `icc = NULL` to estimate ICC from the data.")
  if (!is.null(fpc_N) && length(fpc_N) != 1)
    stop("`fpc_N` must be a single number or NULL (got length = ", length(fpc_N), "). ",
         "Set `fpc_N = NULL` to skip the finite-population correction.")
  if (!is.character(method) || length(method) != 1)
    stop("`method` must be a single character string: ",
         "'wald', 'clopper-pearson', or 'agresti-coull' (got class `",
         class(method)[1], "`, length ", length(method), ").")
  if (!method %in% c("wald", "clopper-pearson", "agresti-coull"))
    stop("`method` must be one of 'wald', 'clopper-pearson', or 'agresti-coull' ",
         "(got '", method, "').")

  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, "). ",
         "Use, e.g., 0.95 for a 95% confidence interval.")
  if (!is.null(icc) && !is.finite(icc))
    stop("`icc` must be a finite number (got ", icc, "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.null(icc) && (icc < 0 || icc > 1))
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "Values outside this range imply negative within-cluster variance, ",
         "which is not possible. To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be a finite positive number representing the total population size ",
         "(got ", fpc_N, "). Set `fpc_N = NULL` to skip the finite-population correction.")
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a finite number (got ", sensitivity, ").")
  if (sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, "). ",
         "A sensitivity of 0 means the test never detects true positives.")
  if (!is.finite(specificity))
    stop("`specificity` must be a finite number (got ", specificity, ").")
  if (specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, "). ",
         "A specificity of 0 means the test always returns a false positive.")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("sensitivity + specificity must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "With se + sp <= 1 the test performs at or below chance and true prevalence ",
         "is not identifiable from apparent prevalence.")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1 (correction = ", round(correction, 4), "). ",
            "The Rogan-Gladen adjustment is numerically unstable -- ",
            "prevalence estimates and CI will be unreliable.")

  n_clusters <- length(n)
  n_total    <- sum(n)

  if (!is.null(fpc_N) && fpc_N < n_total)
    stop("`fpc_N` = ", fpc_N, " is less than the total sample size = ", n_total, ". ",
         "The population must be at least as large as the sample. ",
         "Check your inputs, or set `fpc_N = NULL` to skip the FPC.")

  p_hat <- sum(x) / n_total   # apparent prevalence

  # -----------------------------------------------------------------
  # Design effect / ICC (Kish formula, Module 5)
  # -----------------------------------------------------------------
  n_bar <- mean(n)

  if (is.null(icc)) {
    if (n_clusters < 2 || n_bar == 1) {
      # Single cluster or all clusters of size 1 -- Kish denominator is 0.
      icc_used <- 0
      deff     <- 1
    } else {
      p_i     <- x / n
      var_obs <- stats::var(p_i)
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
    icc_used <- 0
    deff     <- 1
  } else {
    icc_used <- icc
    deff     <- 1 + (n_bar - 1) * icc_used
  }

  n_eff <- n_total / deff

  # Check after computing n_eff: if fpc_N <= n_eff the FPC collapses to 0
  # and the CI is undefined. For SRS (deff=1) this is fpc_N == n_total; for
  # clustered data n_eff < n_total so fpc_N = n_total is still valid.
  if (!is.null(fpc_N) && fpc_N <= n_eff)
    stop("`fpc_N` (", fpc_N, ") is <= the effective sample size n_eff (", round(n_eff, 2),
         "). Sampling variance is zero and the CI is undefined. For SRS this means ",
         "you have surveyed the entire population; set `fpc_N = NULL` if no FPC is needed.")

  # FPC factor (1 when fpc_N is NULL)
  fpc <- if (!is.null(fpc_N)) sqrt((fpc_N - n_eff) / (fpc_N - 1)) else 1

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

  # Fires for the asymmetric methods, and also for "wald" when an endpoint
  # has been clamped to [0, 1] (which breaks its usual symmetry).
  if (moe > 0 && abs(moe_lower - moe_upper) > 0.1 * moe)
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
