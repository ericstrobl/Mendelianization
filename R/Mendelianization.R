#' Run Mendelianization with optional outcome and covariance screening
#'
#' Computes variant-specific canonical coefficients and learning-aware
#' chi-square statistics from a q by m matrix of sample z- or t-statistics.
#' The function can remove redundant outcomes, estimate an adaptive central-null
#' covariance matrix, and calculate the Score of Mendelianism for lead variants.
#'
#' @param Zstat A q by m numeric matrix of sample z- or t-statistics, with rows
#'   representing variants and columns representing outcomes.
#' @param SoM Logical; whether to compute the Score of Mendelianism. Default:
#'   `FALSE`.
#' @param chr_pos If `SoM = TRUE`, a data frame with chromosome and base-pair
#'   position columns. Its rows must correspond to the rows of `Zstat`.
#' @param alpha Genome-wide significance threshold intended for SoM calculations.
#'   The current function body retains this argument but does not pass it to the
#'   SoM helper functions.
#' @param screen_outcomes Logical; whether to iteratively remove outcomes until
#'   the preliminary null-correlation condition number is at most
#'   `max_condition`. Default: `TRUE`.
#' @param screen_Gamma Logical; whether to estimate `Gamma` using the adaptive
#'   central-null estimator. If `FALSE`, the raw genome-wide cross-product
#'   estimate is used. Default: `TRUE`.
#' @param threshold_grid Numeric vector of candidate trimming proportions passed
#'   to `estimate_null_gamma()`.
#' @param max_condition Maximum permitted preliminary null-correlation condition
#'   number when `screen_outcomes = TRUE`. Default: `100`.
#'
#' @return A list with components:
#' \itemize{
#'   \item `Alpha`: raw canonical coefficients, with retained outcomes in rows
#'     and variants in columns.
#'   \item `Alpha_p`: canonical coefficients standardized using the diagonal of
#'     the precision matrix.
#'   \item `pval`: p-values from the learning-aware chi-square tests.
#'   \item `qform`: corresponding quadratic-form statistics.
#'   \item `SoMs`: Score-of-Mendelianism results when requested; otherwise an
#'     empty numeric vector.
#'   \item `chr_pos`: chromosome and position data sorted into genomic order.
#'   \item `Gamma`: estimated null covariance matrix for the retained outcomes.
#'   \item `Omega`: precision matrix used for phenotype learning and inference.
#'   \item `lambda`: selected trimming proportion, or zero when adaptive
#'     covariance screening is disabled.
#'   \item `selected_outcomes`: names of the retained outcomes.
#'   \item `excluded_outcomes`: names of outcomes removed during screening.
#'   \item `condition_number`: final null-correlation condition number from the
#'     outcome-screening step.
#' }
Mendelianization <- function(
    Zstat,
    SoM = FALSE,
    chr_pos = NULL,
    alpha = 5e-8,
    screen_outcomes = TRUE,
    screen_Gamma = TRUE,
    threshold_grid = seq(0, 0.20, length.out = 20),
    max_condition = 100
) {
  # Chromosome and position information are required to calculate localization.
  if (SoM && is.null(chr_pos)) {
    stop("chr_pos must be supplied when SoM = TRUE.")
  }
  
  # Put variants in genomic order and apply the same ordering to the z statistics.
  if (!is.null(chr_pos)) {
    chr_levels <- c(as.character(1:22), "X", "Y")
    ord <- order(
      match(as.character(chr_pos[[1]]), chr_levels),
      chr_pos[[2]]
    )
    chr_pos <- chr_pos[ord, , drop = FALSE]
    Zstat <- Zstat[ord, , drop = FALSE]
  }
  
  # Construct a preliminary covariance estimate for assessing outcome redundancy.
  Gamma_screen <- crossprod(Zstat) / nrow(Zstat)
  if (screen_outcomes){
    # Remove outcomes greedily until the condition-number requirement is met.
    selection <- select_outcomes(Gamma_screen, max_condition = max_condition)
    Zstat <- Zstat[,selection$index,drop = FALSE]
  } else{
    # Retain every outcome while still recording names and conditioning metadata.
    selection <- select_outcomes(Gamma_screen,max_condition = Inf)
  }
  
  # Either estimate the central-null covariance adaptively or use the raw estimate.
  if (screen_Gamma){
    Gamma_fit <- estimate_null_gamma(Zstat,threshold_grid = threshold_grid)
    Gamma <- Gamma_fit$Gamma 
  } else{
    Gamma <- crossprod(Zstat) / nrow(Zstat)
  }
  
  # Invert Gamma with a small diagonal stabilization, learn each variant-specific
  # phenotype, and calculate its learning-aware quadratic association statistic.
  Omega <- chol2inv(chol(Gamma + 1e-8 * diag(diag(Gamma))))
  Alpha <- Omega %*% t(Zstat)
  qform <- colSums(t(Zstat) * Alpha)
  pval <- 1 - pchisq(qform, df = ncol(Omega))
  
  # Standardize coefficients using the diagonal of the precision matrix.
  Alpha_p <- sweep(Alpha,1,sqrt(diag(Omega)),"/")
  
  # Optionally identify lead variants and quantify genome-wide localization.
  SoMs <- c()
  if (SoM) {
    leads <- find_leads(chr_pos, pval)
    SoMs <- Score_of_Mendelianism_tower(Alpha,Gamma,Zstat,leads,chr_pos)
  }
  
  # Return inference results, covariance quantities, and screening metadata.
  list(
    Alpha = Alpha,
    Alpha_p = Alpha_p,
    pval = pval,
    qform = qform,
    SoMs = SoMs,
    chr_pos = chr_pos,
    Gamma = Gamma,
    Omega = Omega,
    lambda = if (screen_Gamma) Gamma_fit$lambda else 0,
    selected_outcomes = selection$outcomes,
    excluded_outcomes = selection$removed,
    condition_number = selection$condition_number
  )
}
