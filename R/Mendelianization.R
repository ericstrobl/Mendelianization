Mendelianization <- function(Zstat, SoM = FALSE, chr_pos = NULL, alpha = 5e-8) {
#' Computes canonical coefficients and associated p-values from a q × m matrix
#' of sample z- (or t-) statistics, with an option to compute a score of
#' Mendelianism for lead variants/towers.
#'
#' @param Zstat A q × m numeric matrix of sample z- (or t-) statistics.
#' @param SoM Logical; whether to compute the Score of Mendelianism. Default: FALSE.
#' @param chr_pos If \code{SoM = TRUE}, a data frame with two columns: chromosome
#'   and base-pair position, where rows match those of Zstat.
#' @param alpha Genome-wide significance threshold (used by SoM).
#'
#' @returns A list with components:
#' \itemize{
#'   \item \code{Alpha}: raw canonical coefficients (m × q).
#'   \item \code{Alpha_p}: interpretable canonical coefficients.
#'   \item \code{pval}: p-values from the chi-square test with m d.f.
#'   \item \code{SoMs}: (optional) score(s) of Mendelianism if computed.
#'   \item \code{meta}: list with \code{Gamma} and \code{Omega} matrices.
#' }

  if (!is.null(chr_pos)) { # sort chr_pos if not already sorted
    chr_levels <- c(as.character(1:22), "X", "Y")
    ord <- order(
      match(as.character(chr_pos[[1]]), chr_levels),
      chr_pos[[2]]
    )
    chr_pos <- chr_pos[ord, , drop = FALSE]
    Zstat   <- Zstat[, ord, drop = FALSE] # make sure Zstat has the same new ordering
  }
  
  Gamma <- crossprod(Zstat) / nrow(Zstat) # positive semi-definite; Proposition 3 and Theorem 2 in the paper
  Omega   <- chol2inv(chol(Gamma + 1e-8*diag(diag(Gamma)))) # small regularization to ensure invertibility
  Alpha <- Omega  %*% t(Zstat) # raw coefficients                      
  qform <- colSums(t(Zstat) * Alpha) # hypothesis test of Expression (3)
  pval = 1-pchisq(qform,df=ncol(Omega)) # asymptotic distribution is chi^2_m
  Alpha_p <- sweep(Alpha, 1, sqrt(diag(Omega)), "/") # Proposition 4
  
  ### COMPUTE SCORE OF MENDELIANISM
  SoMs = c()
  if (SoM){
    leads = find_leads(chr_pos,pval) # find towers with Q statistics
    SoMs = Score_of_Mendelianism_tower(Alpha,Gamma,Z,leads,chr_pos)
  }
  
  list(
    Alpha = Alpha,  
    Alpha_p =  Alpha_p,
    pval = pval, 
    qform = qform,
    SoMs = SoMs,
    chr_pos = chr_pos, #new positions
    Gamma = Gamma,
    Omega = Omega
  )
}


