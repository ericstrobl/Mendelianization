Mendelianization <- function(Zstat, SoM = FALSE, chr_pos = NULL, alpha = 5e-8) {
#' Computes canonical coefficients and associated p-values from a q × m matrix
#' of sample z- (or t-) statistics, with an option to compute a score of
#' Mendelianism for lead variants/towers.
#'
#' @param Zstat A q × m numeric matrix of sample z- (or t-) statistics.
#' @param SoM Logical; whether to compute the Score of Mendelianism. Default: FALSE.
#' @param chr_pos If \code{SoM = TRUE}, a data frame with two columns: chromosome
#'   and (increasing) base-pair position.
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
  
  Gamma <- crossprod(Zstat) / nrow(Zstat) # positive semi-definite
  Omega   <- chol2inv(chol(Gamma + 1e-8*diag(diag(Gamma)))) # small regularization to ensure invertibility
  Alpha <- Omega  %*% t(Zstat)                      # m x q
  qform <- colSums(t(Zstat) * Alpha) 
  pval = 1-pchisq(qform,df=ncol(Omega))
  Alpha_p <- sweep(Alpha, 1, sqrt(diag(Omega)), "/")
  
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
    SoMs = SoMs,
    meta  = list(Gamma = Gamma, Omega = Omega)
  )
}

