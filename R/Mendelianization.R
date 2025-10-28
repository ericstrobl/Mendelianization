Mendelianization <- function(Zstat, SoM = FALSE, chr_pos = NULL, alpha = 5e-8) {
  # Inputs:
  #   Zstat: q x m matrix of sample z-statistics (or t-statistics)
  #   SoM: whether or not to compute the score of Mendelianism, default: TRUE
  #   chr_pos: if SoM=TRUE, then data frame with one column containing chromosomes and another containing sorted positions (increasing) must be provided
  #   alpha: alpha threshold for genome-wide significance, default: 5e-8
  #
  # Outputs: list of
  #   Alpha: raw canonical coefficients
  #   Alpha_p: interpretable canonical coefficients (Proposition 3)
  #   pval: p-values of chi2-m significance test
  #   meta: list of Gamma and Omega (inverse Gamma) matrices
  
  Gamma <- crossprod(Zstat) / nrow(Zstat) # positive semi-definite
  Omega   <- chol2inv(chol(Gamma + 1e-8*diag(diag(Gamma)))) # small regularization to ensure invertibility
  Alpha <- Omega  %*% t(Zstat)                      # m x q
  qform <- colSums(t(Zstat) * Alpha) 
  pval = 1-pchisq(qform,df=ncol(Omega))
  Alpha_p <- sweep(Alpha, 1, sqrt(diag(Omega)), "/")
  
  ### COMPUTE SCORE OF MENDELIANISM
  SoMs = c()
  if (SoM){
    leads = find_towers(chr_pos,pval) # find towers with Q statistics
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
