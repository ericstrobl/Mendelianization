Score_of_Mendelianism_tower <-function(Alpha,Sigma,Zstat,leads,chr_pos,alpha=5e-8){
  
  z_alpha = qnorm(1 - alpha/2) #z-statistic at p = alpha/2
  if (length(leads)==0) return(0)
  
  SoMs = c()
  d = ncol(as.matrix(Alpha))
  
  for (i in 1:length(leads)) {
    if (d > 1){
      a <- Alpha[, leads[i], drop = FALSE]      # alpha_j (m x 1)
    } else{
      a <- Alpha
    }
    denom_i <- sqrt(pmax(as.numeric(t(a) %*% Sigma %*% a), 1e-20))
    
    Z <- as.vector((Zstat %*% a) / denom_i)
    members = find_lead_tower(chr_pos,2*pnorm(-abs(Z)),leads[i])$members
    
    if (length(members>0)){
      Z = pmax(abs(Z)-z_alpha,0)
      SoM <- mean(Z[members]) / (mean(Z[members]) + sum(Z[-members]))
    } else{
      SoM = 0
    }
    
    SoMs = c(SoMs,SoM)
  }
  names(SoMs) = leads

  
  return(SoMs)
}