plot_tower <- function(out, Z, lead, chr_pos){
  # out = output of Mendelianization
  # lead = index of lead variant you want to focus on
  # chr_pos: data frame with one column containing chromosomes and another containing sorted positions (increasing)
  
  a <- out$Alpha[, lead, drop = FALSE]
  denom_i <- sqrt(pmax(as.numeric(t(a) %*% out$Gamma %*% a), 1e-20))
  
  Zn <- as.vector((Z %*% a) / denom_i)
  ps= 2*pnorm(-abs(Zn))
  GWAS <- data.frame(
    CHR = as.numeric(chr_pos[[1]]),
    BP  = as.numeric(chr_pos[[2]]),
    P   = ps,
    SNP = as.numeric(chr_pos[[1]]),
    stringsAsFactors = FALSE
  )
  require(qqman)
  
  manhattan(
    GWAS,
    chr = "CHR", bp = "BP", p = "P",
    col = c("#4C72B0", "#55A868"),
    suggestiveline = FALSE,          # <- remove blue line
    genomewideline = -log10(5e-8),   # keep red line (or set FALSE to remove)
    cex = 0.5, cex.axis = 0.9, las = 2
  )
  

}

