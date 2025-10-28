find_lead_tower <- function(chr_pos, pvals, lead_var, p_thr = 5e-8, kb = 500) {
  ## find tower based on \mathcal{Z}-statistics with gap-and-pad approach

  kb_bp <- kb * 1000L
  chr0  <- chr_pos$chr[lead_var]
  pos0  <- chr_pos$pos[lead_var]
  
  # significant SNPs on the lead chromosome
  sig <- which(chr_pos$chr == chr0 & pvals <= p_thr)
  if (!length(sig)) return(list(lead = lead_var, members = integer(0)))
  
  # order by position and make kb-blocks
  o      <- order(chr_pos$pos[sig]); sig <- sig[o]; pos <- chr_pos$pos[sig]
  block  <- cumsum(c(TRUE, diff(pos) > kb_bp))
  
  # pick the block: if lead is sig, use its block; else the block whose padded span covers lead
  if (lead_var %in% sig) {
    bid <- block[match(lead_var, sig)]
  
    # expand the chosen block by ±kb and collect significant members inside that window
    block_sig <- sig[block == bid]
    w_lo <- min(chr_pos$pos[block_sig]) - kb_bp
    w_hi <- max(chr_pos$pos[block_sig]) + kb_bp
    members <- which(chr_pos$chr == chr0 & chr_pos$pos >= w_lo & chr_pos$pos <= w_hi & pvals <= p_thr)
    
    # pick lead: prefer the given lead if it's in members; else min-p within members
    lead <- if (lead_var %in% members) lead_var else members[which.min(pvals[members])]
  } else{
    lead = NULL
    members = NULL
    w_lo = NULL; w_hi = NULL
  }
  
  list(lead = lead, members = members, window = c(start = w_lo, end = w_hi))
}
