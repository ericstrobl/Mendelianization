find_towers <- function(chr_pos, pvals, p_thr = 5e-8, kb = 500) {
  # find lead variants based on Q statistics with gap-and-pad approach
  
  # significant indices (in the working order)
  sig <- which(pvals <= p_thr)
  if (!length(sig)) return(list(leads = integer(0), members = list()))
  
  leads <- integer(0)
  members <- list()
  
  # work chromosome-by-chromosome using only significant variants
  for (ch in unique(chr_pos$chr[sig])) {
    idx <- sig[chr_pos$chr[sig] == ch]              # row indices (working order) of sig SNPs on this chr
    if (!length(idx)) next
    pos <- chr_pos$pos[idx]                         # their positions (bp), already sorted if df is sorted
    
    # initial blocks: break when gap > kb
    new_block <- c(TRUE, diff(pos) > kb * 1000)
    block_id  <- cumsum(new_block)
    
    # expand each block by ±kb/2 and merge overlapping intervals, to safely include all nearby variants
    starts <- tapply(pos, block_id, min) - kb * 1000
    ends   <- tapply(pos, block_id, max) + kb * 1000
    o <- order(starts); starts <- starts[o]; ends <- ends[o]
    
    mst <- numeric(0); men <- numeric(0)
    for (i in seq_along(starts)) {
      if (!length(mst) || starts[i] > men[length(men)]) {
        mst <- c(mst, starts[i]); men <- c(men, ends[i])
      } else {
        men[length(men)] <- max(men[length(men)], ends[i])
      }
    }
    
    # assign members/leads from merged intervals (still only significant SNPs)
    for (k in seq_along(mst)) {
      # pick significant indices in this chr whose positions fall in the merged interval
      in_int_sig <- idx[pos >= mst[k] & pos <= men[k]]
      if (!length(in_int_sig)) next
      members[[length(members) + 1L]] <- in_int_sig
      leads <- c(leads, in_int_sig[which.min(pvals[in_int_sig])])
    }
  }
  
  return(leads)
}
