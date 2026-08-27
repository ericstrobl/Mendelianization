# Greedily remove redundant outcomes until the retained panel has an
# acceptable condition number or reaches the minimum permitted size.
select_outcomes <- function(Gamma,max_condition = 100,min_outcomes = 2) {
  # Use the supplied outcome names; create generic labels if none are present.
  outcomes <- colnames(Gamma)
  
  if (is.null(outcomes)) outcomes <- paste0("Y", seq_len(ncol(Gamma)))
  
  # Calculate the condition number for a proposed subset of outcomes.
  condition <- function(index) {
    # Convert the selected covariance submatrix to a correlation matrix so
    # that the condition number reflects dependence rather than marginal scale.
    R <- cov2cor(Gamma[index, index, drop = FALSE])
    
    # For a symmetric positive-definite matrix, the condition number is the
    # ratio of its largest to smallest eigenvalue.
    ev <- eigen(R,symmetric = TRUE,only.values = TRUE)$values
    
    # Treat a singular or non-positive-definite subset as infinitely ill-conditioned.
    if (min(ev) <= 0) return(Inf)
    
    max(ev) / min(ev)
  }
  
  # Begin with every outcome and record removals in their elimination order.
  keep <- seq_len(ncol(Gamma))
  removed <- character()
  
  # Continue only while the panel is too ill-conditioned and can still be reduced.
  while (
    length(keep) > min_outcomes &&
    condition(keep) > max_condition
  ) {
    # Evaluate the condition number obtained by removing each retained outcome.
    candidate_condition <- vapply(seq_along(keep), function(j) condition(keep[-j]), numeric(1))
    
    # Greedily remove the outcome whose deletion gives the smallest condition number.
    remove <- which.min(candidate_condition)
    
    # Report the selected outcome and the resulting improvement in conditioning.
    message(
      "Removing ", outcomes[keep[remove]],
      ": condition ", round(condition(keep), 1),
      " -> ", round(candidate_condition[remove], 1)
    )
    
    # Save the removed outcome's name before deleting its index from the panel.
    removed <- c( removed, outcomes[keep[remove]])
    
    keep <- keep[-remove]
  }
  
  # Confirm when screening was requested but the complete panel already met
  # the condition-number requirement.
  if (length(removed) == 0L && is.finite(max_condition)) {
    message(
      "No outcomes removed: condition number ",
      round(condition(keep), 1)
    )
  }
  
  # Return the retained indices and names, removal order, and final condition number.
  list(
    index = keep,
    outcomes = outcomes[keep],
    removed = removed,
    condition_number = condition(keep)
  )
}