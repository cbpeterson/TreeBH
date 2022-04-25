# Function: get_TreeBH_selections -------------------------------------
# Performs multi-level hierarchical selection 
#
# Args:
#   pvals:    Vector of N p-values (some of which may be NA) for most granular
#               level of tree
#   groups:   N x L matrix with groups[n, l] = k means that hypothesis n belongs
#               to the kth group in level l. Final column should just be 1:N
#               since most granular grouping is individual hypotheses
#   q:        Vector of L error rates to be targeted for each level in tree
#   test:     Vector of L-1 combination methods to obtain level l-1 p-values from
#               level l  p-values. Alternatives are "simes" and "fisher". Default is
#               "simes" at each level.
# Returns:
#   sel:      N x L binary matrix where sel[n, l] = 1 for first element of 
#               level l group indicates that group was selected
get_TreeBH_selections <- function(pvals, groups, q, test = "simes") {  
  # Total number of hypotheses at most granular level
  N <- length(pvals)
  
  # Total number of levels
  L <- ncol(groups)
  
  # Input checks
  if (length(q) != L) {
    stop("Must specify target q for each level of hierarchy")
  }
  if (min(q) < 0 || max(q) > 1) {
    stop("Target q must be in the range 0 - 1")
  }
  if (min(pvals, na.rm = TRUE) < 0 || max(pvals, na.rm = TRUE) > 1) {
    stop("P-values must be in the range 0 - 1")
  }
  if (length(pvals) != nrow(groups)) {
    stop("Dimension mismatch between pvals and groups")
  }
  
  # Check that groups are nested i.e. that all members in a group at any level
  # have one parent
  for (cur_level in 2:L) {
    cur_groups <- unique(sort(groups[ , cur_level]))
    for (cur_group in cur_groups) {
      cur_group_inds <- which(groups[, cur_level] == cur_group)
      parent_groups <- groups[cur_group_inds, cur_level - 1]
      if (length(unique(parent_groups)) != 1) {
        stop("Groups must be nested within hierarchy")
      }
    }
  }
  
  # Check that final level is most granular i.e. at level of individual hypotheses
  if (!length(unique(groups[ , L])) == nrow(groups)) {
    stop("Assumption is that lowest level in tree corresponds to individual hypotheses")
  }
  
  # Define selection matrix to be used in return value
  sel <- matrix(0, nrow = N, ncol = L)
  
  # Statistic used to aggregate p-values at each level
  if (length(test) == 1) {
    test <- rep(test, L - 1)
  }
  
  # Define matrix with adjusted target levels
  q_adj <- matrix(NA, nrow = N, ncol = L)
  
  # Set 1st col to unadjusted input value
  q_adj[ , 1] <- q[1]
  
  # Compute aggregated p-values for groups, starting from the bottom of the tree
  group_pvals <- matrix(NA, nrow = N, ncol = L)
  
  # Most granular level is pvals given as input
  group_pvals[ , L] <- pvals
  
  # Now compute aggregated p-values, moving from finer to coarser levels
  aggregate_levels <- (L-1):1
  for (cur_level in aggregate_levels) {
    cur_groups <- unique(sort(groups[ , cur_level]))
    for (cur_group in cur_groups) {
      cur_group_inds <- which(groups[, cur_level] == cur_group)
      cur_pvals <- group_pvals[cur_group_inds, cur_level + 1]
      # Set first p-value in group to aggegrated p-value to represent entire group
      if (test[cur_level] == "simes") {
        group_pvals[cur_group_inds[1], cur_level] <- get_simes_p(cur_pvals)
      } else if (test[cur_level] == "fisher") {
        group_pvals[cur_group_inds[1], cur_level] <- get_fisher_p(cur_pvals)
      } else {
        stop("Options for parameter test are 'simes' and 'fisher'")
      }
    }
  }
  
  # Now perform selection starting from the top level and moving down the tree
  for (cur_level in 1:L) {
    # Stop if there were no selections in previous level
    if (cur_level > 1 && sum(sel[ , cur_level - 1]) == 0) {
      break
    }
    
    # Get hypotheses selected in previous level
    if (cur_level == 1) {
      sel_prev <- 1 
    } else {
      sel_prev <- which(sel[ , cur_level - 1] == 1)
    }
    
    # Iterate over the selected hypotheses
    for (parent_sel in sel_prev) {
      # Identify selected children of current parent
      if (cur_level == 1) {
        parent_group_ind <- 1
        child_inds <- which(!is.na(group_pvals[ , cur_level]))
      } else {
        parent_group_num <- groups[parent_sel, cur_level - 1]
        parent_group_ind <- min(which(groups[ , cur_level - 1] == parent_group_num))
        
        # Take as representative index the first index within each group
        # This is the same as the entries with non-NA p-values given how they 
        # were calculated above
        child_inds <- which(groups[ , cur_level - 1] == parent_group_num &
                              !is.na(group_pvals[ , cur_level]))
      }
      if (length(child_inds) > 1) {
        sel_ind_within_group <- which(qvalue(group_pvals[child_inds, cur_level],
                                             lambda = 0)$qvalue <=
                                        q_adj[parent_group_ind, cur_level])
      } else {
        sel_ind_within_group <- which(group_pvals[child_inds, cur_level] <=
                                        q_adj[parent_group_ind, cur_level])
      }
      sel[child_inds[sel_ind_within_group], cur_level] <- 1
      
      # Compute adjusted target for families of selected children as long as we
      # are not at bottom level
      if (cur_level < L) {
        # Number of rejections for current parent selection
        R_parent_sel <- length(sel_ind_within_group)
        
        # Number of elements in current hypothesis family i.e. for current parent selection
        n_parent_sel <- length(child_inds)
        
        q_adj[child_inds[sel_ind_within_group], cur_level + 1] <-
          q[cur_level + 1] *
          q_adj[parent_sel, cur_level] / q[cur_level] *
          R_parent_sel / n_parent_sel
      }
    }
  }
  
  # Return the matrix of selections at each level
  sel
}
