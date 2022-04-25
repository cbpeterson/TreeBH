# Function: get_simes_p -------------------------------------------------------
# Compute Simes p-value given a list of p-values, assuming that NAs do not count
# toward the total number of tests
# 
# Args:
#   pvals: vector of p-values (which may include NAs)
# 
# Returns:
#   value of the Simes p-value
get_simes_p <- function(pvals) {
  # Don't include NAs
  pvals <- pvals[!is.na(pvals)]
  
  if (length(pvals) == 0) {
    # Return NA if input had length 0 or was all NAs
    NA
  } else {
    # Else compute Simes p
    pvals <- sort(as.numeric(pvals))
    length(pvals) * min(pvals / seq(1:length(pvals)))
  }
}

# Function: get_fisher_p ------------------------------------------------------
# Compute Fisher p-value given a list of p-values, assuming that NAs do not count
# toward the total number of tests
# 
# Args:
#   pvals: vector of p-values (which may include NAs)
# 
# Returns:
#   value of the Fisher p-value
get_fisher_p <- function(pvals) {
  # Don't include NAs
  pvals <- pvals[!is.na(pvals)]
  
  if (length(pvals) == 0) {
    # Return NA if input had length 0 or was all NAs
    NA
  } else {
    # Else compute Fisher p-value
    chisq_df <- 2 * length(pvals)
    pchisq(-2 * sum(log(pvals)), chisq_df, lower.tail = FALSE)
  }
}


# Function: select_pvals ------------------------------------------------------
# Applies BH procedure to vector of p-values, assuming that NAs do not count 
# toward to total number of tests
#
# Args:
#   pvals: vector of p-values for current group/family of hypotheses
#   q_adj: adjusted target for BH procedure
#
# Returns:
#   Binary indicator vector of whether p-value passes selection threshold
select_pvals <- function(pvals, q_adj) {
  is_sel <- rep(NA, length(pvals))
  is_sel[which(is.na(pvals))] <- FALSE
  if (sum(!is.na(pvals)) > 1) {
    is_sel[which(is.na(is_sel))] <- qvalue(pvals[which(is.na(is_sel))], lambda = 0,
                                           fdr.level = q_adj)$significant
  } else {
    is_sel[which(is.na(is_sel))] <- pvals[which(!is.na(pvals))] <= q_adj
  }
  is_sel
}