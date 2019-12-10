# Description:

#' Apply two-sample t-test for the differential abundance analysis of microbiome count data

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}. Should only have 2 groups.

# Outputs:

#' @return t.pval P-values given by the two-sample t test
#' @return t.details Detailed output from the two-sample t test

tTest = function(count.matrix, phenotype){
  # check the phenotype 
  pheno.level = levels(factor(phenotype))
  if(length(pheno.level) >= 3){
    stop("Two sample t-test can handle only two phenotype levels")
  }
  
  # covert to the compositional data 
  RA.matrix = apply(count.matrix, 1, function(x) { 
    s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
  
  # perform two-sample t test 
  RA.matrix = t(RA.matrix)
  t.details = apply(RA.matrix, 2, function(x){t.test(x[phenotype == pheno.level[1]], x[phenotype == pheno.level[2]])})
  t.pval = unlist(lapply(t.details, function(x){x$p.value} ))
  
  return(list(t.pval = t.pval, 
              t.details = t.details))
}


