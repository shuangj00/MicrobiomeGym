# Description:

#' Apply Wilcoxon rank-sum test for the differential abundance analysis of microbiome count data

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}. Should only have 2 groups.

# Outputs:
#' @return A list of detailed wilcoxon rank-sum test output and the p-values

WilcoxonRankSum = function(count.matrix, phenotype){
  # check the phenotype 
  pheno.level = levels(factor(phenotype))
  if(length(pheno.level) >= 3){
    stop("Wilcoxon rank sum can handle only two phenotype levels")
  }
  
  # covert to compositional data
  RA.matrix = apply(count.matrix, 1, function(x) { 
    s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
  
  # perform the test 
  RA.matrix = t(RA.matrix)
  WC.details = apply(RA.matrix, 2, function(x){wilcox.test(x[phenotype == pheno.level[1]], x[phenotype == pheno.level[2]], exact = FALSE)})
  WC.pval = unlist(lapply(WC.details, function(x){x$p.value} ))
  
  return(list(Wilcoxon.pval = WC.pval, 
              Wilcoxon.details = WC.details))
}


