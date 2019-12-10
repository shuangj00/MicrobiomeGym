# Description:

#' Apply Kruskal Wallis test for the differential abundance analysis of microbiome count data

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}. Can have multiple groups

# Outputs:

#' @return KW.pval P-values given by Kruskal-Wallis test
#' @return KW.details Detailed output from Kruskal-Wallis test

KruskalWallis = function(count.matrix, phenotype){
  # covert to compositional data
  RA.matrix = apply(count.matrix, 1, function(x) { 
    s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
  RA.matrix = t(RA.matrix)
  
  # apply KW test 
  KW.details = apply(RA.matrix, 2, function(x){kruskal.test(x, phenotype)})
  KW.pval = unlist(lapply(KW.details, function(x){x$p.value} ))
  
  return(list(KW.pval = KW.pval, 
              KW.details = KW.details))
}


