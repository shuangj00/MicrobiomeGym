# Description:

#' Apply ANOVA for the differential abundance analysis of microbiome count data

# Parameters:

#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}. Can have multiple levels

# Outputs:

#' @return ANOVA.pval P-values given by ANOVA
#' @return ANOVA.details Detailed output from ANOVA

ANOVA = function(count.matrix, phenotype){
  # covert to compositional data 
  RA.matrix = apply(count.matrix, 1, function(x) { 
    s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
  RA.matrix = t(RA.matrix)
  # apply ANOVA
  ANOVA.details = apply(RA.matrix, 2, function(x){aov(x~ phenotype)})
  ANOVA.pval = unlist(lapply(ANOVA.details, function(x){summary(x)[[1]]$`Pr(>F)`[1]} ))
  
  return(list(ANOVA.pval = ANOVA.pval, 
              ANOVA.details = ANOVA.details))
}


