# Description:

#' Apply DESeq2 for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import BiocManager
#' @import DESeq2

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}

# Outputs:

#' @return DESeq2.pval P-values given by DESeq2
#' @return DESeq2.details Detailed output from DESeq2


DESeq2 = function(count.matrix, phenotype){
  # check input
  source("check.input.R")
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  # load library for DESeq2
  if(!require(DESeq2)){
    if(!require(BiocManager)){
      install.packages("BiocManager")
    }
    BiocManager::install("DESeq2")
    library(DESeq2)
  }
 
  # add pseudo count 1
  count.matrix = t(count.matrix + 1)
  
  # fit DESeq2 
  DDS.res = DESeqDataSetFromMatrix(count.matrix , DataFrame(phenotype), ~ phenotype)
  DESeq2.fit = DESeq(DDS.res)
  DESeq2.res = results(DESeq2.fit)
  DESeq2.pval = DESeq2.res[, "pvalue"]
  names(DESeq2.pval) = rownames(DESeq2.res)
  
  return(list(DESeq2.pval = DESeq2.pval, 
              DESeq2.details = DESeq2.res))
}


