# Description:

#' Apply edgeR for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import BiocManager
#' @import edgeR

# Parameters:

#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}

# Outputs:

#' @return edgeR.pval P-values given by edgeR
#' @return edgeR.details Detailed output from edgeR

edgeR = function(count.matrix, phenotype){
  # check input
  source("check.input.R")
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  # load libraries
  if(!require("edgeR")){
    if(!require(BiocManager)){
      install.packages("BiocManager")
    }else{
      BiocManager::install("edgeR")
      library("edgeR")
    }
  }
 
  # fit edgeR
  dge.l = DGEList(counts = t(count.matrix), group = phenotype)
  edgeR.sizefactor = edgeR:::calcNormFactors(dge.l, method = "RLE")
  
  # add pesudo count (1) to calculate size factor
  edgeR.sizefactor[[2]]$norm.factors = edgeR::calcNormFactors(t(count.matrix + 1), method = "RLE")
  if (!all(is.finite(edgeR.sizefactor$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  
  # perform exact test for edgeR
  est.nb = estimateCommonDisp(edgeR.sizefactor)
  est.nb = estimateTagwiseDisp(est.nb)
  edgeR.exactTest = exactTest(est.nb)
  edgeR.pval = edgeR.exactTest$table$PValue
  names(edgeR.pval) = rownames(edgeR.exactTest$table)
  
  return(list(edgeR.pval = edgeR.pval, 
              edgeR.details = edgeR.exactTest))
}


