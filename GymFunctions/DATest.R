# Description:

#' Apply DA test for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import devtools
#' @import DAtest
#' @import stringr

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}
#' @param core.num Number of cores to use for parallel computing. 
#' @param R Number of times to run the tests. Default is 20
#' @param seed Random seed

# Outputs:
#' @return best.method Name of the best method evaluated on the data by the DA test
#' @return best.pvalue p-values given by the best method
#' @return full.pvalue A matrix of the p-values given by all the methods used in the DA test

DATest = function(count.matrix, 
                  phenotype, 
                  core.num = 2, 
                  R = 20, 
                  seed = 123){
  # check input
  source("check.input.R")
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  seed = check.seed(seed)
  set.seed(seed)
  if(core.num > (detectCores() - 1)){
    stop("Number of cores used for parallel computing exceeds the maximum")
  } else if (core.num != floor(core.num)){
    stop("Number of cores used for parallel computing must be an integer")
  }
  if(R != floor(R)){
    stop("Number of times to run the tests must be an integer")
  }
  
  # load libraries
  if(!require(DAtest)){
    if (!requireNamespace("devtools")) {
      install.packages("devtools")
      devtools::install_github("Russel88/DAtest")
    }
  }
  if(!require(stringr)){
    install.packages("stringr", version = "1.4.0")
  }
  
  # perform DA test
  group.idx = as.factor(phenotype)
  rownames(count.matrix) = group.idx
  DAtest = testDA(t(count.matrix), predictor = group.idx, cores = core.num, verbose = FALSE)
  
  # summarize the DA test output and find the best test #
  DAsummary = summary(DAtest)
  best.idx = which.max(DAsummary$AUC)
  DAbest = DAsummary$Method[best.idx]
  DAbest.abbrv = str_extract_all(DAbest, "\\([^()]+\\)")[[1]]
  DAbest.abbrv = substring(DAbest.abbrv, 2, nchar(DAbest.abbrv)-1)

  # p-values (raw) from all the methods evaluated
  allDA.result = allDA(t(count.matrix), predictor = group.idx,  cores = core.num) 
  full.pval = allDA.result$raw
  full.pval[is.na(full.pval)] = 1
  
  # reorder the taxa
  taxa.idx = unlist(lapply(full.pval$Feature, as.character))
  idx.match = base::match(colnames(count.matrix),taxa.idx)
  idx.match = idx.match[!is.na(idx.match)]
  full.pval = full.pval[idx.match, ]
  rownames(idx.match)  = NULL
  best.pval = full.pval[, DAbest.abbrv]
  
  return(list(best.method = DAbest,
              best.pvalue = best.pval,
              full.pvalue = full.pval))
}


