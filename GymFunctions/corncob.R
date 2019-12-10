# Description:

#' Apply corncob (Count Regression for Correlated Observations with the Beta-binomial) for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import devtools
#' @import corncob
#' @import plyr

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}
#' @param count.min An integer greater than or equal to 1, which is the minimum number of nonzero counts required
#'        for each phenotype group in implementing corncob. Taxa with observed count less than \code{count.min} in 
#'        at least one phenotype group would be filtered. 

# Outputs:
#' @return corncob.pval P-values given by the corncob
#' @return removed.taxa Name of the filtered taxa

corncob = function(count.matrix, 
                   phenotype, 
                   count.min = 1){
  # check input
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  
  # load libraries
  if(!require(corncob)){
    if (!requireNamespace("devtools")) {
      install.packages("devtools")
      devtools::install_github("bryandmartin/corncob")
    }
  }
  if(!require(plyr)){
    install.packages("plyr", version = "1.8.4")
  }
  
  # get the taxa names
  if(is.null(colnames(count.matrix))){
    taxa.names = paste0("Taxon", seq(1, ncol(count.matrix)))
  }else{
    taxa.names = colnames(count.matrix)
  }
  
  # filter the sparse taxa by group (required for fitting corncob)
  nonzero.counts = apply(count.matrix, 2, function(x){dft = data.frame(count = x, pheno = phenotype);
  nonzero.count = ddply(dft, c("pheno"), summarise,N = sum(count != 0)); 
  nonzero.count = nonzero.count[, 2]; 
  nonzero.count})
  
  rm.idx = apply(nonzero.counts,2, function(x){ifelse(any(x < count.min), T, F)} )
  count.matrix = count.matrix[, !rm.idx ]
  rm.taxa = taxa.names[rm.idx]
  
  # apply corncob
  group.idx = as.factor(phenotype)
  pval.list = sapply(seq(1,ncol(count.matrix) ), single.lrt, count.matrix = count.matrix, phenotype = group.idx)
  names(pval.list) = taxa.names[!rm.idx]
  
  return(list(corncob.pval = pval.list,
              removed.taxa = rm.taxa))
}


# Description:

#' Solve for the maximum likelihood estimator for individual taxon using corncob (Count Regression for Correlated Observations with the Beta-binomial) for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import corncob

# Parameters:

#' @param idx Index (column) of the taxon in \code{count.matrix}
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}

# Outputs:
#' @return p-values of the likelihood ratio test between the null model and the model with phenotype as a explanatory variable

single.lrt = function(idx,count.matrix, phenotype){
  require(corncob)
  taxon.i = count.matrix[, idx]
  tmp.df = data.frame("W" = taxon.i, "M" = rowSums(count.matrix), group.idx = phenotype)
  mod1 = bbdml(formula = cbind(W, M - W) ~ group.idx,
               phi.formula = ~ group.idx,
               data = tmp.df)
  mod2 = bbdml(formula = cbind(W, M - W) ~ 1,
                phi.formula = ~ 1,
                data = tmp.df)
  
  pval = lrtest(mod1, mod2)
  return(pval)
}

