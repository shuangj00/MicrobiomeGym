# Description:

#' Apply corncob (Count Regression for Correlated Observations with the Beta-binomial) for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import LDM

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}
 

# Outputs:
#' @return ldm.pval P-values given by the corncob
#' @return ldm.details Full output from LDM

LDMtest = function(count.matrix, 
                   phenotype){
  # check input
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  
  # load libraries
  if(!require(LDM)){
    message("Download the R package LDM from http://web1.sph.emory.edu/users/yhu30/software.html")
  }
  
  # get the taxa names
  if(is.null(colnames(count.matrix))){
    taxa.names = paste0("Taxon", seq(1, ncol(count.matrix)))
  }else{
    taxa.names = colnames(count.matrix)
  }
  
  # need to remove the taxa with zeros across all the samples
  rm.idx = apply(count.matrix, 2, function(x){all(x == 0)})
  rm.taxa = taxa.names[rm.idx]
  keep.taxa = taxa.names[!rm.idx]
  
  # apply LDM test
  count.matrix = count.matrix[, !rm.idx]
  otu_tab <<- t(count.matrix)
  group.idx <<- as.character(phenotype)
  fit.ldm = ldm(otu_tab ~ (group.idx) ,
             test.otu=TRUE, 
             test.global=T,
             dist.method="bray",
             n.rej.stop=100)
  
  # remove the temporal variable otu_tab and group.idx
  objs = ls(pos = ".GlobalEnv")
  rm(list = objs[grep("otu_tab", objs)], pos = ".GlobalEnv")
  rm(list = objs[grep("group.idx", objs)], pos = ".GlobalEnv")
  
  # summarize the output 
  pval = as.numeric(fit.ldm$p.otu.freq)
  names(pval) = keep.taxa
  
  return(list(ldm.pval = pval,
              ldm.details = fit.ldm))
}


# Description:

#' Solve for the maximum likelihood estimator for individual taxon using corncob (Count Regression for Correlated Observations with the Beta-binomial) for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import devtools
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

