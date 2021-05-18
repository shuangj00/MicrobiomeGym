# Description:

#' Apply the ZINB-WaVE version of DESeq2 for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import BiocManager
#' @import DESeq2
#' @import SummarizedExperiment
#' @import zinbwave

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}
#' @param count.min An integer greater than or equal to 1, which is the minimum number of nonzero counts required
#'        for each phenotype group in implementing corncob. Taxa with observed count less than \code{count.min} in 
#'        at least one phenotype group would be filtered. 

# Outputs:
#' @return WaVEDESeq2.pval P-values given by WaVEDESeq2
#' @return WaVEDESeq2.details Detailed output of WaVEDESeq2

WaVEDESeq2 = function(count.matrix, phenotype,count.min = 1){
  # check input
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  # load libraries
  pkg.list = c("DESeq2","SummarizedExperiment","zinbwave" )
  for(pkg in pkg.list){
    if(!require(package = pkg, character.only = TRUE)){
      if(!require(BiocManager)){
        install.packages("BiocManager")
      }else{
        BiocManager::install(pkgs = pkg)
        library(pkg, character.only = TRUE)
      }
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
  removed.idx = which(rm.idx)
  
  # convert the raw count table to a "SummarizedExperiment" object
  pheno_idx = paste0("grp_",phenotype)
  pheno.f = factor(pheno_idx)
  rownames(count.matrix) = pheno_idx
  colData = DataFrame(Treatment = pheno.f)
  input.df = SummarizedExperiment(assays = list(counts = t(count.matrix)), colData=colData)
  
  # WaVE to reweight the zeros in the data
  message("Implementing ZINB-WaVE for DESeq2...")
  zinbwave.default = zinbwave(input.df, K = 2, epsilon=1000)
  
  # apply DESeq2
  dds = DESeqDataSet(zinbwave.default, design = ~ Treatment)
  dds = DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
  WaVEDESeq2.res = results(dds)
  
  pval.list = WaVEDESeq2.res$pvalue
  pval.ret = pval.list
  if(sum(rm.idx)!=0){
    for(idx in removed.idx){
      if(idx == 1){
        pval.ret = c(NA, pval.ret)
      }else if(idx == length(taxa.names)){
        pval.ret = c(pval.ret, NA)
      }else{
        pval.ret = c(pval.ret[1:(idx - 1)], NA, pval.ret[idx:length(pval.ret)])
      }
      
    }
  }
  names(pval.ret) = taxa.names
  return(list(WaVEDESeq2.pval = pval.ret, 
              removed.idx = removed.idx,
              removed.taxa = rm.taxa,
              WaVEDESeq2.details = WaVEDESeq2.res))
}


