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

# Outputs:
#' @return WaVEDESeq2.pval P-values given by WaVEDESeq2
#' @return WaVEDESeq2.details Detailed output of WaVEDESeq2

WaVEDESeq2 = function(count.matrix, phenotype){
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
  
  return(list(WaVEDESeq2.pval = WaVEDESeq2.res$pvalue, 
              WaVEDESeq2.details = WaVEDESeq2.res))
}


