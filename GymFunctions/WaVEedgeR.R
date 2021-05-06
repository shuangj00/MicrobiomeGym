# Description:

#' Apply the ZINB-WaVE version of edgeR for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import BiocManager
#' @import edgeR
#' @import SummarizedExperiment
#' @import zinbwave

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}. Must be a vector of two levels

# Outputs:
#' @return WaVEedegR.pval P-values for the taxa by WaVE-edegR
#' @return WaVEedegR.details Detailed output of WaVE-edgeR

WaVEedgeR = function(count.matrix, phenotype){
  # check input
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  # load libraries
  pkg.list = c("edgeR","SummarizedExperiment","zinbwave" )
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
  pheno.f = factor(paste0("grp_",phenotype))
  rownames(count.matrix) = paste0("grp_",phenotype)
  colData = DataFrame(Treatment = pheno.f)
  input.df = SummarizedExperiment(assays = list(counts = t(count.matrix)), colData=colData)
  
  # WaVE to reweight the zeros in the data
  message("Implementing ZINB-WaVE for edgeR...")
  zinbwave.default = zinbwave(input.df, K = 2, epsilon=1000, observationalWeights = TRUE)
  zinbwave.weights = assay(zinbwave.default, "weights")
  
  # apply edgeR
  dge = DGEList(assay(zinbwave.default))
  dge = edgeR::calcNormFactors(dge)
  
  design = model.matrix(~Treatment, data = colData(input.df))
  dge$weights = zinbwave.weights
  dge = estimateDisp(dge, design)
  glmfit = glmFit(dge, design)
  edgeR.lrt = glmWeightedF(glmfit, coef = ncol(glmfit$design))
  WaVEedegR.output = topTags(edgeR.lrt, n = "Inf")$table

  # organize edgeR output
  reorder.idx = base::match(colnames(count.matrix),rownames(WaVEedegR.output))
  WaVEedegR.output = WaVEedegR.output[reorder.idx,]
  WaVEedegR.pval = WaVEedegR.output$PValue
  WaVEedegR.pval[is.na(WaVEedegR.pval)] = 1
  
  return(list(WaVEedegR.pval = WaVEedegR.pval, 
              WaVEedegR.details = WaVEedegR.output))
}


