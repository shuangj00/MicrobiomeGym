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

WaVEedgeR = function(count.matrix, phenotype, count.min = 1){
  # check input
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  # load libraries
  pkg.list = c("edgeR","SummarizedExperiment","zinbwave" )
  for(pkg in pkg.list){
    if(!require(package = pkg, warn.conflicts = FALSE,character.only = TRUE, quietly = TRUE)){
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
  
  pval.list = WaVEedegR.pval
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
  
  return(list(WaVEedegR.pval = pval.ret, 
              removed.idx = removed.idx,
              removed.taxa = rm.taxa,
              WaVEedegR.details = WaVEedegR.output))
}


