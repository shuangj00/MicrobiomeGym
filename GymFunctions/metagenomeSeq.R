# Description:

#' Apply edgeR for the differential abundance analysis of microbiome count data

# Dependencies:

#' @import BiocManager
#' @import metagenomeSeq

# Parameters:
#' @param count.matrix A sample-by-taxon original count matrix. All the entry should be nonnegative counts
#' @param phenotype The phenotype indicator vector for all samples in \code{count.matrix}

# Outputs:
#' @return metagenomeSeq.pval P-values for the taxa by metagenomeSeq
#' @return metagenomeSeq.details Detailed output of metagenomeSeq

metagenomeSeq = function(count.matrix, phenotype){
  # check input
  count.matrix = check.countmatrix(count.matrix)
  phenotype = check.phenotype(phenotype)
  
  # load libraries
  if(!require("metagenomeSeq")){
    if(!require(BiocManager)){
      install.packages("BiocManager")
    }else{
      BiocManager::install("metagenomeSeq")
      library("metagenomeSeq")
    }
  }
 
  # prepare the input for fitting metagenomeSeq
  counts.tab = t(count.matrix)
  rownames(counts.tab) = sprintf("OTU%03d", seq(nrow(counts.tab)))
  
  # convert to AnnotatedDataFrame
  ADF = AnnotatedDataFrame(data.frame(phenotype, row.names = colnames(counts.tab)))
  
  # define dummy 'feature' data for OTUs, using their name Helps with
  # extraction and relating to taxonomy later on.
  TDF = AnnotatedDataFrame(data.frame(OTUname = rownames(counts.tab), row.names =  rownames(counts.tab)))
  
  # create the metagenomeSeq object
  MGS = newMRexperiment(counts = counts.tab, phenoData = ADF, featureData = TDF)
  
  # trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
  MGS = cumNorm(MGS)
  
  # fit metagenomeSeq
  metagenomeSeq.fit = fitFeatureModel(MGS, model.matrix(~phenotype))
  metagenomeSeq.pval = metagenomeSeq.fit@pvalues
  metagenomeSeq.pval[is.na(metagenomeSeq.pval)] = 1
  metagenomeSeq.pval[is.infinite(metagenomeSeq.pval)] = 1
  
  return(list(metagenomeSeq.pval = metagenomeSeq.pval, 
              metagenomeSeq.details = metagenomeSeq.fit))
}


