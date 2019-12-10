# Description
#' run MCMC algorithm to fit the ZINB-DPP model for the input n-by-p count matrix

# Dependencies:
#' @import Rcpp

# Parameters:
#' @param count.matrix A count matrix from metagenomic shotgun sequencing or 16SrRNA sequencing technologies. 
#'                  Columns represent the taxonomic features and rows represents the samples.
#' @param phenotype A phenotype indicator vector for all the samples (should be a vector of 0s if only have
#'                  one phenotype, or a vector of 0 and 1 if have 2 phenotypes. The current method can handle
#'                  at most 2 phenotypes)
#' @param N.mcmc Number of MCMC iterations
#' @param b Shape hyper-parameter for the variance term
#' @param h Scale hyper-parameter for the variance term
#' @param count.min: minimum number of nonzero counts (for each phenotype group) required for a taxon to fit the ZINB-DPP model
#' @param seed: Random seed

# Output: 
#' @return fold.change Fold change (group 1 over group 0) of the normalized abundance for all p taxa
#' @return remove.idx A binary vector of length p that indicates if a taxon is dropped from fitting the ZINB-DPP model
#' @return size.factor A matrix that stores the MCMC output (after burn-in the first half) for the size factor of all the sample
#' @return alpha.matrix A matrix of the posterior mean of {alpha_ij}, i = 1, ..., n, and j = 1, ..., p
#' @return phi A vector of the posterior mean of the dispersion parameter
#' @return H.ppi An n-by-p matrix of the posterior probability of inclusion (PPI) of being a missing value 
#'         for each count in the input matrix
#' @return gamma.ppi A vector of the posterior probability of inclusion (PPI) for each of p taxon to be discriminating between 
#'         patient phenotypes (if we have 2 groups)
#' @return gamma.accept.rate Acceptance rate for updating gamma in the Metropolis–Hastings algorithm
#' @return phi.accept.rate Acceptance rate for updating phi in the Metropolis–Hastings algorithm
#' @return alpha.accept.rate Acceptance rate for updating alpha in the Metropolis–Hastings algorithm

ZINBDPP = function(count.matrix = count.matrix,
                   phenotype = phenotype, 
                   N.mcmc = 10000, 
                   b = 2, 
                   h = 20, 
                   count.min = 2,
                   seed = 123){
  # check input 
  source("check.input.R")
  if(! all(count.matrix == floor(count.matrix)) | any(count.matrix < 0)){
    stop("Elements in the input count matrix must be nonnegative integers")
  }
  if(!is.matrix(count.matrix)){
    count.matrix = as.matrix(count.matrix) 
  }
  phenotype = check.phenotype(phenotype)
  seed = check.seed(seed)
  if(b <= 0 | h <= 0){
    stop("parameter 'b' and 'h' must be strictly positive")
  }
  if(N.mcmc != floor(N.mcmc) ){
    stop("N.mcmc should be a large positive integer to ensure convergence (suggested value: >= 10000)")
  }
  if(count.min != floor(count.min)){
    stop("count.min must be an integer")
  }
  if(count.min <= 1 ){
    stop("count.min must be larger than 1 to fit ZINB-DPP model" )
  }
  if(any(count.min >= table(phenotype))){
    stop("count.min should be smaller than any group size")
  }
  
  # load package
  if (!require(Rcpp)) {install.packages("Rcpp", version = "1.0.2")}
  sourceCpp(file = "fit.ZINBDPP.cpp")
  
  
  # set initial values for the MCMC
  s0 = rep(1, nrow(count.matrix))
  S0 = matrix(1, 1, 1)
  
  
  # fit ZINB-DPP model 
  set.seed(seed)
  message(paste0("Start fitting the ZINB-DPP model with ", N.mcmc, " iterations. This may take a while... \n"))
  MCMC.output  = fitZINBDPP(Y = count.matrix,
                            z = phenotype, 
                            s = s0,
                            iter = N.mcmc, 
                            DPP = TRUE,
                            S = S0, 
                            aggregate = FALSE, 
                            b = b, 
                            h = h, 
                            mincount = count.min,
                            MRF = FALSE,
                            G = S0)
  MCMC.temp = MCMC.output
  rm(MCMC.output)
  MCMC.temp$fold.change = apply(MCMC.temp$fold.change, 2, mean)
  MCMC.temp$size.factor = apply(MCMC.temp$size.factor, 2, mean)
  MCMC.temp$pi = NULL
  MCMC.output = MCMC.temp
  return(MCMC.output)
}
