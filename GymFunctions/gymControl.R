###note: this function has been moved to MicrobiomeGym.R#######################################################################################
#  Store the details of model input for \code{MicrobiomeGym}
#  maintaniner: Shuang Jiang, <shuangj@smu.edu>
##########################################################################################

# Description
#' Performing differential abundancn analysis for microbiome data using models provided
#' by Microbiome Gym 

# Dependencies:
#' @import parallel 

# Parameters:
#' @param method Name of the model provided by Microbiome Gym.
#'               Must be one of the following:
#'               "ANOVA", "tTest", "KruskalWallis", "WilcoxonRankSum",
#'               "DESeq2", "edgeR", "WaVEDESeq2", "WaVEedgeR",
#'               "metagenomeSeq", "DATest", "corncob", "LDMtest","ZINBDPP".
#' @param core.num An integer being greater than or equal 1. Number of cores to use for parallel computing.
#' @param R A positive integer. Number of times to run the tests. 
#' @param N.mcmc A positive integer. Number of MCMC iterations. Suggested value is at least 10,000.
#' @param b Shape hyper-parameter for the variance term for fitting the \code{ZINB-DPP} model.
#' @param h Scale hyper-parameter for the variance term for fitting the \code{ZINB-DPP} model.
#' @param count.min minimum number of none zero counts required for a taxon to be considered in the analysis 
#'                   (only needed for fitting the \code{ZINB-DPP} model or the \code{corncob} model)
#' @param seed A positive integer. Random seed.

# Output: 
#' @return An echo of the parameters specified

gymControl = function(method = "ANOVA",
                      core.num = 2, 
                      R = 20, 
                      N.mcmc = 10000,
                      b = 2, 
                      h = 20,
                      count.min = 2,
                      seed = 123){
  # check the model input 
  if (is.null(method)) {
    stop("Null method name not allowed")
  } else if(length(method) > 1){
    valid.idx = which(method %in% c("ANOVA", "tTest", "KruskalWallis", "WilcoxonRankSum",
                                    "DESeq2", "edgeR", "WaVEDESeq2", "WaVEedgeR",
                                    "metagenomeSeq", "DATest", "corncob", "LDMtest",
                                    "ZINBDPP"))
    valid.name = method[valid.idx]
    if (length(valid.name) == 0) {
      stop("Please provide a valid test name")
    }else if(length(valid.name) == 1){
      method = valid.name
    } else {
      method = valid.name[1]
      message(paste0("Can only implement one method at a time. Now running: ", method) )
    }
  } else if(! (method %in% c("ANOVA", "tTest", "KruskalWallis", "WilcoxonRankSum",
                             "DESeq2", "edgeR", "WaVEDESeq2", "WaVEedgeR",
                             "metagenomeSeq", "DATest", "corncob", "LDMtest",
                             "ZINBDPP") )){
    stop("Please provide a valid test name")
  }
 
  if(method == "DATest"){
    if (!require(parallel)) {install.packages("parallel", version = "3.6.0")}
    core.max = detectCores()
    if(core.num <= 0 | core.num > core.max | core.num!= floor(core.num) | is.null(core.num)){
      stop(paste0("Number of cores used for DAtest should be any integer between 1 - ",core.max ))
    }
    if(R != floor(R) | is.null(R)){
      stop("Number of times to run the DA test must be an integer")
    }
  }
  
  if(method == "ZINBDPP"){
    if(N.mcmc != floor(N.mcmc) ){
      stop("N.mcmc should be a large positive integer to ensure convergence (suggested value: >= 10000)")
    }
    if(b <= 0 | h <= 0){
      stop("parameter 'b' and 'h' must be strictly positive")
    }
    if(count.min <= 1 | count.min != floor(count.min)){
      stop("count.min must be an integer larger than 1 to fit ZINB-DPP model" )
    }
  }
  
  if(method == "corncob"){
    if(count.min < 1 | count.min != floor(count.min)){
      stop("count.min must be an integer greater than or equal to 1" )
    }
  }
    
  source("check.input.R")
  if(is.null(seed)){
    stop("Please provide a random seed")
  }else{
    seed = check.seed(seed)
  }
  
  # return the list for model input
  list(method = method, 
       core.num = core.num, 
       R = R, 
       N.mcmc = N.mcmc,
       b = b,
       h = h,
       count.min = count.min,
       seed = seed)
}
