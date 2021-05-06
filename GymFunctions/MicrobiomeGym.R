##########################################################################################
#  Microbiome Gym
#  maintaniner: Shuang Jiang, <shuangj@smu.edu>
##########################################################################################

# Description
#' Performing differential abundancn analysis for microbiome data using models provided
#' by Microbiome Gym 

# Dependencies:
#' @import (all the dependency packages are included in the method functions)

# Parameters:
#' @param count.matrix A count matrix from metagenomic shotgun sequencing or 16SrRNA sequencing technologies. 
#'                  Columns represent the taxonomic features and rows represents the samples.
#' @param phenotype A phenotype indicator vector for all the samples (should be a vector of 0s if only have
#'                  one phenotype, or a vector of 0 and 1 if have 2 phenotypes. The current method can handle
#'                  at most 2 phenotypes)
#' @param gymControl A list of values that define how this function
#'                   acts. See \code{\link{gymControl}}

# Output: 
#' @return model.name Name of the model used in the differential abundance analysis
#' @return model.details Detailed output of the given model
#' @return count.matrix The input count matrix
#' @return phenotype The input phenotype indicator vector

MicrobiomeGym = function(count.matrix,
                         phenotype,
                         gymControl = gymControl(),
                         includeData = FALSE) {
  # check the model
  method = gymControl$method
  message("Running ", method, " on the input data \n")
  if (method == "ANOVA") {
    source("../GymFunctions/ANOVA.R")
    gym.output = ANOVA(count.matrix = count.matrix, 
                     phenotype = phenotype)
  } else if (method == "tTest") {
    source("../GymFunctions/tTest.R")
    gym.output = tTest(count.matrix = count.matrix, 
                         phenotype = phenotype)
  } else if (method == "KruskalWallis"){
    source("../GymFunctions/KruskalWallis.R")
    gym.output = KruskalWallis(count.matrix = count.matrix, 
                               phenotype = phenotype)
  } else if (method == "WilcoxonRankSum"){
    source("../GymFunctions/WilcoxonRankSum.R")
    gym.output = WilcoxonRankSum(count.matrix = count.matrix, 
                                 phenotype = phenotype)
  } else if (method == "DESeq2"){
    source("../GymFunctions/DESeq2.R")
    gym.output = DESeq2(count.matrix = count.matrix, 
                        phenotype = phenotype)
  } else if (method == "edgeR"){
    source("../GymFunctions/edgeR.R")
    gym.output = edgeR(count.matrix = count.matrix, 
                    phenotype = phenotype)
  } else if (method == "WaVEDESeq2"){
    source("../GymFunctions/WaVEDESeq2.R")
    gym.output = WaVEDESeq2(count.matrix = count.matrix, 
                            phenotype = phenotype)
  } else if (method == "WaVEedgeR"){
    source("../GymFunctions/WaVEedgeR.R")
    gym.output = WaVEedgeR(count.matrix = count.matrix, 
                           phenotype = phenotype)
  } else if (method == "metagenomeSeq"){
    source("../GymFunctions/metagenomeSeq.R")
    gym.output = metagenomeSeq(count.matrix = count.matrix, 
                               phenotype = phenotype)
  } else if (method == "DATest"){
    source("../GymFunctions/DATest.R")
    message(paste0("using ", gymControl$core.num, " cores for DA test") )
    gym.output = DATest(count.matrix = count.matrix, 
                        phenotype = phenotype, 
                        core.num = gymControl$core.num, 
                        R = gymControl$R)
  } else if (method == "corncob"){
    source("../GymFunctions/corncob.R")
    gym.output = corncob(count.matrix = count.matrix, 
                         phenotype = phenotype, 
                         count.min = gymControl$count.min)
  } else if (method == "LDMtest"){
    source("../GymFunctions/LDMtest.R")
    gym.output = LDMtest(count.matrix = count.matrix, 
                         phenotype = phenotype)
  } else if (method == "ZINBDPP"){
    source("../GymFunctions/ZINBDPP.R")
    gym.output = ZINBDPP(count.matrix = count.matrix, 
                         phenotype = phenotype,
                         N.mcmc = gymControl$N.mcmc,
                         b = gymControl$b,
                         h = gymControl$h, 
                         count.min = gymControl$count.min, 
                         seed = gymControl$seed)
  }
  
  version.info = NULL
  if (method %in% c("DESeq2", "edgeR","metagenomeSeq", "corncob")){
    version.info = packageVersion(method)
  }else if (method == "LDMtest"){
    version.info = '2.1'
  }else if (method == "WaVEedgeR"){
    version.info = list(version.edgeR = packageVersion("edgeR"),
                        version.wave = packageVersion("zinbwave"))
  }else if (method == "WaVEDESeq2"){
    version.info = list(version.edgeR = packageVersion("edgeR"),
                        version.wave = packageVersion("zinbwave"))
  }else if (method ==  "DATest"){
    version.info = packageVersion("DAtest")
  }
   
  # add the test name information
  if(includeData){
    gym.output = list(model.name = method, 
                      model.details = gym.output, 
                      count.matrix = count.matrix, 
                      phenotype = phenotype,
                      version.info = version.info)
  }else{
    gym.output = list(model.name = method, 
                      model.details = gym.output,
                      version.info = version.info)
  }
  return(gym.output)
}

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
                      h = 50,
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




plot.ROC = function(response, predictor, method.name = NULL){
  # load packages #
  if (!require(pROC)) {install.packages("pROC", version = "1.15.3", repos="http://cran.us.r-project.org")}
  if (!require(ggplot2)) {install.packages("ggplot2", version = "3.2.1", repos="http://cran.us.r-project.org")}
  if (!require(cowplot)) {install.packages("cowplot", version = "1.0.0", repos="http://cran.us.r-project.org")}
  theme_set(theme_cowplot())
  
  # get ROC curve #
  roc.full = pROC::roc(response = response, predictor = predictor, quiet = T)
  roc.df = data.frame(xlab = 1 - roc.full$specificities,
                      ylab = roc.full$sensitivities)
  auc.val = as.numeric(auc(roc.full))
  # plot ROC curve #
  roc.plt = ggplot(data = roc.df,
                   aes(x=xlab, y=ylab))+
    geom_line(data = roc.df, 
              size = 1,
              color = 'black') +
    coord_fixed()+
    geom_segment(aes(x=0, y=0, xend=1, yend=1), colour="grey", 
                 linetype = "solid", size = 0.1)+
    xlab("False positive rate") + ylab("True positive rate") + 
    ggtitle(method.name) + 
    annotate("text", x = 0.75, y = 0.25,color = 'black', 
             fontface = 2,size = 5,
             label = paste("AUC:" ,sprintf("%0.3f", round(auc.val, digits = 3)))) +
    theme(plot.title = element_text(hjust = 0.5,size=20),
          axis.text=element_text(size=14),
          legend.position = "none",
          axis.title=element_text(size=15),plot.margin = unit(c(t = 0.5,r = 0,b = 0,l = 0), "cm")
    )
  plot(roc.plt)
  
}


# Description: 

#' Check the count matrix input

# Parameters:
#' @param count.matrix A sampe-by-taxon count matrix. Should be the original sequencing counts  

# Output:
#' @return count.matrix A sampe-by-taxon count matrix

check.countmatrix = function(count.matrix){
  if(! all(count.matrix == floor(count.matrix)) | any(count.matrix < 0)){
    stop("Elements in the input count matrix must be nonnegative integers")
  }
  if(class(count.matrix) == "data.frame"){
    count.matrix = data.matrix(count.matrix)
  }
  if( any(apply( count.matrix, 2, function(x){sum(x == 0) == length(x)})) ){
    warning("Some taxa only have 0s across all the samples. Certain differential abundance analysis methods may fail")
  }
  return(count.matrix)
}


# Description: 
#' Check the phenotype indicator vector input

# Parameters:
#' @param phenotype A phenotype indicator vector for all the samples (can have at most 2 types of phenotype)

# Output:
#' @return A vector of 0 if there is only one phenotype, or a binary (0, 1) vector if there are 2 phenotypes

check.phenotype = function(phenotype){
  if(length(unique(phenotype)) > 2 ){
    stop("Can handle at most 2 phenotypes")
  }
  if(length(unique(phenotype)) == 1 && unique(phenotype) != 0 ){
    phenotype = rep(0, length(phenotype))
  }
  if(length(unique(phenotype)) == 2 & any(unique(phenotype) != c(0, 1)) ){
    phenotype.input = phenotype
    phenotype = factor(phenotype, levels = c(unique(phenotype)[1], unique(phenotype)[2]), labels = c(0, 1) )
    phenotype = as.numeric(levels(phenotype))[phenotype]
    cat("samples labeled as", unique(phenotype.input)[1], "is converted to 0, and samples labeled as", unique(phenotype.input)[2], "is converted to 1. \n" )
  }
  return(phenotype)
}


# Description: 
#' Check the random seed input

# Parameters:
#' @param seed A random seed (positive integer)

# Output:
#' @return A random seed

check.seed = function(seed){
  seed.max = .Machine$integer.max
  if(seed < 0 | seed != floor(seed) | seed > seed.max){
    stop("random seed should be a positive integer with range of 1-",seed.max)
  }
  seed = as.numeric(seed)
  return(seed)
}



