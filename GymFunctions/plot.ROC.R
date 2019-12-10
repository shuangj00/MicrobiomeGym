plot.ROC = function(response, predictor, method.name = NULL){
  # load packages #
  if (!require(pROC)) {install.packages("pROC", version = "1.15.3", repos="http://cran.us.r-project.org")}
  if (!require(pROC)) {install.packages("ggplot2", version = "3.2.1", repos="http://cran.us.r-project.org")}
  if (!require(pROC)) {install.packages("cowplot", version = "1.0.0", repos="http://cran.us.r-project.org")}
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