
#' Link regulatory regions to genes
#' 
#' The function associates/links regulatory regions to genes based
#' on regulatory activity and gene expression relationship. 
#' 
#' @param input a GRangesList that contains regulatory activity and 
#' gene expression results. It has to have specific columns, see the
#' output of \code{regActivityAroundTSS} function.
#' @param method one of "pearson","spearman","dcor","elasticnet",
#' "randomForest". Default: "pearson". See Details for more.
#' 
#' @return 
#' 
#' @examples 
#' 
#' @details 
#' 
#' 
linkReg2Gene<-function(input,method="pearson"){
  
}