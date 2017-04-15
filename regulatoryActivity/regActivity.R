
#' Measures regulatory activity over pre-defined regions
#' 
#' The function calculates regulatory activity from histone
#' modification, DNAse or methylation signals.
#' 
#' @param regRegions a GRanges object that contains regulatory regions
#' over which the regulatory activity will be calculated.
#' @param activitySignals a named list of BigWig files. Names correspond to 
#'        unique sample ids or names
#' @param isCovNA (def:FALSE), if this is set to TRUE, uncovered
#' bases are set to NA, this is important when dealing with methylation
#' data.
#' @param summaryOperation "mean"(default),"median" or "sum". This
#' designates which summary operation should be used over the regions
#' @param normalize NULL(default). If set to "quantile" returned activity
#' measures are quantile normalized
#' 
#' @return a GRanges object where its meta-columns correspond
#'         to calculated acvitity measures and column names 
#'         correspond to provided sample ids or names.
#' 
#' @import genomation #NB list any other packages needed
#' 
#' @details           
#'
#' @examples #NB 10 regions and over 2 bw files provide small examples that work on beast
#' 
regActivity<-function(regRegions,activitySignals,
                      isCovNA=FALSE,summaryOperation="mean",
                      normalize=NULL){
  
}