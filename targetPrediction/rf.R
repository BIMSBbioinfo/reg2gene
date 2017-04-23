# randomForest based prediction functions


#' internal RandomForest prediction function
#' 
#' private function that uses randomForest to predict fiven column from the rest
#' and assigns significance to
#' correlation coefficents based on resampling
#' 
#'
#' @param mat matrix of response and predictor variables
#' @param col column number of response variables
#' @param B number of random samples (shuffling) of response variables
#' @param ... further arguments to randomForest, not implemented yet
#'
#' @keywords internal
#' 
#' @example 
#' 
#' 
#' 
glmnetResample<-function(mat,col=1,B=1000,...){
  require(randomForest)
  
  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))
  
  # original coefs in this case importance values as %incMSE
  mod<- randomForest(x = mat[,-col], y = mat[,col],
                     importance=TRUE,na.action=na.omit)
  
  orig=importance(mod)[,1]
  
  #coefs from resampling
  coefs=matrix(0.0,ncol=ncol(mat[,-col]),nrow=(B))
  
  
  # calculate coeff/importance for resampled Ys
  for(i in 1:B){
    mod<- randomForest(x = mat[,-col], y = Ys[[i]],
                       importance=TRUE,na.action=na.omit)
    
    coefs[i,]=importance(mod)[,1]
  }
  
  # calculate p-vals 
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)
  
  # arrange output format and return
  return(pvals)
}


