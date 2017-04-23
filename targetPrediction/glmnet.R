# glmnet based prediction functions


#' internal function that uses glmnet for prediction
#' 
#' private function that uses glmnet to predict fiven column from the rest
#' and assigns significance to
#' correlation coefficents based on resampling
#' @param mat matrix of response and predictor variables
#' @param col column number of response variables
#' @param B number of random samples (shuffling) of response variables
#' @param ... further arguments to glmnet, not implemented yet
#' @keywords internal
#' @importFrom glmnet glmnet
#' 
#' @example 
#' 
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#' 
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' glmnetResample(scale(m),col=1,B=1000)
glmnetResample<-function(mat,col=1,B=1000,...){
  require(glmnet)
  
  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))
  
  # original coefs
  mod<- cv.glmnet(x = mat[,-col], y = mat[,col],
                  standardize=FALSE,
                  nfolds=5,alpha=0.5)
  
  orig=coef(mod,s="lambda.min")[-1,1]
  
  #coefs from resampling
  coefs=matrix(0.0,ncol=ncol(mat[,-col]),nrow=(B))
  
  
  # calculate coeff for resampled Ys
  for(i in 1:B){
    
    mod<- cv.glmnet(x = mat[,-col], y = Ys[[i]],
                    standardize=FALSE,
                    nfolds=5,alpha=0.5)
    
    coefs[i,]=coef(mod,s="lambda.min")[-1,1]
  }
  
  # calculate p-vals 
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)
  
  # arrange output format and return
  return(pvals)
}

  
