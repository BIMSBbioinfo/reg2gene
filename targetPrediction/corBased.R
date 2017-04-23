# correlation and distance correlation(dcor) based methods to measure association
# dcor is used to measure linear and non-linear associations

source("predUtils.R")

#' Correlation between a col matrix and a vector
#' 
#' A convenience function that computes the the correlation of a given vector
#' with all the other column vectors of a matrix
#' @keywords internal
#' @example 
corMat<-function(y,mat,method="pearson"){
  
  apply(mat,2,function(x,y,method){cor(x,y,method=method)},
        y=y,method=method)
}

#' Distance correlation between a col matrix and a vector
#' 
#' a convenience function that computes the the distance correlation of a given 
#' vector (y) with all the other columns of a matrix (mat)
#'
#' @keywords internal
#' @example 
dcorMat<-function(y,mat){
  require(energy)
  apply(mat,2,function(x,y){dcor(x,y)},
        y=y)
}

#' Correlation with resampling p-values
#' 
#' private function that calculate correlation of given
#' column matrix with a column and assigns significance to
#' correlation coefficents based on resampling
#' 
#' @param mat 
#' @param method
#' @param col
#' @param B number of random samples (shuffling) of given column
#' @keywords internal
#' @example 
#' 
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#' 
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' corResample(m,method="pearson",col=1,B=1000)
#'
corResample<-function(mat,method="pearson",col=1,B=1000){
  
  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))
  
  # original coefs
  orig=corMat(y=mat[,col],mat=mat[,-col],method)
  
  #coefs from resampling
  coefs=matrix(0.0,ncol=ncol(mat[,-col]),nrow=(B))

  
  # calculate coeff for resampled Ys
  for(i in 1:B){
    
    coefs[i,] = corMat(Ys[[i]],mat=mat[,-col],method)
  }
  
  # calculate p-vals 
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)
  
  # arrange output format and return
  return(pvals)
}


#' Distance correlation with resampling p-values
#'
#' private function that calculates distance correlation of given
#' column matrix with a vector and assigns significance to
#' correlation coefficents based on resampling
#' @param mat 
#' @param col
#' @param B number of random samples (shuffling) of given column
#' @keywords internal
#' @example 
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#' 
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' dcorResample(m,col=1,B=1000) 
dcorResample<-function(mat,col=1,B=1000){
  
  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))
  
  # original coeff
  orig=dcorMat(y=mat[,col],mat=mat[,-col])
  
  #coefs
  coefs=matrix(0.0,ncol=ncol(mat[,-col]),nrow=(B))

  
  # calculate coeff for resampled Ys
  for(i in 1:B){
    
    coefs[i,] = dcorMat(Ys[[i]],mat=mat[,-col])
  }
  
  # calculate p-vals 
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)
  
  # arrange output format and return
  return(pvals)
}