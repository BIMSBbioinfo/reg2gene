# correlation and distance correlation(dcor) based methods to measure association
# dcor is used to measure linear and non-linear associations

source("targetPrediction/predUtils")

# a convenience function that computes the the correlation of a given column
# with all the other columns
corMat<-function(y,mat,method="pearson"){
  
  apply(mat,2,function(x,y,method){cor(x,y,method=method)},
        y=y,method=method)
}

# a convenience function that computes the the distance correlation of a given 
# column with all the other columns
dcorMat<-function(y,mat){
  require(energy)
  apply(mat,2,function(x,y){dcor(x,y)},
        y=y)
}

# private function that calculate correlation of given
# column with other columns and assigns significance to
# correlation coefficents based on resampling
# @param mat 
# @param method
# @param col
# @param B number of random samples (shuffling) of given column
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

# private function that calculats distance correlation of given
# column with other columns and assigns significance to
# correlation coefficents based on resampling
# @param mat 
# @param col
# @param B number of random samples (shuffling) of given column
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