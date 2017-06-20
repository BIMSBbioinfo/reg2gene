# this contains utility functions that will be needed for all
# prediction functions



#' Estimate P-values from resampling statistics using Gamma distribution
#' 
#' The function estimates gamma distribution based P-values from the null
#' distribution obtained from resampling, a gamma distribution is fit
#' to the null distribution and p-values are calculated based on that fitted
#' distribution
#' 
#' @param vals a vector or matrix of column vectors. Each column corresponds
#' to coefficients obtained from resampling fit
#' @param orig a single value or a vector of original coefficients obtained from the
#' original matrix 
#' @param abs if TRUE, use absolute values of coefficients to fit a distribution
#' @param add an numeric value defaults to 0. if more than 0, add that value to
#' orig and vals. This way gamma distribution can model null distributions that
#' have negative values.
#' 
#' @importFrom fitdistrplus fitdist
#' @keywords internal
estimateGammaPval<-function(vals,orig,abs=TRUE,add=0){
  
  orig2=(orig) # keep original values
  
  # if needed take absolute values
  if(abs){
    vals=abs(vals)
    orig2=abs(orig2)
  }
  if(add){ # adding numbers might be needed, gamma can only modep positive values
    vals=vals+add
    orig=orig+add
  }
  
  if(class(vals) != "matrix"){
  
    if(any(is.na(vals))){
      result=c(coefs=NA,pval=NA,p2=NA)
    
    }else{
      param <-  fitdistrplus::fitdist(vals,"gamma",method="mme") # fit gamma
      # p-val based on gamma
      pval=1-pgamma(orig2,param$estimate[1],param$estimate[2]) 
    
      # if 0 p-val returned set it to minimum possible or another option is
      # 2.2e-16 
      if(pval==0){pval=1-pgamma(max(vals),param$estimate[1],param$estimate[2])}
    
      # result vector
      result=c(coefs=orig,pval=pval,p2=sum(vals >= orig2)/length(vals))}
    
  }else{
    
    # same as above it is just when the input is a matrix
    result=mapply(function(x,orig2,orig){
      if(any(is.na(x)) | all(x==0)){
        c(coefs=NA,pval=NA,p2=NA)
        
      }else{
      param <-  fitdistrplus::fitdist(x,"gamma",method="mme")
      
      pval=1-pgamma(orig2,param$estimate[1],param$estimate[2])
      if(pval==0){pval=1-pgamma(max(x),param$estimate[1],param$estimate[2])}
      c(coefs=orig,pval=pval,p2=sum(x >= orig2)/length(x) )
      }
    },split(t(vals),1:ncol(vals)),orig2,orig)
    
  }
}

#' decide if the input colum matrix has zero variation
#' 
#' The input column matrix represent gene expression and regulatory region
#' activities. If it has zero variation for gene expression 
#' 
#' @param mat a column matrix that contains gene expression values and 
#'            regulatory region activities. 
#' @param col  column number of response variables
#' @example 
#' 
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#' 
#' @keywords internal
#' @author Altuna Akalin
zeroVar<-function(mat,col=1){
  
  # flag 0 variation in gene expression
  if(sd(mat[,col])==0){
    return(TRUE)
  }
  FALSE
}

#' decide if the input colum from a matrix has many zeros
#' 
#' The input column matrix represent gene expression and regulatory region
#' activities. decides if the genes
#' 
#' @param mat a column matrix that contains gene expression values and 
#'            regulatory region activities. 
#' @param col  column number of response variables
#' @example 
#' 
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#' 
#' @keywords internal
#' @author Altuna Akalin
manyZeros<-function(mat,col=1){
  # if many values are zero cross-validation or subsampling based techniques
  # will fail
  if( sum(mat[,col]==0)/nrow(mat) > 0.9 ){
    return(TRUE)
  }
  FALSE
}


