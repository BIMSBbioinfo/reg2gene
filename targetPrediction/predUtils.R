# this contains utility functions that will be needed for all
# prediction functions



# function estimates gamma distribution based P-values from the null
# distribution obtained from resampling
estimateGammaPval<-function(vals,orig,abs=TRUE,add=0){
  
  
  if(abs){
    vals=abs(vals)
  }
  if(add){
    vals=vals+add
  }
  
  if(class(vals) != "matrix"){
    param <-  fitdistrplus::fitdist(vals,"gamma") # fit gamma
  
    pvals=1-pgamma(orig,param$estimate[1],param$estimate[2]) # p-val based on gamma
    
    # result vector
    result=c(coefs=orig,pvals=pvals,p2=sum(vals >= orig)/length(vals))
  }else{
    
    # same as above it is just when the input is a matrix
    result=mapply(function(x,orig){
      param <-  fitdistrplus::fitdist(x,"gamma")
      
      pvals=1-pgamma(orig,param$estimate[1],param$estimate[2])
      c(coefs=orig,pvals=pvals,p2=sum(x >= orig)/length(x) )
    },split(t(vals),1:ncol(vals)),orig)
    
  }
}

