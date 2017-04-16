# this contains utility functions that will be needed for all
# prediction functions

# function resamples the expression vector
# returns a list of response variable with shuffled rows for given column
# B is number of random samples 
resampleMat<-function(mat,col=1,B=1000){
  lapply(1:B,sample(mat[,col],nrow(mat)))
}

# function estimates gamma distribution based P-values from the null
# distribution obtained from resampling
estimateGammaPval<-function(){
  
}