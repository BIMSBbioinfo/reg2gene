# correlation and distance correlation(dcor) based methods to measure association
# dcor is used to measure linear and non-linear associations

source("targetPrediction/predUtils")

corMat<-function(mat,method="pearson",col=1){
  
  apply(mat[,-col],2,function(x,y,method){cor(x,y,method=method)},
        y=mat[,col],method=method)
}


dcorMat<-function(mat,method="pearson",col=1){
  require(energy)
  apply(mat[,-col],2,function(x,y){dcor(x,y)},
        y=mat[,col])
}

corResample<-function(mat,method="pearson",col=1,B=1000){
  
}

dcorResample<-function(mat,method="pearson",col=1,B=1000){
  
}