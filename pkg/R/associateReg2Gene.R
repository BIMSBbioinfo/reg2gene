#' Associate regulatory regions to genes in CHUNKS
#'
#' The function runs \code{\link{associateReg2Gene}} function for chunks of 
#' genes to avoid reaching a memory quota in R
#'
#' @param input a GRangesList that contains regulatory activity and
#' gene expression results. It has to have specific columns, see the
#' output of \code{\link{associateReg2Gene}} function.
#' @param out.dir (character) Path to the directory where intermediate results
#' will be stored.
#' @param method one of "pearson","spearman","dcor","elasticnet",
#' "randomForest". Default: "pearson". 
#' @param tag a string to be used to append to results from methods. This string
#'  will be attached to "n","coefs","pval" columns. Default NULL.
#' @param scale if TRUE (default) the the values for gene expression and
#' regulatory activity will be scaled to have 0 mean and unit variance using
#' \code{\link[base]{scale}} function.
#' @param cores number of cores to be used
#' @param B number of randomizations, default 1000. This procedure
#' is used to estimate
#' P-values for coeficients returned by different methods.
#' @param chunks Default: 1 - everything is analyzed together. If >1 then
#' regulatory regions are divided into N chunks (useful for >>N of genes/
#' regulatory regions to avoid reaching R memory quota)
#' @param saveTag (character, defult: "ChunkAnalysis") Unique naming for created
#' .rds files which correspond to different gene chunks. This character is used 
#' to pull together all .rds files from one round of analysis.
#' @param ... further arguments to methods, not implemented yet
#' @param remove.chunks (default T). To remove all chunks of .rds data created
#' while running a function
#' @return a GRanges object containing every potential association and
#' between a regulatory region and TSS, and the estimated association statistic
#' , its P-value and Q-value.
#' @examples
#' \dontrun{
#' regAcFile=system.file("extdata", "sampleRegActivityAroundTSS.rds",
#'                       package = "reg2gene")
#' input=readRDS(regAcFile)
#' b=associateReg2Gene(input)
#'
#' }
#'
#' @details
#' The function implements \code{\link{associateReg2Gene}} function, but 
#' regulatory regions are dividied into predefined number of chunks; 
#' anaysis is performed separately for each chunk; and results are saved in 
#' predefined output directory. After analysis, all results are pulled
#' together, Q-value is calculated, and all intermediate files are removed.
#'
#' @author Inga Patarcic
#'
#' @importFrom doMC registerDoMC
#' @import foreach
#' @import stringr 
#' @importFrom qvalue qvalue
#' @export
associateReg2GeneChunks <- function(input,
                                    out.dir,
                                    method="pearson",
                                    tag=NULL,
                                    saveTag="ChunkAnalysis",
                                    scale=TRUE,
                                    cores=1,
                                    B=1000,
                                    chunks=1,
                                    remove.chunks=T,
                                    ...){
  
  # save warnings, output,etc. 
  sink(paste0(out.dir,"log.txt"), append=TRUE, split=TRUE)
  
  # splitting in gene chunks
  Split.factor <- split(1:length(input), sort(1:length(input)%%chunks))
  
  
  # run associateReg2Gene for gene chunks and save them separately
  lapply(Split.factor,function(x){
    print(x[1])
    
    # adjusting for the fact that code breakes from time to time
    
    counter=1
      while(counter<=10){
        Res_associateReg2Gene <- try(associateReg2Gene(input[x],
                                                       method=method,
                                                       B=B,
                                                       tag=tag,
                                                       scale=scale,
                                                       cores=cores),silent=TRUE)
        if(!is(Res_associateReg2Gene, 'try-error')) break
        if(is(Res_associateReg2Gene, 'try-error')) {counter=counter+1}
      }
    
    # report error even after 10 iterations of trying
    if (is(Res_associateReg2Gene, 'try-error')&counter>10){print(paste0(out.dir,
                                                          saveTag,x[1],
                                                    ".rds was not successful"))}
    
    saveRDS(Res_associateReg2Gene,paste0(out.dir,saveTag,x[1],".rds"))})
  
  
  # input all associateReg2Gene gene chunks in correct order
  GeneChunks <- sort(list.files(out.dir,pattern = saveTag,full.names = T,
                                recursive=T))
  GeneOrder <- order(as.numeric(str_replace(str_extract(GeneChunks,
                                                  "[0-9]*.rds"),".rds",""))) 
  
  GeneEnhancerResults <- do.call(c,lapply(GeneChunks[GeneOrder], 
                                          function(x) readRDS(x)))
  
  if (typeof(GeneEnhancerResults)=="list"){
    GeneEnhancerResults <- do.call(c,GeneEnhancerResults)}
  
  # calculating appropriate Qvalue for all reported Pvalues
  # assumption Pvalue is always stored in the 6th column
  
  QValue <- qvalue(mcols(GeneEnhancerResults)[,6])$qvalues
  mcols(GeneEnhancerResults)[,7] <- QValue
  
  # remove all input files
  if (remove.chunks==T){lapply(GeneChunks,function(x){system(paste0("rm ",x))})}
  
  return(GeneEnhancerResults)
  
  sink()
  
  
} 

#' Associate regulatory regions to genes
#'
#' The function associates/links regulatory regions to genes based
#' on regulatory activity and gene expression relationship.
#'
#' @param input a GRangesList that contains regulatory activity and
#' gene expression results. It has to have specific columns, see the
#' output of \code{\link{regActivityAroundTSS}} function.
#' @param method one of "pearson","spearman","dcor","elasticnet",
#' "randomForest". Default: "pearson". See Details for more.
#' @param tag a string to be used to append to results from methods. This string
#'  will be attached to "n","coefs","pval" columns. Default NULL.
#' @param scale if TRUE (default) the the values for gene expression and
#' regulatory activity will be scaled to have 0 mean and unit variance using
#' \code{\link[base]{scale}} function.
#' @param cores number of cores to be used
#' @param B number of randomizations, default 1000. This procedure
#' is used to estimate
#' P-values for coeficients returned by different methods.
#' @param chunks. Default: 1 - everything is analyzed together. Useful for >>N 
#' of  genes/regulatory regions R memory 
#' @param ... further arguments to methods, not implemented yet
#' @return a GRanges object containing every potential association and
#' between a regulatory region and TSS, and the estimated association statistic
#' , its P-value and Q-value.
#'
#' @examples
#' \dontrun{
#' regAcFile=system.file("extdata", "sampleRegActivityAroundTSS.rds",
#'                       package = "reg2gene")
#' input=readRDS(regAcFile)
#' b=associateReg2Gene(input)
#'
#' }
#'
#' @details
#' The function implements four methods to associate regulatory
#' activity to genes: correlation (pearson or spearman), distance
#' correlation ("dcor"), elasticnet and random forests. Simple linear regression on
#' scaled data is equivalent to Pearson correlation. Distance correlation is
#' a metric for statistical dependence between two variables. It can capture
#' non-linear relationships different from correlation.
#' Elasticnet and randomForests are methods that
#' can capture additivity between regulatory activities and their relationship
#' to gene expression. For all the methods, P-values are estimated by shuffling
#' the gene expression vector \code{B} times and getting a null distribution
#' for the estimated statistic, and comparing the original statistic to the null
#' distribution. Estimated statistics are
#' correlation coefficients for correlation and distance correlation methods,
#' Regression coefficients for elasticnet and gini importance scores for random
#' forests.
#'
#' @author Altuna Akalin
#'
#' @importFrom doMC registerDoMC
#' @import foreach 
#' @importFrom qvalue qvalue
#' @export
associateReg2Gene<-function(input,method="pearson",tag=NULL,scale=TRUE,
                            cores=1,B=1000,...){

  # drop NULL genes in the list
  nulls=which(sapply(input,is.null))
  if(length(nulls)>0){
      input=input[-c(nulls)]
  }

  # set up cores for multicore shit
  if(cores>2){
    registerDoMC(cores)
  }

  # decide which prediction/association method you want to use
  if(method=="pearson"){
    model.func=quote(corResample(mat,scale,method="pearson",col=1,B=B) )
  }else if(method == "spearman"){
    model.func=quote( corResample(mat,scale,method="pearson",col=1,B=B) )
  }else if(method == "dcor"){
    model.func=quote( dcorResample(mat,scale,col=1,B=B) )
  }else if(method == "elasticnet"){
    model.func=quote( glmnetResample(mat,scale,col=1,B=B) )
  }else if(method == "randomForest"){
    model.func=quote( rfResample(mat,scale,col=1,B=B) )
  }

  # call the function in a foreach loop
  res=foreach(gr=input) %dopar%
  {
          

      mat=t(as.matrix(mcols(gr)[,-which(names(mcols(gr)) %in% c("name","name2",
                                                          "featureType"))]))
        colnames(mat) <-  gr$name
    
        
      # only work with complete cases, hope there won't be
      # any columns with only NAs
      # mat=mat[complete.cases(mat),]
      
  # se complete.cases (remove cells, remove enhancers, remove random NA)  
      mat <- stepwise.complete.cases(mat)


      # adjusting for the case where all rows and columns are filtered OUT
     
      x  <- data.frame(matrix(NA,nrow=length(gr$name[-1]),ncol=4,byrow = T,
                              dimnames=list(gr$name[-1],
                                       c("n","coefs","pval","qval"))))
            
      
      # necessary for genes with one enhancer that is filtered out due to NA 
          if (is.matrix(mat)==F) {x <- x
          
          # get coef and pvalues, add number of samples to the result
          # adjusted for data with only one enhancer tested 
   }else if (nrow(mat)==2){x <- data.frame(cbind(n=nrow(mat),
                                                  t(eval(model.func))))
                        colnames(x) <- c("n","coefs","pval","pval2")
      
    }else if ((nrow(mat)>2)&(ncol(mat)==length(gr$name))){
              x <- cbind(n=nrow(mat),t(eval(model.func)))
              colnames(x) <- c("n","coefs","pval","pval2")
            
    }else if ((nrow(mat)>2)&(ncol(mat)!=length(gr$name))){
        # adjusting for the case where just some enhancers are filtered any(NA) 
        x[colnames(mat)[-1],] <- cbind(n=nrow(mat), t(eval(model.func)))
          
    }
      
   return(x)
      
    }
      
  

  # combine stats with Granges
  comb.res=do.call("rbind.data.frame",res)[,1:4]
  rownames(comb.res)=NULL # needed for GRanges DataFrame

  # add qvalue calculations
  comb.res <- qvaluCal(comb.res)
    
  
  # if a tag for column names given, add it
  if(!is.null(tag)){
    colnames(comb.res)=paste(colnames(comb.res),tag,sep=".")
  }

  # a function to create input ranges to output ranges
  # later p-values, effect sizes etc will appended to this
  # GRanges object
  gr2=grlist2gr(input)
  mcols(gr2)=cbind(mcols(gr2),DataFrame(comb.res ) )

  # return result
  gr2
}




#########-------------------------------------------------###########
########### utility functions for association prediction ############

# this part contains utility functions that will be needed for all
# prediction functions


#' convert GRangesList to GRanges with associated regulatory regions (internal)
#'
#' The function converts output of \code{\link{regActivityAroundTSS}} (a
#' GRangesList) to a GRanges object by rearranging the rows and columns so
#' that output GRanges object mainly contains gene name and position with
#' an additional column "reg" which contains the regulatory region coordinates
#'
#' @param grlist a GRangesList output from \code{\link{regActivityAroundTSS}}
#'
#' @examples
#'
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#'
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' glmnetResample(scale(m),col=1,B=1000)
#'
#' @keywords internal
#' @author Altuna Akalin
grlist2gr<-function(grlist){

  unlist(endoapply(grlist, function(x){

    reg=granges(x[x$featureType != "gene",])
    gene=x[x$featureType == "gene",c("name" ,"name2")]

    grpairs=rep(gene,length(reg))
    grpairs$reg=reg
    
    # ordering meta-data to have correct input for benchmark()
    mcols(grpairs) <-  mcols(grpairs)[c("reg","name","name2")]
    
    grpairs
  }))
}



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
#' @examples
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
#' @examples
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

#' stepwise complete cases method
#'
#' Removes NA but in stepwise order: 1st cells that have NA for at least
#' 3/4 of enhancer regions, then the same thing is applied for cell types
#' then individual NAs are removed
#'
#' @param mat a column matrix that contains gene expression values and
#'            regulatory region activities.
#'
#' @keywords internal
#' @author Inga Patarcic
stepwise.complete.cases <- function(x){
  
  similarity=0.75
  
  remove.cc <- function(x,similarity){
    n.na <- sum(is.na(x))
    if (n.na>=similarity*length(x)){return(F)}
    if (n.na<similarity*length(x)){return(T)}
  }
  
  # remove problematic cells, then enhancers then remove remaining NAs 
  x <- x[apply(x,1,remove.cc,similarity=similarity),]
  x <- x[,apply(x,2,remove.cc,similarity=similarity)]
  
  if (is.matrix(x)){x <- x[,apply(x,2,remove.cc,similarity=0.00001)]}
  
  return(x)
  
}

#### END OF utility functions for association prediction ####
#########------------------------------------------------############



#########-------------------------------------------------###########
########### association prediction functions

#' Correlation between a col matrix and a vector
#'
#' A convenience function that computes the the correlation of a given vector
#' with all the other column vectors of a matrix
#' @keywords internal
corMat<-function(y,mat,method="pearson"){

  if (!is.vector(mat)){
  return(apply(mat,2,function(x,y,method){cor(x,y,method=method)},
        y=y,method=method))}
  
  if (is.vector(mat)){return(cor(mat,y,method=method))}
  
}

#' Distance correlation between a col matrix and a vector
#'
#' a convenience function that computes the the distance correlation of a given
#' vector (y) with all the other columns of a matrix (mat)
#'
#' @keywords internal
#' @importFrom energy dcor
dcorMat<-function(y,mat){
  #require(energy)
  if (!is.vector(mat)){return(apply(mat,2,
                                    function(x,y){energy::dcor(x,y)},y=y))}
  
  if (is.vector(mat)){return(dcor(mat,y))}
}

#' Correlation with resampling p-values
#'
#' private function that calculate correlation of given
#' column matrix with a column and assigns significance to
#' correlation coefficents based on resampling
#'
#' @param mat matrix of response and predictor variables
#' @param col column number of response variable
#' @param method pearson or spearman
#' @param B number of random samples (shuffling) of given column
#' @keywords internal
#' @examples
#'
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#'
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' corResample(m,method="pearson",col=1,B=1000)
#'
corResample<-function(mat,scale,method="pearson",col=1,B=1000){

  # decide if the matrix needs to be dropped and NA returned due
  # to low or zero variation in gene expression
  if( zeroVar(mat) | manyZeros(mat) ){
    return(matrix(NA,ncol=ncol(mat)-1,nrow=3,
                  dimnames=list(c("coefs","pval","p2"),1:(ncol(mat)-1) ))
    )
  }
  
  # scale if necessary
  if(scale){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }

  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))

  # original coefs
  # if only 1 enhancer available
  # if more than one enhancer available
  orig=corMat(y=mat[,col],mat=mat[,-col],method)

  #coefs from resampling (adjusted for 1 or more enhancers)
  Nenh <- ncol(mat[,-col])
    if (is.null(Nenh)){coefs=matrix(0.0,ncol=1,nrow=(B))}
    if (!is.null(Nenh)){coefs=matrix(0.0,ncol=Nenh,nrow=(B))}

  # calculate coeff for resampled Ys
  for(i in 1:B){
    #cat(i,"\n")
    coefs[i,] = corMat(Ys[[i]],mat=mat[,-col],method)
  }

  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}


#' Distance correlation with resampling based p-values
#'
#' private function that calculates distance correlation of given
#' column matrix with a vector and assigns significance to
#' correlation coefficents based on resampling
#' @param mat matrix of response and predictor variables
#' @param col column number of response variable
#' @param B number of random samples (shuffling) of given column
#' @keywords internal
#' @examples
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#'
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' dcorResample(m,col=1,B=1000)
dcorResample<-function(mat,scale,col=1,B=1000){

  # decide if the matrix needs to be dropped and NA returned due
  # to low or zero variation in gene expression
  if( zeroVar(mat) | manyZeros(mat) ){
    return(matrix(NA,ncol=ncol(mat)-1,nrow=3,
                  dimnames=list(c("coefs","pval","p2"),1:(ncol(mat)-1) ))
    )
  }

  if(scale){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }
  
  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))

  # original coeff
  orig=dcorMat(y=mat[,col],mat=mat[,-col])

  #coefs (adjusted for 1 or more enhancers)
  Nenh <- ncol(mat[,-col])
  if (is.null(Nenh)){coefs=matrix(0.0,ncol=1,nrow=(B))}
  if (!is.null(Nenh)){coefs=matrix(0.0,ncol=Nenh,nrow=(B))}
  


  # calculate coeff for resampled Ys
  for(i in 1:B){
    #cat(i,"\n")
    coefs[i,] = dcorMat(Ys[[i]],mat=mat[,-col])
  }

  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}



#' internal function that uses glmnet for prediction
#'
#' private function that uses glmnet to predict fiven column from the rest
#' and assigns significance to coefficents based on resampling
#' @param mat matrix of response and predictor variables
#' @param col column number of response variable
#' @param B number of random samples (shuffling) of response variables
#' @param ... further arguments to glmnet, not implemented yet
#' @keywords internal
#' @importFrom glmnet glmnet
#'
#' @examples
#'
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#'
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' glmnetResample(scale(m),col=1,B=1000)
glmnetResample<-function(mat,scale,col=1,B=1000,...){
  #require(glmnet)

  # decide if the matrix needs to be dropped and NA returned due
  # to low or zero variation in gene expression
  if( zeroVar(mat) | manyZeros(mat) | ncol(mat)<=2){
    return(matrix(NA,ncol=ncol(mat)-1,nrow=3,
                  dimnames=list(c("coefs","pval","p2"),1:(ncol(mat)-1) ))
    )
  }
  
  
  if(scale){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }

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
  # catching errors if in resampling all 0s selected 
  
  tryCatch({
  for(i in 1:B){

    mod<- cv.glmnet(x = mat[,-col], y = Ys[[i]],
                    standardize=FALSE,
                    nfolds=5,alpha=0.5)

    coefs[i,]=coef(mod,s="lambda.min")[-1,1]
  }
  },error=function(coefs){coefs <- matrix(NA,ncol=ncol(mat[,-col]),nrow=(B))})

  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}


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
#' @importFrom ranger ranger
#' @examples
#' m=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
#'
#' m=apply(m, 2,as.numeric) # change to numeric matrix
#' rfResample(scale(m),col=1,B=1000)
#'
rfResample<-function(mat,scale,col=1,B=1000,...){
  #require(randomForest)
  #require(ranger)

  # decide if the matrix needs to be dropped and NA returned due
  # to low or zero variation in gene expression
  if( zeroVar(mat) | manyZeros(mat) | ncol(mat)<=2){
    return(matrix(NA,ncol=ncol(mat)-1,nrow=3,
                  dimnames=list(c("coefs","pval","p2"),1:(ncol(mat)-1) ))
    )
  }

  if(scale){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }

  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))
  # original coefs in this case importance values as %incMSE
  #mod<- randomForest(x = mat[,-col], y = mat[,col],
  #                   importance=TRUE,na.action=na.omit)
  #  orig=importance(mod)[,2] # gini index
  
  colnames(mat) <- c("gene",rep("enh",ncol(mat)-1))
  
  mod <- ranger( "gene ~ .", data = data.frame(mat), importance="impurity",
                num.trees = 500,write.forest = FALSE,
                num.threads=1)
  orig=mod$variable.importance # gini index
  #coefs from resampling
  coefs=matrix(0.0,ncol=ncol(mat[,-col]),nrow=(B))


  # calculate coeff/importance for resampled Ys
  # handling case where by intense sampling occasionally only zeros are selected
  tryCatch({
  for(i in 1:B){
    #mod<- randomForest(x = mat[,-col], y = Ys[[i]],
    #                   importance=TRUE,na.action=na.omit)

    #coefs[i,]=importance(mod)[,2]

    mod <-ranger( y ~ ., data = data.frame(y=Ys[[i]],mat[,-col]),
                  importance="impurity",
                  num.trees = 500,write.forest = FALSE)
    coefs[i,]=mod$variable.importance # gini index

  }
  },error=function(coefs){coefs <- matrix(NA,ncol=ncol(mat[,-col]),nrow=(B))})

  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}


#' internal qvalue calculation function
#'
#' private function calculates qvalue for 2 categories of pvalue
#' @importFrom qvalue qvalue
#' @keywords internal

qvaluCal <- function(comb.res){
  
    # add qvalues for pval1&2 - adjusted for filtered results
  
  qvalR <- matrix(NA,ncol=2,nrow=nrow(comb.res),
                  dimnames =  list(1:nrow(comb.res),c("qval", "qval2")))
      comb.res=cbind(comb.res,qvalR)
    
    # add qval for pval1&2  
      test <- !all(is.na(comb.res[,3])) # test if any NA present
    if (test){comb.res$qval=qvalue::qvalue(comb.res[,3])$qvalues}
      test <- !all(is.na(comb.res[,4]))  # test if any NA present
    if (test) {comb.res$qval2=qvalue::qvalue(comb.res[,4])$qvalues}

  return(comb.res)
  }
  
  ########### END OF association prediction functions
#########-------------------------------------------------###########
