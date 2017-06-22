
#' Associate regulatory regions to genes
#' 
#' The function associates/links regulatory regions to genes based
#' on regulatory activity and gene expression relationship. 
#' 
#' @param input a GRangesList that contains regulatory activity and 
#' gene expression results. It has to have specific columns, see the
#' output of \code{regActivityAroundTSS} function.
#' @param method one of "pearson","spearman","dcor","elasticnet",
#' "randomForest". Default: "pearson". See Details for more.
#' @param tag a string to be used to append to results from methods. This string
#'  will be attached to "n","coefs","pval" columns. Default NULL.
#'  
#' @return 
#' 
#' @examples 
#' 
#' @details 
#' 
#' @author Altuna Akalin
#' 
associateReg2Gene<-function(input,method="pearson",tag=NULL){
  
  # a function to create input ranges to output ranges
  # later p-values, effect sizes etc will appended to this
  # GRanges object
  gr=grlist2gr(input)
  
  # decide which method you want to apply and call that
  # function
  if(method=="pearson"){
    res=mclapply(input,function(gr){
      # get matrix of activity and expressions
      mat=t(as.matrix(mcols(gr)[,-(1:3)]))
      
      # get coef and pvalues, add number of samples to the result
      cbind(n=nrow(mat),t(corResample(mat,method="pearson",col=1,B=1000)))
      })
  }else if(method == "spearman"){
    res=mclapply(input,function(gr){
      # get matrix of activity and expressions
      mat=t(as.matrix(mcols(gr)[,-(1:3)]))
      cbind(n=nrow(mat),t(corResample(mat,method="pearson",col=1,B=1000)))
    })
  }else if(method == "dcor"){
    res=mclapply(input,function(gr){
      # get matrix of activity and expressions
      mat=t(as.matrix(mcols(gr)[,-(1:3)]))
      cbind(n=nrow(mat),t(dcorResample(mat,col=1,B=1000)))
    })
  }else if(method == "elasticnet"){
    stop("elasticnet not implemented yet")
  }else if(method == "randomForest"){
    stop("randomForest not implemented yet")
  }
  
  # combine stats with Granges
  comb.res=do.call("rbind",res)[,1:3]
  rownames(comb.res)=NULL
  if(!is.null(tag)){
    colnames(comb.res)=paste(colnames(comb.res),tag,sep=".")
  }
  mcols(gr)=cbind(mcols(gr),DataFrame(comb.res) )
  
  # return result
  gr
}

#' convert GRangesList to GRanges with associated regulatory regions (internal)
#' 
#' The function converts output of \code{\link{regActivityAroundTSS}} (a
#' GRangesList) to a GRanges object by rearranging the rows and columns so
#' that output GRanges object mainly contains gene name and position with
#' an additional column "reg" which contains the regulatory region coordinates 
#' 
#' @param grlist a GRangesList output from \code{\link{regActivityAroundTSS}}
#' 
#' @example 
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
  
  grpairs
  }))
}
