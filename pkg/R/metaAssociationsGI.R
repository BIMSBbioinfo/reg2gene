# functions for meta analysis

#' meta-analysis for regulatory region and gene associations
#'
#' The function combines association P-values and coefficients from different 
#' data sets using Fisher's method and weigthed averaging respectively. 
#' It is useful to combine datasets produced by different research groups while
#' avoiding problems such as batch effects. 
#' For example, it can be used to combine Roadmap and Blueprint 
#' \code{\link{associateReg2Gene}} gene~enhancer single-model information 
#' obtained using across cell RNA-Seq signals and CHiP-Seq H3K27ac tracks.
#' Aggregating single-model information (gene-enhancer pairs) by means of 
#' meta-analysis across different data-sources should increase the statistical 
#' power, improve the precision and accuracy of estimates and altogether 
#' produce more robust and reproducible results. After,combining P-values, 
#' q-values are calculated using the \code{\link[qvalue]{qvalue}} function.
#'
#' @param associations A list of \code{\link[InteractionSet]{GInteractions}}
#'  objects outputed from \code{\link{associateReg2Gene}} function.
#'
#'
#' @return returns a \code{\link[InteractionSet]{GInteractions}} object with an 
#' updated statistics from the meta-analysis. The output will be similar to
#' the output of \code{\link{associateReg2Gene}} function, it
#' will have the same meta-data columns in the same order.
#'
#' @examples # creating datasets
#' 
#' require(GenomicRanges)
#' require(InteractionSet)
#' 
#' gr2 <- gr <- GRanges(seqnames=rep("chr1",3),IRanges(1:3,3:5))
#'    x <- 1:5
#'    y <- 2:6
#'    z <- 10:14
#'    a <- rep(0,length(x))
#'    
#'    
#'    GeneInfo <- as.data.frame(matrix(c(rep("gene",3),rep("regulatory",6)),
#'                ncol = 3,byrow = TRUE),stringsAsFactors=FALSE)
#'                colnames(GeneInfo) <- c("featureType","name","name2")
#'  
#'  mcols(gr) <- DataFrame(cbind(GeneInfo,rbind(x,y,z)))
#'  mcols(gr2) <- DataFrame(cbind(GeneInfo,rbind(x,y,a)))
#'  
#'  # create associateReg2Gene output objects, GInteractions will all 
#'  # output results
#'  
#'  AssocObject <- reg2gene::associateReg2Gene(gr)
#'  AssocObject2 <- reg2gene::associateReg2Gene(gr2)
#'  
#'  # input for meta-analysis is list of such objects
#'  
#'  associations <- list(AssocObject,AssocObject2)
#'  names(associations) <- c("AssocObject","AssocObject2")
#'  
#'  # Run metaA
#'  metaAssociations(associations)
#'
#' @importFrom  data.table data.table
#' @importFrom  stats pchisq
#' @import InteractionSet
#' 
#' @author Altunislav Akalinski
#' 
#' @export
metaAssociations <- function(associations){

  # check the input GRanges structure and if something is wrong return error
  
  if(!checkGrlStrMeta(associations)){
    
    stop("\ncheck if input 'associations' object has the correct structure:\n",
         "\n a)It must be a list of GInteractions,:\n",
         "\n b) and each GInteractions must have:name,name2,n,coefs,pval,pval2,
               qval,qval2 columns\n")
  }

  # create an key-stats table where each key is an association
  # and stats are the association statistics for each data set in the input
  # object
  # for associations, using this table we will create intermediate
  # tables to do meta-analysis
  # this is needed because we can't assume each GRanges object will
  # have the same order or same number of elements, there has to be
  # a merge operation
  
  aTbl = assocTable(associations)

  # create p-values matrix and do meta analysis with Fisher's method
  pCols=seq(2,ncol(aTbl),by=3)+2
  
  p=as.matrix(aTbl[,pCols,with = FALSE]) # this the p-value matrix
  
  comb=-2 * rowSums(log(p),na.rm=TRUE) # this combines p-values
  
  pval=1-pchisq(comb,df=2*rowSums(!is.na(p))) # this calculates new P-values
                                              # using Chi-sq test

  # create sample number and effect size matrices and
  # do weighted averaging for the effect sizes
  coefs=as.matrix(aTbl[,(pCols-1),with = FALSE])
  
  ns=as.matrix(aTbl[,(pCols-2),with = FALSE])
  
  coefs2=rowSums(ns*coefs,na.rm=TRUE)/rowSums(ns,na.rm=TRUE) # average coeffs

  #recreate the GInteractions from the key
  
  df=do.call("rbind",strsplit(aTbl$ky,"||",fixed=T))
  
  gr <- GInteractions(GRanges(df[,1]),
                      GRanges(df[,4]),
                      name=df[,2],
                      name2=df[,3],
                      n=rowSums(ns,na.rm=T),
                      coefs=coefs2,
                      pval=pval,
                      qval=NA)
  
  # calculate q-value

  qval <- try(qvalue::qvalue(gr$pval)$qvalues,silent=T)
  
  if(!inherits(qval, 'try-error')){gr$qval <- qval}
  
  #return the GInteractions 
  gr
}

#----- private functions ------

# this function creates a table of keys and stats from
# a list which contains GInteractions  as the
# output of associateReg2Gene or metaAssociate functions.
# this is needed to process the associations that
# exist in many datasets and merge them into the same table with
# stats from the statistical tests
#' @author Altunislav Akalinski
#' @keywords internal
assocTable<-function(grlist){

  # make keys and stats for each GRanges
  # keys are locations for TSS,regulatory regions and names
  dt.list=lapply(grlist,function(x)

    data.table::data.table(ky=paste(as.character(first(x)),
                                    x$name,
                                    x$name2,
                                    as.character(second(x)),sep="||"),
                                    n=x$n,
                                    coefs=x$coefs,
                                    pval=x$pval,
                                    key="ky")

  )

  # merge tables
  # data.table merge is faster, that's why we rely on it
  Reduce(function(...) merge(..., by = "ky"),
         dt.list)
}



# the function checks if the input GInteractions  in the list have
# the right columns, it should have "name","name2","coefs","pval",...
# columns.
#' @author Altunislav Akalinski
#' @keywords internal
checkGrlStrMeta <- function(associations){

  ExpectedAssResults <- c("name","name2","n","coefs","pval","pval2",
                          "qval","qval2")
 
  all(c(sapply(associations,function(x){
 
       is.reg=all(ExpectedAssResults%in%colnames(mcols(x)))
    
       is.gr=(class(x) == "GInteractions")
    
    all(is.reg,is.gr)}
    
  ), (class(associations)=="list") ))

  }



# not used
# This function maps a vector of p-values to the open unit interval
# (that is, it moves them away from 0 and 1).
# Needed for input to p-value combining function.
# @keywords internal
# open01<-function (p, B)
# {
#   (p + 1/(2 * B))/(1 + 1/B)
# }

# not used here for reference
# combine p-vals using Fisher's method
# @keywords internal
# @example
# fisherComb<-function(p, B)
# {
#   if(missing(B)){
#     comb=-2 * sum(log(p))
#   }else{
#     comb=-2 * sum(log(open01(p, B)))
#   }
#   1-pchisq(comb,df=2*length(p))
# }


