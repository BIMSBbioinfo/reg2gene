# functions for meta analysis

#' meta-analysis for regulatory region and gene associations
#'
#' The function combines association P-values and coefficients from
#' different data
#' sets using Fisher's method and weigthed averaging respectively. After,
#' combining P-values, q-values are calculated using the
#' \code{\link[qvalue]{qvalue}} function.
#'
#' @param associations A GRangesList where each element is a GRanges
#' object output from \code{\link{associateReg2Gene}} function.
#'
#'
#' @return returns a GRanges object with updated statistics
#' from the meta-analysis. The output will be similar to
#' the output of \code{\link{associateReg2Gene}} function, it
#' will have the same meta-data columns in the same order.
#'
#' @examples
#' \dontrun{
#' #
#' assoc=system.file("extdata", "sampleAssociationsGrlist.rds",
#'                       package = "reg2gene")
#'
#' meta=metaAssociations(assoc)
#' }
#'
#' @importFrom  data.table data.table
#'
#' @export
#'
#' @author Altunislav Akalinski
metaAssociations<-function(associations){

  # check the input GRanges structure and if something is wrong return error
  if(!checkGrlStrMeta(associations)){
    stop("\ncheck if input 'associations' object has the correct structure:\n",
         "1) It must be a GRangesList,
          2) where each GRanges has a column called 'reg'\n",
         "and that column points to another GRanges\n",
         "3) GRanges objects must have 'coefs' and 'pval' columns")
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
  p    =as.matrix(aTbl[,pCols,with = FALSE]) # this the p-value matrix
  comb=-2 * rowSums(log(p),na.rm=TRUE) # this combines p-values
  pval=1-pchisq(comb,df=2*rowSums(!is.na(p))) # this calculates new P-values
                                              # using Chi-sq test

  # create sample number and effect size matrices and
  # do weighted averaging for the effect sizes
  coefs=as.matrix(aTbl[,(pCols-1),with = FALSE])
  ns=as.matrix(aTbl[,(pCols-2),with = FALSE])
  coefs2=rowSums(ns*coefs,na.rm=TRUE)/rowSums(ns,na.rm=TRUE) # average coeffs

  #recreate the GRanges from the key
  df=do.call("rbind",strsplit(aTbl$ky,"||",fixed=T))
  gr=GRanges(df[,1],name=df[,2],name2=df[,3],reg=GRanges(df[,4]),
             n=rowSums(ns,na.rm=T),coefs=coefs2,pval=pval)

  # calculate q-value
  gr$qval=qvalue(gr$pval)$qvalues

  #return the GRanges
  gr
}

#----- private functions ------

# this function creates a table of keys and stats from
# a GRangesList which contains GRanges as the
# output of associateReg2Gene or metaAssociate functions.
# this is needed to process the associations that
# exist in many datasets and merge them into the same table with
# stats from the statistical tests
assocTable<-function(grlist){

  # make keys and stats for each GRanges
  # keys are locations for TSS,regulatory regions and names
  dt.list=lapply(grlist,function(x)

    data.table(ky=paste(as.character(x),x$name,x$name2,
                        as.character(x$reg),sep="||"),
               n=x$n,coefs=x$coefs,pval=x$pval,key="ky")

  )

  # merge tables
  # data.table merge is faster, that's why we rely on it
  Reduce(function(...) merge(..., by = "ky"),
         dt.list)
}



# the function checks if the input GRanges in the GRangesList have
# the right columns, it should have"name","name2","reg","coefs","pval"
# columns. And reg column points to a GRanges object
# returns TRUE if all checks out
checkGrlStrMeta<-function(associations){

  all(c(sapply(associations,function(x){
    is.reg=all(c("name","name2","reg","coefs","pval")  %in% colnames(mcols(x)))
    #is.reg=sum(str_detect(colnames(mcols(x)),"name|name2|reg|coefs|pval"))==5
    is.gr=(class(x$reg) == "GRanges")
    all(is.reg,is.gr)}
  ), (class(associations)=="GRangesList") ))
}



# not used
# This function maps a vector of p-values to the open unit interval
# (that is, it moves them away from 0 and 1).
# Needed for input to p-value combining function.
# @keywords internal
open01<-function (p, B)
{
  (p + 1/(2 * B))/(1 + 1/B)
}

# not used here for reference
# combine p-vals using Fisher's method
# @keywords internal
# @example
fisherComb<-function(p, B)
{
  if(missing(B)){
    comb=-2 * sum(log(p))
  }else{
    comb=-2 * sum(log(open01(p, B)))
  }
  1-pchisq(comb,df=2*length(p))
}


