# voteAssociations()

#' Majority vote decision for (regulatory region)-gene associations
#'
#' The function calculates regularatory region to gene associations based on
#' a majority vote. Multiple associations output by different methods and
#' data sets can be merged this way, and only the associations that have
#' support in \code{vote.threshold} fraction of the datasets will be retained.
#'
#' @param associations a named \code{GRangesList} where each element is a set of
#' associations returned by \code{\link{associateReg2Gene}} or
#' \code{\link{metaAssociations}}.
#' @param cutoff.stat (character,"pval" is default). 
#' Which statistics to filter:"qval" or "pval"
#' @param cutoff.val q-value cutoff that will be used to filter elements in the
#' input \code{GRangesList}. This argument will only take effect when the
#' \code{GRanges} object has a numeric "qval" column. If the input object
#' lacks this column, every association in the object will be treated as a valid
#' association prediction.
#' @param vote.threshold A value between 0 and 1, designates the threshold
#' needed for fraction of votes necessary to retain an association. Defaults to
#' 0.5, meaning fraction of votes should be  greater than or equal to 0.5 to
#' retain association.
#'
#' @return A \code{GRanges} object that contains votes for associations from
#' the voting procedure. The object will contain meta-columns for individual
#' votes and vote stastics. The object can be further filtered to obtain
#' the desired level of votes using \code{[]} or equivalent methods.
#' @import data.table data.table
#' @examples
#' \dontrun{
#' #
#' assoc=system.file("extdata", "sampleAssociationsGrlist.rds",
#'                       package = "reg2gene")
#'
#' voted=voteAssociations(assoc, cutoff.val=0.05,vote.threshold=0.5)
#'
#' # get all the vote counts
#' voted=voteAssociations(assoc, cutoff.val=0.05,vote.threshold=0)
#' }
#' @author Altuna Akalin
#' @export
voteAssociations<-function(associations,
                           cutoff.stat="pval",
                           cutoff.val=0.05,
                           vote.threshold=0.5){

  # check the input GRanges structure and if something is wrong return error
  if(!checkGrlStr(associations)){
    stop("\ncheck if input 'associations' object has the correct structure:\n",
         "it must be a GRangesList where each GRanges has a column called 'reg'\n",
         "and that column points to another GRanges")
  }
  
  # function that removes tags from column names when necessarry
  associations <- checkTagNames(associations=associations,
                                cutoff.column=cutoff.stat)

  # create Table from after merge
  dtl=assocTable2(grlist=associations,
                  cutoff.column=cutoff.stat,
                  qval=cutoff.val)

  q=as.matrix(dtl[,-1,with = FALSE]) # make presence/absence matrix for assoc.
  q[is.na(q)]=0

  # calculate votes
  votes=rowSums(q,na.rm=TRUE)

  #recreate the GRanges from the key
  df=do.call("rbind",strsplit(dtl$ky,"||",fixed=T))
  gr=GRanges(df[,1],name=df[,2],name2=df[,3],reg=GRanges(df[,4]),
             votes=votes)


  # arrange GRanges and return
  
  gr=gr[gr$votes>=(vote.threshold*ncol(q)),]
  return(gr)
}



#voteAssociations(associations=a, cutoff.stat="pval",cutoff.val=0.05,vote.threshold=0.5)
#----- private functions ------

#' this function creates a table of keys and pvals from
#' a GRangesList which contains GRanges as the
#' output of associateReg2Gene or metaAssociate functions.
#' this is needed to process the associations that
#' exist in many datasets and merge them into the same table with
#' stats from the statistical tests
#' 
#' @param grlist 
#' @param cutoff.column (character,"pval" is defalut). 
#' Which statistics to filter:"qval" or "pval"
#' @param qval q-value cutoff that will be used to filter elements in the
#' input \code{GRangesList}. This argument will only take effect when the
#' \code{GRanges} object has a numeric "qval" column. If the input object
#' lacks this column, every association in the object will be treated as a valid
#' association prediction.
#' 
#' @keywords internal
#' @author Altuna Akalin
assocTable2<-function(grlist,
                      qval,
                      cutoff.column){

require(data.table)  
  # make keys and stats for each GRanges
  # keys are locations for TSS,regulatory regions and names
  dt.list=lapply(grlist,function(x,qval){

    # adjusting for filtering for pval or qval
    if (cutoff.column=="qval"){
    
            if("qval" %in% colnames(mcols(x))){
              x$qval[is.na(x$qval)]=1#remove NAs
              x=x[x$qval<qval,]
            }}
    
    if (cutoff.column=="pval"){
      
      if("pval" %in% colnames(mcols(x))){
        x$pval[is.na(x$pval)]=1#remove NAs
        x=x[x$pval<qval,]
      }}
    
            
            # adjusting for the fact that no qvalue or pvalue pass threshold
            if (length(x)==0){stop("No association identified!")}
            
            if (length(x)!=0){
            data.table(ky=paste(as.character(x),x$name,x$name2,
                                as.character(x$reg),sep="||"),
                       association=1,key="ky")
              }
            },qval=qval)

  
  # merge tables
  # data.table merge is faster, that's why we rely on it
  Reduce(function(...) merge(...,all=TRUE, by = "ky"),
         dt.list)
}

#' the function checks if the input GRanges in the GRangesList have
#' the right columns, it should have at least a reg column as a GRanges object
#' returns TRUE if GRanges has a reg column pointing to a GRanges and
#' if associations is a GRangesList object
#' @param associations a named \code{GRangesList} where each element is a set of
#' associations returned by \code{\link{associateReg2Gene}} or
#' \code{\link{metaAssociations}}.
#' @keywords internal
#' @author Altuna Akalin
checkGrlStr<-function(associations){

  all(c(sapply(associations,function(x){
  is.reg=("reg" %in% colnames(mcols(x)) )
  is.gr=(class(x$reg) == "GRanges")
  all(is.reg,is.gr)}
  ), (class(associations)=="GRangesList") ))

  

  }


#' f() that removes tags if added to qval or pval columns
#' 
#' @param associations a named \code{GRangesList} where each element is a set of
#' associations returned by \code{\link{associateReg2Gene}} or
#' \code{\link{metaAssociations}}.
#' @param cutoff.column (character,"pval" is defalut). 
#' Which statistics to filter:"qval" or "pval"
#' 
#' @keywords internal
#' @author Inga Patarcic
checkTagNames <- function(associations,cutoff.column){

  associations <- lapply(associations,function(x){
  
        tag.added <- names(mcols(x))[stringr::str_detect(names(mcols(x)),
                                                         cutoff.column)]
        
        tag.added <- str_replace(tag.added,cutoff.column,"")
        
        if  (tag.added!="") {names(mcols(x)) <- str_replace(names(mcols(x)),
                                                            tag.added,"")}
        
        return(x)
        
      })

return(associations)
  
}
