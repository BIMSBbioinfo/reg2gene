# voteAssociations()

#' Majority vote decision for (regulatory region)-gene associations
#'
#' The function calculates regularatory region to gene associations based on
#' a majority vote. Multiple associations output by different methods and
#' data sets can be merged this way, and only the associations that have
#' support in \code{vote.threshold} fraction of the datasets will be retained.
#'
#' @param associations A list of \code{\link[InteractionSet]{GInteractions}}
#'  objects outputed from \code{\link{associateReg2Gene}} or
#' \code{\link{metaAssociations}}.
#' @param cutoff.stat (character,"pval" is default). Which statistics to 
#' filter:"qval" or "pval"
#' @param cutoff.val q-value cutoff that will be used to filter elements in the
#' input list. This argument will only take effect when the 
#' \code{\link[InteractionSet]{GInteractions}} object has a numeric "qval" 
#' column. If the input object lacks this column, every association in the 
#' object will be treated as a valid association prediction.
#' @param vote.threshold A value between 0 and 1, designates the threshold
#' needed for fraction of votes necessary to retain an association. Defaults to
#' 0.5, meaning fraction of votes should be  greater than or equal to 0.5 to
#' retain association.
#' 
#' @details Firstly, function selects POSITIVES (statistically associated 
#' gene~enhancer pairs) for each result of \code{\link{associateReg2Gene}} 
#' analysis that wants to be combined by  majority voting (for example results 
#' of H3K4me1 and H3K27ac). Assessing statistically associated gene~enhancer 
#' pairs has been done by filtering the statistics (cutoff.stat) of the elements
#' of the input list (gene~enhancer pairs) based on the defined cutoff value 
#' (cutoff.val). 
#'
#' @return A \code{\link[InteractionSet]{GInteractions}} object that contains 
#' votes for associations from the voting procedure. The object will contain 
#' meta-columns for individual votes and vote stastics. The object can be 
#' further filtered to obtain the desired level of votes using \code{[]} or 
#' equivalent methods.
#' @importFrom data.table data.table
#' @import InteractionSet
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
#'  # Run voteAssociations
#'  voteAssociations(associations, 
#'                   cutoff.stat="pval",
#'                   cutoff.val=0.05,
#'                   vote.threshold=0.5)
#'                   
#'  voteAssociations(associations,
#'                   cutoff.stat="pval",
#'                   cutoff.val=0.05,
#'                   vote.threshold=0.51)
#' @author Altuna Akalin
#' @export
voteAssociations<-function(associations,
                           cutoff.stat="pval",
                           cutoff.val=0.05,
                           vote.threshold=0.5){

  # check the input GRanges structure and if something is wrong return error
  if(!checkGrlStr(associations)){
    
    stop("\ncheck if input 'associations' object has the correct structure:\n",
         "\n a)It must be a list of GInteractions,:\n",
         "\n b) and each GInteractions must have:name,name2,n,coefs,pval,pval2,
               qval,qval2 columns\n")
  }
  
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

  gr <- GInteractions(GRanges(df[,1]),
                      GRanges(df[,4]),
                      name=df[,2],
                      name2=df[,3],
                      votes=votes)
  
  
  # arrange GRanges and return
  
  gr=gr[gr$votes>=(vote.threshold*ncol(q)),]
  
  return(gr)
}





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
#' @importFrom data.table data.table
#' 
#' @keywords internal
#' @author Altuna Akalin
assocTable2<-function(grlist,
                      qval,
                      cutoff.column){

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
              data.table::data.table(ky=paste(as.character(first(x)),
                                              x$name,
                                              x$name2,
                                              as.character(second(x)),sep="||"),
                       association=1,key="ky")
              }
            },qval=qval)

  
  # merge tables
  # data.table merge is faster, that's why we rely on it
  Reduce(function(...) merge(...,all=TRUE, by = "ky"),
         dt.list)
}

# the function checks if the input GInteractions  in the list have
# the right columns, it should have "name","name2","coefs","pval",...
# columns.
#' @author Altunislav Akalinski
#' @keywords internal
checkGrlStr<-function(associations){
  
  ExpectedAssResults <- c("name","name2","n","coefs","pval","pval2",
                          "qval","qval2")
  
  all(c(sapply(associations,function(x){
    
    is.reg=all(ExpectedAssResults%in%colnames(mcols(x)))
    
    is.gr=(class(x) == "GInteractions")
    
    all(is.reg,is.gr)}
    
  ), (class(associations)=="list") ))
  
}

