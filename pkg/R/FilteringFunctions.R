.Exp_filtering = function(INPUT){
  
  if (any(INPUT!=0) &
      
      (sum(INPUT==0,na.rm=T)/length(INPUT) < 0.9) &
      
      (var(as.numeric(INPUT),na.rm=T)!=0))
    
    #&(as.logical(min(complete.cases(INPUT)))))
    
  {return(TRUE)}else{
    return(FALSE)
  }
  
}


#' Filters genes and enhancer regions prior modelling
#' 
#' The function that takes as an input individual GRanges object stored in GRangesList from
#' \code{regActivityAroundTSS} and returns a vector of TRUE,FALSE entries corresponding to
#' members of GRange object and whether they pass filtering or not.  
#' 
#' @param TSSmodel a GRanges object, member of the GRangesList from \code{regActivityAroundTSS}
#'   
#' @return a vector of TRUE,FALSE entries which indicates whether ranges in the input
#' GRanges object passed filtering or not.
#'  
#' @import GenomicRanges
#'  
#' @details It checks whether genes and enhancer region have:a) an expression 
#' or enhancer activity eqaul to 0 across all samples(cell types) or 90% of expression values or
#' enhancer activity is equal to 0; b) variance across samples is different from 0.
#' FALSE value indicates that these ranges should filtered prior the modelling procedure; 
#' TRUE indicates that it is ok to proceed with this ranges. 
#' !IMPORTANT NOTE! If the first entry is FALSE, than that gene should be removed from the analysis
#' 
#' 
#' @example 
#' 
#' GR_exp_toy <- GRanges(rep("chr1",4),IRanges(1:4,2:5))
#' mcols(GR_exp_toy) <-  matrix(rep("test",16),nrow=4,
#'                              dimnames = list(1:4,c("gene.indicator","featureType","name","name2")))
#' mcols(GR_exp_toy) <- cbind(mcols(GR_exp_toy),data.frame(matrix(c(rep(0,10),c(rep(0,9),1),
#' c(rep(0,3),rep(1,7)),c(rep(1,9),NA)),byrow = T,nrow=4),stringsAsFactors = F))
#' FilterPerGeneModel(GR_exp_toy)
#' 
#' 
#' @export
FilterPerGeneModel <- function(TSSmodel){
  
  MetaData <- data.frame(mcols(TSSmodel),stringsAsFactors = F)
  
  MetaData <- MetaData[,!colnames(MetaData)%in%c("gene.indicator","featureType","name","name2")]
  
  return(apply(MetaData,1,.Exp_filtering))
  
}


















#' #' Selects Reg2gene entries which could be benchmarked with Benchmark dataset
#' #'
#' # #' The function that takes as an input results of \code{linkReg2Gene} or
#' #' any other modelling procedure implemente in \code{reg2gene} package,
#' # #' such as \code{MetaReg2Gene} and predefined Benchmark dataset.
#' #' For input Reg2Gene object it adds a column with info whether reported interactions is
#' #' benchmarked or not. By default it reportes how many times Reg2Gene interactions is
#' #' observed in the benchmark dataset. If binary is set to TRUE, then only logical vector of
#' #' TRUE (overlapping benchmark dataset at least once) and FALSE(not overlapping benchmark
#' #' dataset at all) entries is reported.
#' #'
#' #'
#' # #'  @param Reg2GenePreFilter a GRanges object output from \code{linkReg2Gene} or \code{MetaReg2Gene}
#' #'  function or a character matrix containing at least six columns named:
#' #'  "chr1","start1","end1","chr2","start2","end2" and where "chr1","start1","end1" corresponds
#' #'  to regulatory regions whereas "chr2","start2","end2" correpond to TSS of a gene.
#' #'
#' #'  @param Benchmark character matrix with at least three columns named:
#' #'  "chr1","start1","end1" which store information about coordinates
#' #'  of the regions which we want to compare and 1) one additional column
#' #'  with gene names, or 2)  "chr2","start2","end2" columns
#' #'
#' #'
#' #' @param byReg2Gene character (def:NULL), name of the column where gene names are stored
#' #' in Reg2Gene object. Optional and useful if filtering should be performed based on the overlap
#' #' between regions and gene names.
#' #' @param byBenchmark character (def:NULL), name of the column where gene names are stored
#' #' in Benchmark object. Optional and useful if filtering should be performed based on the overlap
#' #' between regions and gene names.
#' 
#' #library(plyr)
#' #library(stringr)
#' #library(GenomicRanges)
#' #library(reg2gene)
#' 
#' 
#' Filter_PreBench <- function(Reg2GenePreFilter,
#'                             Benchmark,
#'                             byReg2Gene=NULL,
#'                             byBenchmark=NULL){
#' 
#'       if (class(Reg2GenePreFilter)=="GRanges") { Reg2GenePreFilter <- .ReArrangeGR(GR=Reg2GenePreFilter)}
#' 
#' 
#'       # setting both Benchmark coordinates on top of one another
#' 
#'          MinNames <- c("chr1","start1","end1")
#'               Benchmark.AllRegions <- rbind(Benchmark[,MinNames],setNames(Benchmark[,str_replace(MinNames,"1","2")],MinNames))
#'           # filtering overlap of Coord1
#'               Reg2GenePreFilter.Coord1 <- OverlapRegions(Reg2GenePreFilter,Benchmark.AllRegions)
#' 
#'           # filtering overlap of Coord2, after filtering of overlap Coord1 already performed
#'                   Reg2GenePreFilter.Coord1 <- rename(Reg2GenePreFilter.Coord1,"Reg1_chr2"="chr1","Reg1_start2"="start1","Reg1_end2"="end1")
#'                   Reg2GenePreFilter.Coord2 <- OverlapRegions(Reg2GenePreFilter.Coord1,Benchmark.AllRegions)
#' 
#'           #rearranging the format
#'                     Reg2GeneFiltered <- Reg2GenePreFilter.Coord2[,1:ncol(Reg2GenePreFilter)]
#'                           colnames(Reg2GeneFiltered) <- colnames(Reg2GenePreFilter)
#'                     Reg2GeneFiltered <- unique(Reg2GeneFiltered)
#' 
#'       return(Reg2GeneFiltered)
#' 
#' 
#'           }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #Filter_PreMetA <- function()
