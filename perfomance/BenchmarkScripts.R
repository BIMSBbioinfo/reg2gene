.ReArrangeGR <- function(GR){
  
  # function that rearranges GRanges object as DF if necesarry
  DF <- as.data.frame(GR)
  DF <- rename(DF, c("seqnames"="chr2", "start"="start2","end"="end2",
                     "reg.seqnames"="chr1","reg.start"="start1","reg.end"="end1"))
  Min_Column_names <- c("chr1","start1","end1","chr2","start2","end2")
  # column name rearrangemet
  DF <- (DF[,c(Min_Column_names,colnames(DF)[!colnames(DF)%in%Min_Column_names])])
  
  return(DF)
}


# 1) TO DO: allow GRanges as an input for both parameters
# 2) rewrite foo grnages object
# 3) add filtering based on gene names
# 4) crate small examples, check all functions




#' Benchmarks Reg2gene models using Benchmark data
#' 
#' The function that takes as an input results of \code{linkReg2Gene} or
#' any other modelling procedure implemente in \code{reg2gene} package,
#' such as \code{MetaReg2Gene} and predefined Benchmark dataset.
#' For input Reg2Gene object it adds a column with info whether reported interactions is 
#' benchmarked or not. By default it reportes how many times Reg2Gene interactions is
#' observed in the benchmark dataset. If binary is set to TRUE, then only logical vector of
#' TRUE (overlapping benchmark dataset at least once) and FALSE(not overlapping benchmark 
#' dataset at all) entries is reported.
#' 
#' 
#'  @param Reg2Gene a GRanges object output from \code{linkReg2Gene} or \code{MetaReg2Gene}
#'  function or a character matrix containing at least six columns named:
#'  "chr1","start1","end1","chr2","start2","end2" and where "chr1","start1","end1" corresponds
#'  to regulatory regions whereas "chr2","start2","end2" correpond to TSS of a gene.
#'  
#'  @param Benchmark character matrix with at least three columns named:
#'  "chr1","start1","end1" which store information about coordinates
#'  of the regions which we want to compare and 1) one additional column 
#'  with gene names, or 2)  "chr2","start2","end2" columns
#'   
#' 
#'  @param binary (def:FALSE) how many times Reg2Gene interactions is
#'  observed in the benchmark dataset. If TRUE, reports if overlap with
#'  benchmark dataset is observed at least once).
#' 
#' @param byReg2Gene character (def:NULL), name of the column where gene names are stored 
#' in Reg2Gene object. Optional and useful if filtering should be performed based on the overlap
#' between regions and gene names. 
#' @param byBenchmark character (def:NULL), name of the column where gene names are stored 
#' in Benchmark object. Optional and useful if filtering should be performed based on the overlap
#' between regions and gene names.
#' 
#' @return dataframe equal to the input dataframe, but a column with info which reg2gene entry is
#' benchmarked or not is added.
#' 
#' @import GenomicRanges
#' @import stringr
#' @import plyr
#' 
#' 
#' @details GRanges object or dataframe as an output of linkReg2Gene or Reg2GeneMA 
#' and Benchmark dataset is inputed. Threshold can be entered. First, Reg2Gene object
#' is thresholded based on selected column and threshold value, then ComplexOverlap between
#' Reg2Gene and Benchmark object is performed
#' 
#' 
#' 
#' @examples 
#'  
#'  @export


  BenchMarkReg2Gene <- function(Reg2Gene,
                               Benchmark,
                               binary=FALSE,
                               byReg2Gene=NULL,
                               byBenchmark=NULL) {
   
   # rearrangement into DF
                       if (class(Reg2Gene)=="GRanges") { Reg2Gene <- .ReArrangeGR(GR=Reg2Gene)}
                   # test if input is DF   
                       if (class(Reg2Gene)!="data.frame") {stop("Reg2Gene input should be either dataframe or GRanges object")}
                   # test column names reported
                   Min_Column_names <- c("chr1","start1","end1","chr2","start2","end2")
                     if (min(Min_Column_names%in%colnames(Reg2Gene))==0) {stop("Reg2Gene colnames should be chr1,start1,end1,
                                                                           chr2,start2,end2,columnTOthreshold at least!")}
                   
                   
   
  # Benchmarking           

        # ComplexBenchmarking   
           if (min(Min_Column_names%in%colnames(Benchmark))!=0) {
                  Reg2GeneBenchOverlap <- ComplexOverlaps(Reg1=Reg2Gene,
                                                         Reg2=Benchmark)}
              
        # Region Benchmark + column 
           if (min(Min_Column_names%in%colnames(Benchmark))==0) {
             Reg2GeneBenchOverlap <- OverlapRegions(Reg1=Reg2Gene,
                                                    Reg2=Benchmark,
                                                    byReg1 = byReg2Gene,
                                                    byReg2 = byBenchmark)}
                                                   
 
        # counting benchmark overlaps
             Reg2Gene.Coord <- c("Reg1_chr1","Reg1_start1","Reg1_end1","Reg1_chr2","Reg1_start2","Reg1_end2")
             Reg2GeneBenchmarked.FREQ <- count(Reg2GeneBenchOverlap, vars = Reg2Gene.Coord)
             colnames(Reg2GeneBenchmarked.FREQ)[1:6] <- str_replace(Reg2Gene.Coord,"Reg1_","")
             
             # adding column with counted overlaps          
             Reg2GeneBenchmarked <- join(Reg2Gene,Reg2GeneBenchmarked.FREQ)
                  Reg2GeneBenchmarked$freq[is.na(Reg2GeneBenchmarked$freq)] <- 0
          # counts or binary reported   
            if (binary==TRUE){ Reg2GeneBenchmarked$freq <- as.logical(Reg2GeneBenchmarked$freq)}
       
  return(Reg2GeneBenchmarked)
        
}
  
  
       
  
  
  
         
#' ConfusionTable statistics of Benchmarked Reg2Gene Object
#' 
#' A function takes as an input results of benchmarking procedure for 
#' Reg2Gene objects \code{BenchMarkReg2Gene} and outputs
#' confusion table statistics.
#' 
#' 
#' @param BenchMarkedReg2Gene A result of \code{BenchMarkReg2Gene}.
#' !WARNING! It is advised to reduce the number of true negatives by including only
#' Reg2gene entries that could be potentially benchmarked. THUS run \code{Filter_PreBench}
#' prior running \code{BenchMarkReg2Gene})
#' 
#' @param thresholdID character (def:NULL) name which indicates a column which 
#' will be filtered
#' 
#' @param thresholdValue numeric (def:0.05) A value of a threshold. Everything below 
#' that thresold will be consider as a positive, an is set to be equal to 1.
#' Otherwise 0.
#' 
#' @param statistics character (def "ConfusionMatrix"). Currentlty, "ConfusionMatrix" or
#' "PPV" are options. If "ConfusionMatrix" is requested then output of 
#' \code{\link[caret]{confusionMatrix}} is reported. Otherwise, if "PPV" 
#' is requested then output of \code{\link[caret]{posPredValue}}
#' 
#' @return Reports either positive predictive value or statictics based on 
#' confusion matrix \code{\link[caret]{confusionMatrix}}
#' 
#' @import GenomicRanges
#' @import stringr
#' @import caret::confusionMatrix
#' @import caret::posPredValue
#' 
#' @details Reports statistics based on the confusion matrix. There are four different categories
#' in the confusion matrix: 
#' 1) TP (true positive): reg2gene entry that was reported to be associated (reported gene-enhancer
#' statistics lower than a predefined threshold) and is benchmarked.
#' 2) FP (false positive): reg2gene entry that was reported to be associated (reported gene-enhancer
#' statistics lower than a predefined threshold) BUT does not overlap benchmark dataset
#' 3) TN (true negative): reg2gene entry that was NOT reported to be associated (reported gene-enhancer
#' statistics NOT lower than a predefined threshold) AND is NOT benchmarked.
#' 4) FN (false negative): reg2gene entry that was NOT reported to be associated (reported gene-enhancer
#' statistics NOT lower than a predefined threshold) BUT is benchmarked.
#' 
#' @examples 
#'      
#' @export



ConfusionMatrixReg2Gene <- function(BenchMarkedReg2Gene,
                                    thresholdID=NULL,
                                    thresholdValue=0.05,
                                    statistics="PPV") {
  
   # thresholding BenchMarkedReg2Gene object
         
         if (is.null(thresholdID)){stop("Cannot run thresholding, no column defined")}
         
            PredictedEntries <- BenchMarkedReg2Gene[,thresholdID] <= thresholdValue
            BenchMarkedEntries <- BenchMarkedReg2Gene$freq
            
            Confusion.matrix <- table(factor(PredictedEntries,levels = c(T,F)),
                                      factor(BenchMarkedEntries,levels = c(T,F)))
          
        # carefull ordering
            Confusion.matrix <- Confusion.matrix[c("TRUE","FALSE"),c("TRUE","FALSE")]
            
            if (statistics=="PPV"){ CM.Statistics <-  posPredValue(Confusion.matrix)}
            if (statistics=="ConfusionMatrix"){ CM.Statistics <- confusionMatrix(Confusion.matrix)}
            
           return(CM.Statistics)
            
            
}

