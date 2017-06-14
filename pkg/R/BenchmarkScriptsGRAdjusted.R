#' Benchmarks Reg2gene models using Benchmark data
#' 
# #' The function that takes as an input results of \code{linkReg2Gene} or
#' any other modelling procedure implemente in \code{reg2gene} package,
# #' such as \code{MetaReg2Gene} and predefined Benchmark dataset.
#' For input Reg2Gene object it adds a column with info whether reported interactions is 
#' benchmarked or not. By default it reportes how many times Reg2Gene interactions is
#' observed in the benchmark dataset. If binary is set to TRUE, then only logical vector of
#' TRUE (overlapping benchmark dataset at least once) and FALSE(not overlapping benchmark 
#' dataset at all) entries is reported.
#' 
#' 
#' @param Reg2Gene a GRanges object output from \code{linkReg2Gene} or \code{MetaReg2Gene}
#' @param Benchmark a GRanges object
#'   
#' 
#' @param binary (def:FALSE) how many times Reg2Gene interactions is
#'  observed in the benchmark dataset. If TRUE, reports if overlap with
#'  benchmark dataset is observed at least once).
#' 
#' 
#' @return a GRanges object which corresponds to Reg2Gene object
#' but with column added which reports how many times interaction was overlapped
#' or whether reported interaction was overlapped
#' 
#' 
#' 
#' @details GRanges object or dataframe as an output of linkReg2Gene or Reg2GeneMA 
#' and Benchmark dataset is inputed. ComplexOverlap between Reg2Gene and Benchmark 
#' object is performed, and for each input pair it is reported whether this pair is
#' benchmarked or not,and how many times (if binary=F) 
#' 
#' 
#' @examples BenchMarkReg2Gene(GRReg1_toy ,GRReg2_toy )
#' BenchMarkReg2Gene(GRReg2_toy,GRReg1_toy )
#' 
#'  
#'@export
BenchMarkReg2Gene <- function(Reg2Gene,
                               Benchmark,
                               binary=FALSE){
                   
  
  
  ############################### 
  # Benchmarking           

        # ComplexBenchmarking   
      Reg2GeneBenchOverlap <- ComplexOverlaps(Reg2Gene,Benchmark)
 
          # counts or binary reported   
            OverlapVector <- rep(0,length(Reg2Gene))
            DuplicatedBenchmarks <- table(Reg2GeneBenchOverlap$Coor1Coord2PAIR)
            OverlapVector[as.integer(names(DuplicatedBenchmarks))] <- DuplicatedBenchmarks
          
        if (binary==TRUE){ OverlapVector <- as.logical(OverlapVector)}
           
            Reg2Gene$BenchmarkO <- OverlapVector
            
            
  return(Reg2Gene)
        
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
#' @examples ConfusionMatrixReg2Gene(BenchMarkedReg2Gene=BenchmarkedReg2Gene_toy,thresholdID = "PValues",thresholdValue = 0.05)
#' 
#' ConfusionMatrixReg2Gene(BenchmarkedReg2Gene_toy,thresholdID = "PValues",thresholdValue = 0.3, statistics = "ConfusionMatrix")
#' 
#' 
#'      
#' @export
ConfusionMatrixReg2Gene <- function(BenchMarkedReg2Gene,
                                    thresholdID=NULL,
                                    thresholdValue=0.05,
                                    statistics="PPV") {
  
  require(caret)
   # thresholding BenchMarkedReg2Gene object
         
         if (is.null(thresholdID)){stop("Cannot run thresholding, no column defined")}
         
            PredictedEntries <- mcols(BenchMarkedReg2Gene)[,thresholdID] <= thresholdValue
            BenchMarkedEntries <- mcols(BenchMarkedReg2Gene)[,"BenchmarkO"]
            
            Confusion.matrix <- table(factor(PredictedEntries,levels = c(T,F)),
                                      factor(BenchMarkedEntries,levels = c(T,F)))
          
        # carefull ordering
            Confusion.matrix <- Confusion.matrix[c("TRUE","FALSE"),c("TRUE","FALSE")]
            
            if (statistics=="PPV"){ CM.Statistics <-  posPredValue(Confusion.matrix)}
            if (statistics=="ConfusionMatrix"){ CM.Statistics <- confusionMatrix(Confusion.matrix)}
            
           return(CM.Statistics)
            
            
}





