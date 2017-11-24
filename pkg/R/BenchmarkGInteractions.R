#####################################################
#-----------------Benchmark f()---------------------#
#####################################################

#' Benchmarks reg2gene models using benchmark data
#' 
#' The function that takes as input results of \code{\link{associateReg2Gene}} 
#' or any other modelling procedure implemente in \code{reg2gene} package,
#' and predefined benchmark dataset as 
#' \code{\link[InteractionSet]{GInteractions}} object.
#' This function adds a metadata column with info about benchmarking success -
#' whether tested regions are benchmarked or not.  
#' By default it reportes how many times interactions is observed in the 
#' benchmark dataset. If binary is set to TRUE, then vector of 0' and 1's is 
#' reported (1 - overlapping benchmark dataset at least once) and 0 (not
#' overlapping benchmark dataset at all). 
#' 
#' 
#' @param reg2Gene a \code{\link[InteractionSet]{GInteractions}} object output
#' from \code{\link{associateReg2Gene}}). Usually, 1st GRanges object, or anchor1
#' corresponds to the enahncer location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. 
#' @param benchData a \code{\link[InteractionSet]{GInteractions}} object output
#' from \code{\link{associateReg2Gene}}) or a list of GInteractions object. 
#' Both regions are used in the benchmarking procedure. 
#' This object stores benchmarking informations eg interacting region
#' coordinates from techniques such as HiC,eQTL studies...
#' @param binary (def:FALSE) how many times reg2Gene interactions is observed in
#' the benchmark dataset(s). If TRUE, reports if overlap with benchmark dataset 
#' is observed at least once).
#' @param nCores possible to be runned in parallel. Argument for mclapply f();
#' how many cores to use.
#' @param ignore.strand argument to be passed to
#' \code{\link[IRanges]{findOverlaps}}. When set to TRUE, the strand 
#' information is ignored in the overlap analysis.
#' @param ... further arguments to methods, not implemented yet
#' 
#' 
#' @return GInteractions object with added benchmark results metadata [Bench 
#' column].Each column metadata column corresponds to one benchmark dataset 
#' analyzed if input is list() 
#' Values can be either 0/1 (not/benchmarked) or 0-n (how many times each 
#' gene-enhancer pair is benchmared).
#' 
#' 
#' @details GInteractions objects - an output of \code{\link{associateReg2Gene}}
#' [or a list of such objects] and benchmark dataset are overlaped. 
#' \code{\link[InteractionSet]{linkOverlaps}} between reg2Gene 
#' and benchmark object is performed, and for each input pair it is reported 
#' whether this pair is benchmarked or not,and how many times (if binary=F).
#' Criss-cross overlap of interacting regions is performed; thus is anchor1 from
#' benchmark dataset is overlapping anchor2 from tested dataset, than anchor2
#' from benchmark dataset needs to overlap anchor1 from tested datased, or
#' vice-versa.  
#' @author Inga Patarcic
#' 
#' @import GenomicRanges
#' @import InteractionSet
#' 
#' @examples
#' require(GenomicRanges)
#' require(InteractionSet)
#'    
#'    reg2Gene <- GInteractions(GRReg1_toy,GRReg1_toy$reg)
#'    benchData <- GInteractions(GRReg2_toy,GRReg2_toy$reg)
#'    
#' benchmarkGI(reg2Gene,
#'             benchData,
#'             binary=FALSE)
#'             
#' benchmarkGI(reg2Gene,
#'             benchData,
#'             binary=TRUE)             
#' 
#'    
#' ##################   
#' # example for list:
#' 
#' # NOTE: anchor1.Bench1Exp & anchor1.Bench2Exp are expected/precalculated
#' # values for this benchmark example
#'  
#'   benchDataList <- list(benchData,reg2Gene)
#'   names(benchDataList) <- c("benchData1","benchData2")
#' 
#' 
#' reg2GeneB <- benchmarkGI(reg2Gene,
#'                         benchDataList,
#'                         ignore.strand=TRUE,
#'                         binary=FALSE,
#'                         nCores = 1)               
#'             
#'  # Checking what happends when anchor1&anchor2 both overlap only one region
#'  # in benchmark dataset?  OK, They are not benchmarked...
#'  
#'    benchmarkGI(reg2Gene[1],
#'                benchData[1])      
#'
#'  # WARNING!
#'  # However, be careful with benchmarkin anchors that overlap, because they
#'  # will be reported to be benchmarked.  
#'    
#'    benchmarkGI(reg2Gene[5],
#'                benchData)           
#'                                     
#' @export      
benchmarkGI <- function(reg2Gene,
                        benchData,
                        binary=FALSE,
                        nCores=1,
                        ignore.strand=FALSE,
                        ...) {
  
  if (class(benchData)=="GInteractions"){
    
    
    BenchmarkRes <- benchmarkGIsimple(reg2Gene=reg2Gene,
                                      benchData=benchData,
                                      binary=binary,
                                      ignore.strand=ignore.strand)
    
    return(BenchmarkRes)
    
  }
  
  
  
  if (class(benchData)=="list"){
    
    
    BenchmarkRes <- parallel::mclapply(benchData, function(x) {
      
      benchmarkGIsimple(reg2Gene=reg2Gene,
                        benchData=x,
                        binary=binary,
                        ignore.strand=ignore.strand)},
      
      mc.cores = nCores)
    
    # pooling results together
    
    BenchResults <- DataFrame(lapply(BenchmarkRes,function(x) x$Bench))
    
    mcols(reg2Gene) <- c(mcols(reg2Gene),BenchResults)
    
    
    return(reg2Gene)
    
    
    
  }  
  
  
}



#' Benchmarks help function 
#' 
#' Benchmark reg2gene models using benchmark data but for only one GInteraction
#' object [not lists]
#' 
#' Read description for \code{benchmarkGI}
#' 
#' @author Inga Patarcic
#' @keywords internal
benchmarkGIsimple <- function(reg2Gene,
                              benchData,
                              ignore.strand,
                              binary) {
  
  # perform Overlaps 
  OverlapGI <- linkOverlaps(reg2Gene,
                            first(benchData),
                            second(benchData),
                            ignore.strand=ignore.strand)
  
  
  # getting rows that are confirmed by benchmark
  
  benchRow <- OverlapGI$query[OverlapGI$subject1==OverlapGI$subject2]
  
  # reporting benchmark results
  
  reg2Gene$Bench <- rep(0,length(reg2Gene))
  
  if (binary==TRUE){ reg2Gene$Bench[unique(benchRow)] <- 1}
  
  if (binary==FALSE){ 
    # add counts   
    reg2Gene$Bench[as.integer(names(table(benchRow)))] <- table(benchRow)
    
  }
  
  
  return(reg2Gene)
  
}





 ####################################
 # filtering f()

#' Preadjustiong for high number of true negatives in regulatoryReg-TSS pairs
#' 
#' Function that eliminates all regulatory region-TSS pairs [anchor1 and anchor2 
#' from reg2Gene] that do not overlap with any benchmark anchor1 or anchor2 
#' location.
#' Eg it selects reg2gene regions only when regulatory region and TSS overlap 
#' with at least one benchmarking region.  This is useful to improve the 
#' confusion matrix statistics reported by \code{confusionMatrix} by eliminating
#' the high number of true negatives. TN are very abundant in the reg2gene 
#' dataset since benchmark dataset usually covers much smaller regions of
#' the genome (method limitations)
#'
#' @param reg2Gene a \code{\link[InteractionSet]{GInteractions}} object output
#' from \code{\link{associateReg2Gene}}). Usually, 1st GRanges object,or anchor1
#' corresponds to the enahncer location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. 
#' @param benchData a \code{\link[InteractionSet]{GInteractions}} object output
#' from \code{\link{associateReg2Gene}}) or a list of GInteractions object. 
#' Both regions are used in the benchmarking procedure. 
#' This object stores benchmarking informations eg interacting region
#' coordinates from techniques such as HiC,eQTL studies...
#' @param nCores possible to be runned in parallel. Argument for mclapply f();
#' how many cores to use.
#' @param ignore.strand argument to be passed to
#' \code{\link[IRanges]{findOverlaps}}. When set to TRUE, the strand 
#' information is ignored in the overlap analysis.
#' @param ... further arguments to methods, not implemented yet
#' 
#' 
#' @return GInteractions object with added benchmark results metadata [Filter 
#' column].Each column metadata column corresponds to one benchmark dataset 
#' analyzed if input is list() 
#' Values are either 0/1 (possible or not to benchmark).
#' 
#' @details All benchmark regions that can be confirmed by any combination of 
#  enh/gene pairs [anchor1 or anchor2 form reg2Gene object] is obtained. Then
#' selected unique union of anchor1 or anchor2 form reg2Gene object is used
#' as anchor1-anchor2 pairs that can be benchmarked. Reasoning, if present in
#' this set anchor1 or anchor2 regions form reg2Gene object necessarily need to
#' have other member of the pair overlapping somewhere in benchmark dataset.
#' 
#' @examples require(GenomicRanges)
#' require(InteractionSet)
#'    
#'    reg2Gene <- GInteractions(GRReg1_toy,GRReg1_toy$reg)
#'    benchData <- GInteractions(GRReg2_toy,GRReg2_toy$reg)
#'    
#'  filterPreBenchGI(reg2Gene,
#'                   benchData,
#'                   nCores=1)
#'                   
#' ##################   
#' # example for list:
#' 
#' # NOTE: anchor1.Bench1Exp & anchor1.Bench2Exp are expected/precalculated
#' # values for this benchmark example
#'  
#'   benchDataList <- list(benchData,reg2Gene)
#'   names(benchDataList) <- c("fData1","fData2")
#'   reg2GeneBF.list <-filterPreBenchGI(reg2Gene,
#'                                      benchDataList,
#'                                      nCores = 1)
#'  
#'  ####################################
#'  
#'  # Checking what happends when anchor1&anchor2 both overlap only one region
#'  # in benchmark dataset?  OK, They are filtered OUT...
#'  
#'    filterPreBenchGI(reg2Gene[1],
#'                     benchData[1],
#'                     nCores = 1) 
#'     
#'   
#' @export
filterPreBenchGI <- function(reg2Gene,
                             benchData,
                             nCores=1,
                             ignore.strand=FALSE,
                             ...) {
  
  if (class(benchData)=="GInteractions"){
    
    
    filterRes <- filterPreBenchGIsimple(reg2Gene=reg2Gene,
                                        benchData=benchData,
                                        ignore.strand=ignore.strand)
    
    return(filterRes)
    
  }
  
  
  
  if (class(benchData)=="list"){
    
    
    filterRes <- parallel::mclapply(benchData, function(x) {
      
      filterPreBenchGIsimple(reg2Gene=reg2Gene,
                             benchData=x,
                             ignore.strand=ignore.strand)},
      
      mc.cores = nCores)
    
    # pooling results together
    
    filterResults <- DataFrame(lapply(filterRes,function(x) x$Filter))
    
    
    mcols(reg2Gene) <- c(mcols(reg2Gene),filterResults)
    
    
    return(reg2Gene)
    
    
    
  }  
  
  
}

#' Filtering help function 
#' 
#' Filter reg2gene models using benchmark data but for only one GInteraction
#' object [not lists]
#' 
#' Read description for \code{filterPreBenchGI}
#' 
#' @author Inga Patarcic
#' @keywords internal 
filterPreBenchGIsimple <- function(reg2Gene,
                          benchData,
                          ignore.strand=FALSE){


         BenchCovered <- linkOverlaps(benchData,
                                     first(reg2Gene),
                                     second(reg2Gene),
                                     ignore.strand=ignore.strand)
                 
         
         filterRow <- unique(c(BenchCovered$subject1,
                             BenchCovered$subject2))
         
         reg2Gene$Filter <- rep(0,length(reg2Gene))
         
         reg2Gene$Filter[filterRow] <- 1
         
  return(reg2Gene)
         
         
         }




#####################################################
#-----------------ConfusionMatrix f()---------------------#
#####################################################


#' ConfusionTable statistics for the benchmarked reg2Gene object
#' 
#' A function takes as an benchmarked reg2gene object and outputs
#' confusion table statistics.
#' 
#' 
#' @param reg2GeneBench GInteractions object with added benchmark 
#' \code{benchmarkGI} and OPTIONAL filerPreBench \code{filterPreBenchGI} 
#' metadata. However,prior this analysis, it is advised to reduce the number of 
#' true negatives by including only reg2gene entries that could be potentially
#' benchmarked. THUS run \code{filterPreBenchGI} prior 
#' running \code{benchmarkGI})
#' 
#' @param benchCol character (default "Bench"), results of benchmarking 
#' procedure; eg a column name of the metadata which indicates the column where 
#' the result of benchmarking procedure is stored. A vector of 0's and 1's is
#' expected.
#' 
#' @param prefilterCol character (default NULL), results of prefiltering 
#' procedure; eg a column name of the metadata which indicates the column where 
#' the result of prefiltering procedure is stored. A vector of 0's and 1's is
#' expected.
#' 
#' @param thresholdID character (def:NULL) name which indicates a column where
#' statistics of the modelling procedure is stored. This column is filtered such
#' that everything below predefined threshold is considered to be statistically 
#' significant association (set to be equal to 0), whereas everything above that
#' threshold is 0".
#' 
#' @param thresholdValue numeric (def:0.05) A value of a threshold. 
#' 
#' @param statistics character (def "ConfusionMatrix"). Currentlty, 
#' "ConfusionMatrix" or "PPV" are options. 
#' 
#' @return integer vector or a list. If "ConfusionMatrix" is requested then the
#' output is list with following elements:TP,FP,TN,FN,Specificity, Accuracy, PPV,
#' NPV,F1. Otherwise only PPV is reported.
#' 
#' 
#' @details Reports statistics based on the confusion matrix. 
#' There are four different categories in the confusion matrix: 
#' TP = (number of) true positive: reg2gene entry that was reported to be 
#' associated (reported gene-enhancer statistics lower than a predefined
#' threshold) and was benchmarked. 
#' FP = (number of) false positive: reg2gene entry that was reported to be
#' associated (reported gene-enhancer statistics lower than a predefined 
#' threshold) BUT was not overlapped with benchmark dataset
#' FN = (number of) false negative: reg2gene entry that was NOT reported to
#' be associated (reported gene-enhancer statistics NOT lower than a predefined 
#' threshold) BUT is benchmarked
#' TN = (number of) true negative: reg2gene entry 
#' that was NOT reported to be associated (reported gene-enhancer
#' statistics NOT lower than a predefined threshold) AND is NOT benchmarked.
#' If no benchmark and no statistically significant data is entered, then 0 is 
#' reported.
#' 
#' The formulas used in this function are: \deqn{Sensitivity = TP/(TP+FN)} 
#' \deqn{Specificity = TN/(TN+FP)} 
#' \deqn{Accuracy = (TP+TN)/(TP+FN+FP+TN)}
#' \deqn{PPV = TP/(TP+FP)} 
#' \deqn{NPV = TN/(TN+FN)} 
#' \deqn{F1 = (2*TP)/((2*TP)+FP+FN)}
#' 
#' @examples require(GenomicRanges)
#' require(InteractionSet)
#'    
#'    reg2GeneBench <- GInteractions(GRReg1_toy,GRReg1_toy$reg)
#'    reg2GeneBench$PValue <- reg2GeneBench$PValue <- seq(0, 1, 
#'    length.out = length(GRReg1_toy))
#' 
#'
#' confusionMatrix(reg2GeneBench,
#'                 thresholdID = "PValue",
#'                 thresholdValue = 0.05,
#'                 benchCol = "anchor1.Bench1Exp",
#'                 prefilterCol = "anchor1.Filter1Exp",
#'                 statistics = "ConfusionMatrix")
#'                 
#' confusionMatrix(reg2GeneBench,
#'                 thresholdID = "PValue", 
#'                 thresholdValue = 0.05,
#'                 benchCol = "anchor1.Bench1Exp",
#'                 prefilterCol = "anchor1.Filter1Exp",
#'                 statistics = "PPV")
#'                
#'                
#'      
#' @export
confusionMatrix <- function(reg2GeneBench,
                            benchCol=NULL,
                            prefilterCol=NULL,
                            thresholdID=NULL,
                            thresholdValue=0.05,
                            statistics="PPV") {
  
  
           if (is.null(thresholdID)){
              
              stop("Cannot run thresholding, no column defined")
              
              }
  
           if (!any(stringr::str_detect(colnames(mcols(reg2GeneBench)),
                                         thresholdID))){
              
              stop("thresholdID column not identified")
              
           }
  
  
          if (!any(stringr::str_detect(colnames(mcols(reg2GeneBench)),
                                       benchCol))){
            
              stop("benchCol column not detected")
            
          }  
  
  
  
  # filter results of benchmarking procedure by prefiltering column
  
          if (!is.null(prefilterCol)){ 
            
            rows.filt <- as.logical(mcols(reg2GeneBench)[,prefilterCol])
            
            reg2GeneBench <- reg2GeneBench[rows.filt]
          
            }
  
  
  ##################################
  # thresholding reg2GeneBench object
  
  
  predictedEntries <- mcols(reg2GeneBench)[,thresholdID] <= thresholdValue
  
  benchmarkedEntries <- as.logical(mcols(reg2GeneBench)[,benchCol])
  
  
  
              if (all(c(predictedEntries==F,benchmarkedEntries==F))){return(0)}
  
  
  
  
  # calculating TP,FN,TN,FP adjusted for 0 observance
            TP <- sum(which(predictedEntries)%in%which(benchmarkedEntries))
            TN <- sum(which(!predictedEntries)%in%which(!benchmarkedEntries))
            FP <- sum(which(predictedEntries)%in%which(!benchmarkedEntries))
            FN <- sum(which(!predictedEntries)%in%which(benchmarkedEntries))
            
  
  
              if (statistics=="PPV"){ 
                
                return(TP/(TP+FP))
                
                }
  
              if (statistics=="ConfusionMatrix"){ 
    
                      statisticsCM <- list()
                      
                      statisticsCM$TP <- TP;  statisticsCM$FP <- FP
                      statisticsCM$TN <- TN;  statisticsCM$FN <- FN
                      statisticsCM$PPV <- TP/(TP+FP)
                      statisticsCM$NPV <- TN/(TN+FN)
                      statisticsCM$Accuracy <- (TP+TN)/(TP+FN+FP+TN)
                      statisticsCM$F1 <- (2*TP)/((2*TP)+FP+FN)
    
                return(statisticsCM)
  
                  
                }
  
  
}

 




