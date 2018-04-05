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
#' @param forceByName (def:FALSE) force benchmark data to have an equal gene 
#' coordinates as reg2gene if gene names overlap. IMPORTANT! Gene coordinates 
#' are necessarilly a second anchor of the input reg2gene and benchData objects,
#' and column with gene names needs to be called "name".
#' @param mc.cores possible to be runned in parallel. Argument for mclapply f();
#' how many cores to use.
#' @param ignore.strand argument to be passed to
#' \code{\link[IRanges]{findOverlaps}}. When set to TRUE, the strand 
#' information is ignored in the overlap analysis.
#' @param preFilter (def:FALSE).  If TRUE, additional columns are added to the
#' input reg2Gene object (additionally to the Bench column that is reported by 
#' default) that
#' store info whether tested regions have any potential to be benchmarked.
#' Meaning, if all regulatory region-TSS pairs [anchor1 and anchor2 from
#' reg2Gene] do not overlap with any benchmark anchor1 or anchor2 location they
#' will be reported to be 0 (or no potential to be benchmarked at all), 
#' otherwise it is 1 (possible to be benchmarked).E.g. it selects reg2Gene 
#' regions only when both regulatory region and TSS  have overlapping regions 
#' somewhere in the benchmarking set; across all benchmark anchor pairs, but 
#' not necessarily overlapping regions of the same benchmark pair.
#' This info is
#' important to a priori remove high number of true negatives in 
#' regulatoryReg-TSS pairs, before running \code{\link{confusionMatrix}} since
#' TN are very abundant in the reg2gene dataset since benchmark dataset usually
#' covers much smaller regions of the genome (method limitations)
#' @param ... further arguments to methods, not implemented yet
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
#' 
#' Additionally, details for preFilter option:
#' All benchmark regions that can be confirmed by any combination of 
#' enh/gene pairs [anchor1 or anchor2 form reg2Gene object] is obtained. Then
#' selected unique union of anchor1 or anchor2 form reg2Gene object is used
#' as anchor1-anchor2 pairs that can be benchmarked. Reasoning, if present in
#' this set anchor1 or anchor2 regions form reg2Gene object necessarily need to
#' have other member of the pair overlapping somewhere in benchmark dataset. 
#' @author Inga Patarcic
#' 
#' @import GenomicRanges
#' @import InteractionSet
#' 
#' @examples # Creating testing and benchmarking dataset
#' require(GenomicRanges)
#' require(InteractionSet)
#'    
#'    reg2Gene <- GInteractions(GRReg1_toy,GRReg1_toy$reg)
#'    benchData <- GInteractions(GRReg2_toy,GRReg2_toy$reg)
#'    
#' benchmarkAssociations(reg2Gene,
#'             benchData,
#'             binary=FALSE)
#'             
#' benchmarkAssociations(reg2Gene,
#'             benchData,
#'             binary=TRUE)             
#' 
#' # add prefilter
#' 
#' benchmarkAssociations(reg2Gene,
#'             benchData,
#'             binary=TRUE,
#'             preFilter=TRUE) 
#'             
#'  # forceByName argument           
#'            
#'     reg2Gene$name <- reg2Gene$anchor1.name
#'     benchData$name <- benchData$anchor1.name
#'  
#'              benchmarkAssociations(reg2Gene,
#'                       benchData,
#'                       binary=TRUE,
#'                       forceByName = TRUE)
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
#' reg2GeneB <- benchmarkAssociations(reg2Gene,
#'                         benchDataList,
#'                         ignore.strand=TRUE,
#'                         binary=FALSE,
#'                         mc.cores = 1)               
#'     
#'  # forceByName = T can work for benchmark lists as well
#'  
#'             benchDataList <- list(benchData,benchData[1:5])
#'             names(benchDataList) <- c("BL1","BL2")
#' 
#' 
#'   benchmarkAssociations(benchData=benchDataList,
#'                         reg2Gene = reg2Gene,
#'                         forceByName = TRUE)
#'                             
#'                             
#'  # Checking what happends when anchor1&anchor2 both overlap only one region
#'  # in benchmark dataset?  OK, They are not benchmarked...
#'  
#'    benchmarkAssociations(reg2Gene[1],
#'                benchData[1])      
#'
#'  # WARNING! 
#'  # However, one need to be careful when benchmarking anchors that overlap 
#'  # within test set (eg enhancer overlaps gene region), 
#'  # because these regions will be benchmarked.
#'    
#'    benchmarkAssociations(reg2Gene[5],
#'                benchData)           
#'                                     
#' @export      
benchmarkAssociations <- function(reg2Gene,
                        benchData,
                        preFilter=FALSE,
                        binary=FALSE,
                        forceByName=FALSE,
                        mc.cores=1,
                        ignore.strand=FALSE,
                        ...) {
  
  
  
  if (class(benchData)=="GInteractions"){
    
    if (forceByName==T){
      
      # force benchData to have equal gene coordinates as reg2gene 
        benchData <- forceByname(benchData = benchData,reg2Gene=reg2Gene)
      
    }
    
    BenchmarkRes <- benchmarkAssociationssimple(reg2Gene=reg2Gene,
                                      benchData=benchData,
                                      binary=binary,
                                      ignore.strand=ignore.strand)
    
    
    if (preFilter==TRUE){
      
      # add info about filtering procedure
      
            filterRes <- filterPreBenchGIsimple(reg2Gene=reg2Gene,
                                          benchData=benchData,
                                          ignore.strand=ignore.strand)
      
            BenchmarkRes$Filter <- filterRes$Filter
    }
    
       return(BenchmarkRes)
    
  }
  
  
  
  if (class(benchData)=="list"){
    
    
    if (forceByName==T){
      
      # force benchData to have equal gene coordinates as reg2gene 
      benchData <- lapply(benchDataList, function(x){
                            return(forceByname(benchData=x,
                                        reg2Gene=reg2Gene))
                          })
      
    }
    
    
    BenchmarkRes <- parallel::mclapply(benchData, function(x) {
      
      benchmarkAssociationssimple(reg2Gene=reg2Gene,
                        benchData=x,
                        binary=binary,
                        ignore.strand=ignore.strand)},
      
      mc.cores = mc.cores)
    
    # pooling results together
    
    BenchResults <- DataFrame(lapply(BenchmarkRes,function(x) x$Bench))
    
    mcols(reg2Gene) <- c(mcols(reg2Gene),BenchResults)
    
    
    if (preFilter==TRUE){
      # add info about filtering if requested
        filterRes <- parallel::mclapply(benchData, function(x) {
        
        filterPreBenchGIsimple(reg2Gene=reg2Gene,
                               benchData=x,
                               ignore.strand=ignore.strand)},
        
        mc.cores = mc.cores)
      
      # pooling results together
      
      filterResults <- DataFrame(lapply(filterRes,function(x) x$Filter))
      
      # adjusting names for filtering
      names(filterResults) <- paste0("Filter_",names(filterResults))
      
      mcols(reg2Gene) <- c(mcols(reg2Gene),filterResults)
      
    }
    
    
    return(reg2Gene)
    
    
    
  }  
  
  
}



#' Benchmarks help function 
#' 
#' Benchmark reg2gene models using benchmark data but for only one GInteraction
#' object [not lists]
#' 
#' Read description for \code{benchmarkAssociations}
#' 
#' @author Inga Patarcic
#' @keywords internal
benchmarkAssociationssimple <- function(reg2Gene,
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
#' \code{benchmarkAssociations} and OPTIONAL  \code{filterPreBenchGI} 
#' metadata. However,prior this analysis, it is advised to reduce the number of
#' true negatives by including only reg2gene entries that could be potentially
#' benchmarked.To get that info run \code{filterPreBenchGI} and add filter 
#' column with metadata from this function to reg2GeneBench GInteractions 
#' object prior running \code{benchmarkAssociations})
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
#'    reg2GeneBench$PValue <- seq(0, 1,length.out = length(GRReg1_toy))
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
                            benchCol="Bench",
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

 
#' Benchmarks help function - forceByName
#' 
#' Forces benchmark input object to use gene coordinates of reg2Gene object if
#' gene names overlap. It is crucial that benchData & reg2Gene both contain
#' meta-data named "name" -> where gene names are stored, 
#' and that anchor2 in both objects corresponds to gene object. 
#' 
#' Read description for \code{benchmarkAssociations}
#' 
#' @author Inga Patarcic
#' @keywords internal
forceByname <- function(benchData,reg2Gene){
  
  # identifying gene names present in benchmark object and object used to 
  # benchmark
  
      MatchPairs <- match(benchData$name,reg2Gene$name)
      
      MatchPairs <- cbind(which(complete.cases(MatchPairs)),
                          MatchPairs[complete.cases(MatchPairs)])
      
  # extracting gene GRanges object
      
      tmpbenchData <- second(benchData)
  # matching    
      tmpbenchData[MatchPairs[,1]] <- second(reg2Gene)[MatchPairs[,2]]
   
      tmpelementsMeta <- elementMetadata(benchData)
  # back to GInteractions object    
      benchData <- GInteractions(first(benchData),
                                 tmpbenchData)
      
      elementMetadata(benchData) <- tmpelementsMeta

      return(benchData)
}
