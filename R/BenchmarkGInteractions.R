#####################################################
#-----------------Benchmark f()---------------------#
#####################################################

#' Benchmarks interactions reg2gene models using benchmark data
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
#' @param interactions a \code{\link[InteractionSet]{GInteractions}} object 
#' output from \code{\link{associateReg2Gene}}). 
#' Usually, 1st GRanges object, or anchor1
#' corresponds to the enahncer location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. 
#' @param benchInteractions a \code{\link[InteractionSet]{GInteractions}} object output
#' from \code{\link{associateReg2Gene}}) or a list of GInteractions object. 
#' Both regions are used in the benchmarking procedure. 
#' This object stores benchmarking informations eg interacting region
#' coordinates from techniques such as HiC,eQTL studies...
#' @param binary (def:FALSE) how many times reg2Gene interactions is observed in
#' the benchmark dataset(s). If TRUE, reports if overlap with benchmark dataset 
#' is observed at least once).
#' @param forceByName (def:FALSE) force benchmark data to have an equal gene 
#' coordinates as interactions if gene names overlap. IMPORTANT! Gene
#' coordinates are necessarilly a second anchor of the input interactions and 
#' benchInteractions objects,and column with gene names needs to be called "name".
#' @param mc.cores possible to be runned in parallel. Argument for mclapply f();
#' how many cores to use.
#' @param ignore.strand argument to be passed to
#' \code{\link[IRanges]{findOverlaps}}. When set to TRUE, the strand 
#' information is ignored in the overlap analysis.
#' @param preFilter (def:FALSE).  If TRUE, additional columns are added to the
#' input interactions object (additionally to the Bench column that is reported
#' by default) that
#' store info whether tested regions have any potential to be benchmarked.
#' Meaning, if all regulatory region-TSS pairs [anchor1 and anchor2 from
#' interactions] do not overlap with any benchmark anchor1 or anchor2 location 
#' they will be reported to be 0 (or no potential to be benchmarked at all), 
#' otherwise it is 1 (possible to be benchmarked).E.g. it selects interactions 
#' regions only when both regulatory region and TSS  have overlapping regions 
#' somewhere in the benchmarking set; across all benchmark anchor pairs, but 
#' not necessarily overlapping regions of the same benchmark pair.
#' This info is
#' important to a priori remove high number of true negatives in 
#' regulatoryReg-TSS pairs, before running \code{\link{confusionMatrix}} since
#' TN are very abundant in the interactions dataset since benchmark dataset 
#' usually covers much smaller regions of the genome (method limitations)
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
#' \code{\link[InteractionSet]{linkOverlaps}} between interactions 
#' and benchmark object is performed, and for each input pair it is reported 
#' whether this pair is benchmarked or not,and how many times (if binary=F).
#' Criss-cross overlap of interacting regions is performed; thus is anchor1 from
#' benchmark dataset is overlapping anchor2 from tested dataset, than anchor2
#' from benchmark dataset needs to overlap anchor1 from tested datased, or
#' vice-versa. 
#' 
#' Additionally, details for preFilter option:
#' All benchmark regions that can be confirmed by any combination of 
#' enh/gene pairs [anchor1 or anchor2 form interactions object] is obtained. 
#' Then
#' selected unique union of anchor1 or anchor2 form interactions object is used
#' as anchor1-anchor2 pairs that can be benchmarked. Reasoning, if present in
#' this set anchor1 or anchor2 regions form interactions object necessarily need
#' to have other member of the pair overlapping somewhere in benchmark dataset. 
#' @author Inga Patarcic
#' 
#' @import GenomicRanges
#' @import InteractionSet
#' 
#' @examples # Creating testing and benchmarking dataset
#' require(GenomicRanges)
#' require(InteractionSet)
#'    
#'    interactions <- GInteractions(GRReg1_toy,GRReg1_toy$reg)
#'    benchInteractions <- GInteractions(GRReg2_toy,GRReg2_toy$reg)
#'    
#' benchmarkInteractions(interactions,
#'             benchInteractions,
#'             binary=FALSE)
#'             
#' benchmarkInteractions(interactions,
#'             benchInteractions,
#'             binary=TRUE)             
#' 
#' # add prefilter
#' 
#' benchmarkInteractions(interactions,
#'             benchInteractions,
#'             binary=TRUE,
#'             preFilter=TRUE) 
#'             
#'  # forceByName argument           
#'            
#'     interactions$name <- interactions$anchor1.name
#'     benchInteractions$name <- benchInteractions$anchor1.name
#'  
#'              benchmarkInteractions(interactions,
#'                       benchInteractions,
#'                       binary=TRUE,
#'                       forceByName = TRUE)
#'                 
#' ##################   
#' # example for list:
#' 
#' # NOTE: anchor1.Bench1Exp & anchor1.Bench2Exp are expected/precalculated
#' # values for this benchmark example
#'  
#'   benchInteractionsList <- list(benchInteractions,interactions)
#'   names(benchInteractionsList) <- c("benchInteractions1",
#'   "benchInteractions2")
#' 
#' 
#' interactionsB <- benchmarkInteractions(interactions,
#'                         benchInteractionsList,
#'                         ignore.strand=TRUE,
#'                         binary=FALSE,
#'                         mc.cores = 1)               
#'     
#'  # forceByName = T can work for benchmark lists as well
#'  
#'             benchInteractionsList <- list(benchInteractions,
#'             benchInteractions[1:5])
#'             names(benchInteractionsList) <- c("BL1","BL2")
#' 
#' 
#'   benchmarkInteractions(benchInteractions=benchInteractionsList,
#'                         interactions = interactions,
#'                         forceByName = TRUE)
#'                             
#'                             
#'  # Checking what happends when anchor1&anchor2 both overlap only one region
#'  # in benchmark dataset?  OK, They are not benchmarked...
#'  
#'    benchmarkInteractions(interactions[1],
#'                benchInteractions[1])      
#'
#'  # WARNING! 
#'  # However, one need to be careful when benchmarking anchors that overlap 
#'  # within test set (eg enhancer overlaps gene region), 
#'  # because these regions will be benchmarked.
#'    
#'    benchmarkInteractions(interactions[5],
#'                benchInteractions)           
#'                                     
#' @export      
benchmarkInteractions <- function(interactions,
                        benchInteractions,
                        preFilter=FALSE,
                        binary=FALSE,
                        forceByName=FALSE,
                        mc.cores=1,
                        ignore.strand=FALSE,
                        ...) {
  
  
  
  if (class(benchInteractions)=="GInteractions"){
    
    if (forceByName==T){
      
      # force benchInteractions to have equal gene coordinates as interactions 
        benchInteractions <- forceByname(benchInteractions = benchInteractions,
                                 interactions=interactions)
      
    }
    
    BenchmarkRes <- benchmarkInteractionssimple(interactions=interactions,
                                      benchInteractions=benchInteractions,
                                      binary=binary,
                                      ignore.strand=ignore.strand)
    
    
    if (preFilter==TRUE){
      
      # add info about filtering procedure
      
            filterRes <- filterPreBenchGIsimple(interactions=interactions,
                                          benchInteractions=benchInteractions,
                                          ignore.strand=ignore.strand)
      
            BenchmarkRes$Filter <- filterRes$Filter
    }
    
       return(BenchmarkRes)
    
  }
  
  
  
  if (class(benchInteractions)=="list"){
    
    
    if (forceByName==T){
      
      # force benchInteractions to have equal gene coordinates as interactions 
      benchInteractions <- lapply(benchInteractionsList, function(x){
                            return(forceByname(benchInteractions=x,
                                        interactions=interactions))
                          })
      
    }
    
    
    BenchmarkRes <- parallel::mclapply(benchInteractions, function(x) {
      
      benchmarkInteractionssimple(interactions=interactions,
                        benchInteractions=x,
                        binary=binary,
                        ignore.strand=ignore.strand)},
      
      mc.cores = mc.cores)
    
    # pooling results together
    
    BenchResults <- DataFrame(lapply(BenchmarkRes,function(x) x$Bench))
    
    mcols(interactions) <- c(mcols(interactions),BenchResults)
    
    
    if (preFilter==TRUE){
      # add info about filtering if requested
        filterRes <- parallel::mclapply(benchInteractions, function(x) {
        
        filterPreBenchGIsimple(interactions=interactions,
                               benchInteractions=x,
                               ignore.strand=ignore.strand)},
        
        mc.cores = mc.cores)
      
      # pooling results together
      
      filterResults <- DataFrame(lapply(filterRes,function(x) x$Filter))
      
      # adjusting names for filtering
      names(filterResults) <- paste0("Filter_",names(filterResults))
      
      mcols(interactions) <- c(mcols(interactions),filterResults)
      
    }
    
    
    return(interactions)
    
    
    
  }  
  
  
}



#' Benchmarks help function 
#' 
#' Benchmark interactions models using benchmark data but for only one 
#' GInteraction object [not lists]
#' 
#' Read description for \code{benchmarkInteractions}
#' 
#' @author Inga Patarcic
#' @keywords internal
benchmarkInteractionssimple <- function(interactions,
                              benchInteractions,
                              ignore.strand,
                              binary) {
  
  # perform Overlaps 
  OverlapGI <- linkOverlaps(interactions,
                            first(benchInteractions),
                            second(benchInteractions),
                            ignore.strand=ignore.strand)
  
  
  # getting rows that are confirmed by benchmark
  
  benchRow <- OverlapGI$query[OverlapGI$subject1==OverlapGI$subject2]
  
  # reporting benchmark results
  
  interactions$Bench <- rep(0,length(interactions))
  
  if (binary==TRUE){ interactions$Bench[unique(benchRow)] <- 1}
  
  if (binary==FALSE){ 
    # add counts   
    interactions$Bench[as.integer(names(table(benchRow)))] <- table(benchRow)
    
  }
  
  
  return(interactions)
  
}




#' Filtering help function 
#' 
#' Filter interactions models using benchmark data but for only one GInteraction
#' object [not lists]
#' 
#' Read description for \code{filterPreBenchGI}
#' 
#' @author Inga Patarcic
#' @keywords internal 
filterPreBenchGIsimple <- function(interactions,
                          benchInteractions,
                          ignore.strand=FALSE){


         BenchCovered <- linkOverlaps(benchInteractions,
                                     first(interactions),
                                     second(interactions),
                                     ignore.strand=ignore.strand)
                 
         
         filterRow <- unique(c(BenchCovered$subject1,
                             BenchCovered$subject2))
         
         interactions$Filter <- rep(0,length(interactions))
         
         interactions$Filter[filterRow] <- 1
         
  return(interactions)
         
         
         }




#####################################################
#-----------------ConfusionMatrix f()---------------------#
#####################################################


#' ConfusionTable statistics for the benchmarked interactions object
#' 
#' A function takes as an benchmarked interactions object and outputs
#' confusion table statistics.
#' 
#' 
#' @param interactionsBench GInteractions object with added benchmark 
#' \code{benchmarkInteractions} and OPTIONAL  \code{filterPreBenchGI} 
#' metadata. However,prior this analysis, it is advised to reduce the number of
#' true negatives by including only interactions entries that could be 
#' potentially
#' benchmarked.To get that info run \code{filterPreBenchGI} and add filter 
#' column with metadata from this function to interactionsBench GInteractions 
#' object prior running \code{benchmarkInteractions})
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
#' output is list with following elements:TP,FP,TN,FN,Specificity, Accuracy,PPV,
#' NPV,F1. Otherwise only PPV is reported.
#' 
#' 
#' @details Reports statistics based on the confusion matrix. 
#' There are four different categories in the confusion matrix: 
#' TP = (number of) true positive: interactions entry that was reported to be 
#' associated (reported gene-enhancer statistics lower than a predefined
#' threshold) and was benchmarked. 
#' FP = (number of) false positive: interactions entry that was reported to be
#' associated (reported gene-enhancer statistics lower than a predefined 
#' threshold) BUT was not overlapped with benchmark dataset
#' FN = (number of) false negative: interactions entry that was NOT reported to
#' be associated (reported gene-enhancer statistics NOT lower than a predefined 
#' threshold) BUT is benchmarked
#' TN = (number of) true negative: interactions entry 
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
#'    interactionsBench <- GInteractions(GRReg1_toy,GRReg1_toy$reg)
#'    interactionsBench$PValue <- seq(0, 1,length.out = length(GRReg1_toy))
#' 
#'
#' confusionMatrix(interactionsBench,
#'                 thresholdID = "PValue",
#'                 thresholdValue = 0.05,
#'                 benchCol = "anchor1.Bench1Exp",
#'                 prefilterCol = "anchor1.Filter1Exp",
#'                 statistics = "ConfusionMatrix")
#'                 
#' confusionMatrix(interactionsBench,
#'                 thresholdID = "PValue", 
#'                 thresholdValue = 0.05,
#'                 benchCol = "anchor1.Bench1Exp",
#'                 prefilterCol = "anchor1.Filter1Exp",
#'                 statistics = "PPV")
#'                
#'                
#'      
#' @export
confusionMatrix <- function(interactionsBench,
                            benchCol="Bench",
                            prefilterCol=NULL,
                            thresholdID=NULL,
                            thresholdValue=0.05,
                            statistics="PPV") {
  
  
           if (is.null(thresholdID)){
              
              stop("Cannot run thresholding, no column defined")
              
              }
  
           if (!any(stringr::str_detect(colnames(mcols(interactionsBench)),
                                         thresholdID))){
              
              stop("thresholdID column not identified")
              
           }
  
  
          if (!any(stringr::str_detect(colnames(mcols(interactionsBench)),
                                       benchCol))){
            
              stop("benchCol column not detected")
            
          }  
  
  
  
  # filter results of benchmarking procedure by prefiltering column
  
          if (!is.null(prefilterCol)){ 
            
            rows.filt <- as.logical(mcols(interactionsBench)[,prefilterCol])
            
            interactionsBench <- interactionsBench[rows.filt]
          
            }
  
  
  ##################################
  # thresholding interactionsBench object
  
  
  predictedEntries <- mcols(interactionsBench)[,thresholdID] <= thresholdValue
  
  benchmarkedEntries <- as.logical(mcols(interactionsBench)[,benchCol])
  
  
  
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
#' Forces benchmark input object to use gene coordinates of interactions object 
#' if
#' gene names overlap. It is crucial that benchInteractions 
#' & interactions both contain
#' meta-data named "name" -> where gene names are stored, 
#' and that anchor2 in both objects corresponds to gene object. 
#' 
#' Read description for \code{benchmarkInteractions}
#' 
#' @author Inga Patarcic
#' @keywords internal
forceByname <- function(benchInteractions,interactions){
  
  # identifying gene names present in benchmark object and object used to 
  # benchmark
  
      MatchPairs <- match(benchInteractions$name,interactions$name)
      
      MatchPairs <- cbind(which(complete.cases(MatchPairs)),
                          MatchPairs[complete.cases(MatchPairs)])
      
  # extracting gene GRanges object
      
      tmpbenchInteractions <- second(benchInteractions)
  # matching    
      tmpbenchInteractions[MatchPairs[,1]] <- second(interactions)[MatchPairs[,2]]
   
      tmpelementsMeta <- elementMetadata(benchInteractions)
  # back to GInteractions object    
      benchInteractions <- GInteractions(first(benchInteractions),
                                 tmpbenchInteractions)
      
      elementMetadata(benchInteractions) <- tmpelementsMeta

      return(benchInteractions)
}
