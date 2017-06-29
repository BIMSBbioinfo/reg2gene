#' Function that swaps gene and regulatory region locations
#' 
#' Function that swaps 2 GRanges objects(one stored as a meta-data within the 
#' other).Thus,reg meta-data stores info about gene locations, wheres main 
#' ranges store ranges previously stored as a reg meta-data 
#' 
#' @param reg GRanges object with gene locations and meta-data which stores
#' info about regulatory region locations (as a GRanges object) 
#' @param reg.col  character, a meta-data column name which contains regulatory
#' region locations 
#' @author Inga Patarcic
#' @keywords internal
switchReg = function(reg, regCol){
  vals = values(reg)
  cols = vals[,colnames(vals) != regCol]
  reg2 = values(reg)[[regCol]]
  
  
  values(reg2) = cbind(S4Vectors::DataFrame(granges(reg)), cols)
  colnames(values(reg2))[1] = regCol
  reg2
  
}


#' Overlap of interacting regions   
#' 
#' The function that compares two sets of interacting regions in the genome 
#' and reportes interacting regions from the first datasets that is present/
#' confirmed by the second dataset of interacting regions. 
#'       
#' @param reg1 - GRanges object with at least one meta-data column which stores
#' the second GRanges object. Usually,one GRanges object corresponds to the gene 
#' locations, preferably TSS location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. Meta-data column name should
#' be "reg".
#' @param reg2 - GRanges object with at least one meta-data column which stores
#' the second GRanges object.  Meta-data column name should be "reg".
#'   
#' @return GRanges object with interacting regions from the first datasets 
#' that was present/confirmed by the second dataset of interacting regions. 
#' Meta-data from both input GRanges objects (storing info about interacting 
#' regions in the genome) are included in the final GRanges object, and 
#' meta-data column names are altered such that they indicate their origin (
#' eg from Reg1 object or Reg2 object).
#'  
#' @details Function that performs a simple findOverlapsoverlap twice - 
#' between 2 sets of coordinates stored.  It returns rows in which both 
#' combinations of overlaps are TRUE. 
#' @author Inga Patarcic
#' @keywords internal
overlapRegions <- function(reg1,reg2){
  
  # reg1-reg1 overlap            
  reg1.reg1.Overlap <- data.frame(findOverlaps(reg1,reg2))
  
  # adjusting names
  colnames(mcols(reg1)) <- paste0("reg1_",colnames(mcols(reg1)))
  colnames(mcols(reg2)) <- paste0("reg2_",colnames(mcols(reg2)))
  
  # adding Overlapped regions and mcols on top of those we want to test       
  reg1reg2 <- reg1[reg1.reg1.Overlap$queryHits,]
  reg1reg2$reg2coord1 <- (reg2[reg1.reg1.Overlap$subjectHits,])
  mcols(reg1reg2) <- cbind(mcols(reg1reg2),
                           mcols(reg2[reg1.reg1.Overlap$subjectHits,]))
  
  ########## 
  # reg1-reg1 overlap            
  reg2.reg2.Overlap <- data.frame(findOverlaps(reg1reg2$reg1_reg,
                                               reg1reg2$reg2_reg))
  
  # if overlap is duplicated then reg1 overlaps reg1, and reg2 overlaps 
  # reg2 -> confirmed
  AlloverlapRegions <- which(reg2.reg2.Overlap$queryHits == 
                               reg2.reg2.Overlap$subjectHits)
  InputRows.Overlapped <- reg2.reg2.Overlap[AlloverlapRegions,"queryHits"]
  reg1reg2 <- reg1reg2[InputRows.Overlapped,]
  reg1reg2$Coor1Coord2PAIR <- reg1.reg1.Overlap$queryHits[InputRows.Overlapped]
  
  
  return(reg1reg2)
  
}



#' Criss-cross overlap of interacting regions 
#' 
#' A function that compares four input regions,eg two interacting regions
#' from one dataset with two interacting regions from the other dataset
#' in the criss-cross manner. THUS, the input orientation of the 
#' interacting regions is NOT strict: the first GRanges object of the region1
#' from one dataset will be compared to both; the first GRanges object and the
#' second GRanges object of region2(and its pair is compared in the opposite 
#' direction). 
#'       
#' @param reg1 - GRanges object with at least one meta-data column which stores
#' the second GRanges object. Usually,one GRanges object corresponds to the gene 
#' locations, preferably TSS location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. Meta-data column name should
#' be "reg".
#' @param reg2 - GRanges object with at least one meta-data column which stores
#' the second GRanges object.  Meta-data column name should be "reg".
#' @param reg.col  character, a meta-data column name which contains regulatory
#' region locations for the reg2 object.
#'
#' @return GRanges object with interacting regions from the first datasets 
#' that was present/confirmed by the second dataset of interacting regions. 
#' Meta-data from both input GRanges objects (storing info about interacting 
#' regions in the genome) are included in the final GRanges object, and 
#' meta-data column names are altered such that they indicate their origin (
#' eg from Reg1 object or Reg2 object).
#'  
#' @details Function performs criss-cross overlap between 2 GRanges objects with
#' 2 pairs of GRanges each. It performs reg2gene::OverlapRegions 
#' overlap between  (reg1:Coord1-reg2:Coord1)&(reg1:Coord2-reg2:Coord2) and 
#' vice-versa (reg1:Coord1-reg2:Coord2)&(reg1:Coord2-reg2:Coord1).
#' Switch is achieved by swaping coordinated of reg2 (eg positions of gene and
#' regulatory region locations stored as GRanges object within GRanges object 
#' are switched).
#' @author Inga Patarcic
#' @keywords internal
complexOverlaps <- function(reg1,reg2, reg.col){
  
  # criss-cross overlap when orientation of overlaps is unknown
  reg11reg22 <- overlapRegions(reg1=reg1,reg2=reg2)
  # switch names in reg1 dataframe
  
  reg12reg21 <- overlapRegions(reg1=reg1,reg2=switchReg(reg2,regCol=reg.col))
  
  
  return(c(reg11reg22,reg12reg21))
  
  
}


#' Benchmarks reg2gene models using benchmark data
#' 
#  The function that takes as an input results of \code{associatereg2Gene} or
#' any other modelling procedure implemente in \code{reg2gene} package,
#' and predefined benchmark dataset as GRanges object.
#' This function adds a metadata column with info about benchmarking success -
#' whether tested regions are benchmarked or not.  
#' By default it reportes how many times interactions is observed in the 
#' benchmark dataset. If binary is set to TRUE, then logical vector of TRUE 
#' (overlapping benchmark dataset at least once) and FALSE (not overlapping
#' benchmark dataset at all) is reported.
#' 
#' 
#' @param reg2Gene a GRanges object output (from \code{associatereg2Gene}).
#' GRanges object should have at least one meta-data column which stores
#' the second GRanges object. Usually,one GRanges object corresponds to the gene 
#' locations, preferably TSS location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. Meta-data column name should
#' be "reg".
#' @param benchReg a GRanges object with at least one meta-data column which 
#' stores the second GRanges object named "reg".
#' @param regCol character (default "reg"), a column name of meta-data object.
#' Indicates the column where the location of the 2nd pair of GRanges is stored,
#' (usually regulatory region locations)
#' @param binary (def:FALSE) how many times reg2Gene interactions is observed in
#' the benchmark dataset. If TRUE, reports if overlap with benchmark dataset is
#' observed at least once).
#' 
#' @return GRanges object which is equall to the reg2Gene input object except 
#' that a meta-data column reporting how many times interaction was overlapped
#' or whether reported interaction was overlapped is added and 
#' named "BenchmarkO"
#' 
#' 
#' @details 
#GRanges objects - an output of \code{associateReg2Gene} and 
#' benchmark dataset is inputed are overlaped. ComplexOverlap between reg2Gene 
#' and benchmark object is performed, and for each input pair it is reported 
#' whether this pair is benchmarked or not,and how many times (if binary=F).
#' ComplexOverlap means that criss-cross overlap of interacting regions is 
#' performed; eg four input regions are compared: two interacting regions
#' from one dataset with two interacting regions from the other dataset
#' in the criss-cross manner. Thus, the input orientation of the 
#' interacting regions is not strict: the first pair of the region1
#' from one dataset will be compared to both; the first pair and the
#' second pair of region2(and its pair is compared in the opposite 
#' direction).  
#' @author Inga Patarcic
#' 
#' @examples benchmark(GRReg1_toy,GRReg2_toy)
#' benchmark(GRReg1_toy,GRReg2_toy,binary=TRUE)
#' @export
benchmark <- function(reg2Gene,
                      benchReg,
                      regCol="reg",
                      binary=FALSE){
    
  ############################### 
  # benchmarking           

        # Complexbenchmarking   
      reg2GeneBenchOverlap <- complexOverlaps(reg2Gene,
                                              benchReg,
                                              regCol)
 
          # counts or binary reported   
        OverlapVector <- rep(0,length(reg2Gene))
        Duplicatedbenchmarks <- table(reg2GeneBenchOverlap$Coor1Coord2PAIR)
  OverlapVector[as.integer(names(Duplicatedbenchmarks))] <- Duplicatedbenchmarks
          
        if (binary==TRUE){ OverlapVector <- as.logical(OverlapVector)}
           
            reg2Gene$BenchmarkO <- OverlapVector
            
            
  return(reg2Gene)
        
}
  
  
#' ConfusionTable statistics for the benchmarked reg2Gene object
#' 
#' A function takes as an benchmarked reg2gene object and outputs
#' confusion table statistics.
#' 
#' 
#' @param reg2GeneBench Granges object which contains a metacolumn(type:logical) 
#' with the results of benchmarking procedure. It can be a benchmarked 
#' reg2Gene object resulted from \code{benchmark}.
#' !WARNING! Prior this analysis, it is advised to reduce the number of 
#' true negatives by including only reg2gene entries that could be potentially
#' benchmarked. THUS run \code{Filter_PreBench} prior running \code{benchmark})
#' 
#' @param benchCol character (default "benchmarkO"), a column name of metadata
#' that indicates the column where the result of benchmarking procedure is 
#' stored (as vector of TRUE and FALSE entries).
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
#' @return Integer or list. If "ConfusionMatrix" is requested then the output is 
#' list with following values:TP,FP,TN,FN,Specificity, Accuracy, PPV, NPV,F1.
#' Otherwise only PPV is reported.
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
#' 
#' The formulas used in this function are: \deqn{Sensitivity = TP/(TP+FN)} 
#' \deqn{Specificity = TN/(TN+FP)} 
#' \deqn{Accuracy = (TP+TN)/(TP+FN+FP+TN)}
#' \deqn{PPV = TP/(TP+FP)} 
#' \deqn{NPV = TN/(TN+FN)} 
#' \deqn{F1 = (2*TP)/((2*TP)+FP+FN)}
#' 
#' @examples confusionMatrix(BenchMarkedReg2Gene_toy,thresholdID = "PValue",
#' thresholdValue = 0.05,benchCol = "BenchmarkO",statistics = "ConfusionMatrix")
#' confusionMatrix(BenchMarkedReg2Gene_toy,thresholdID = "PValue",
#'   thresholdValue = 0.5,benchCol = "BenchmarkO",statistics = "PPV")
#' 
#'      
#' @export
confusionMatrix <- function(reg2GeneBench,
                            benchCol="BenchmarkO",
                            thresholdID=NULL,
                            thresholdValue=0.05,
                            statistics="PPV") {
  
   
  if (is.null(thresholdID)){stop("Cannot run thresholding, no column defined")}
  if (!any(stringr::str_detect(colnames(mcols(reg2GeneBench)),thresholdID))){
    stop("thresholdID column not identified")}  
  
  if (!any(stringr::str_detect(colnames(mcols(reg2GeneBench)),benchCol))){
    stop("benchCol column not detected")}  
  
 
  
    
  ##################################
  # thresholding reg2GeneBench object
  predictedEntries <- mcols(reg2GeneBench)[,thresholdID] <= thresholdValue
  benchmarkedEntries <- mcols(reg2GeneBench)[,"BenchmarkO"]
            
  Confusion.matrix <- table(predictedEntries,benchmarkedEntries)
  
            TP <- Confusion.matrix["TRUE","TRUE"]
            TN <- Confusion.matrix["FALSE","FALSE"]
            FN <- Confusion.matrix["FALSE","TRUE"]
            FP <- Confusion.matrix["TRUE","FALSE"]
            
  
  
  if (statistics=="PPV"){ CM.Statistics <- return(TP/(TP+FP))}
  
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

