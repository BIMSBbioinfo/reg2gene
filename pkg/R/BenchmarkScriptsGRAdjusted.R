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
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
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
overlapRegions <- function(reg1,reg2, strandness=FALSE){
  
  # reg1-reg1 overlap            
  reg1.reg1.Overlap <- data.frame(findOverlaps(reg1,reg2, 
                                               ignore.strand=strandness))

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
                                               reg1reg2$reg2_reg,
                                               ignore.strand=strandness))
  
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
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
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
complexOverlaps <- function(reg1,reg2, reg.col,strandness=FALSE){
  
  # criss-cross overlap when orientation of overlaps is unknown
  reg11reg22 <- overlapRegions(reg1=reg1,reg2=reg2, strandness=strandness)
  # switch names in reg1 dataframe
  
  reg12reg21 <- overlapRegions(reg1=reg1,reg2=switchReg(reg2,regCol=reg.col),
                               strandness=strandness)
  
  
  return(c(reg11reg22,reg12reg21))
  
  
}


#' Benchmarks reg2gene models using benchmark data
#' 
#'  The function that takes as an input results of \code{associatereg2Gene} or
#' any other modelling procedure implemente in \code{reg2gene} package,
#' and predefined benchmark dataset as GRanges object.
#' This function adds a metadata column with info about benchmarking success -
#' whether tested regions are benchmarked or not.  
#' By default it reportes how many times interactions is observed in the 
#' benchmark dataset. If binary is set to TRUE, then logical vector of TRUE 
#' (overlapping benchmark dataset at least once) and FALSE (not overlapping
#' benchmark dataset at all) is reported. Optionally, it can report a vector
#' or a dataframe of the result of benchmarking procedure.
#' 
#' 
#' @param reg2Gene a GRanges object output (from \code{associatereg2Gene}).
#' GRanges object should have at least one meta-data column which stores
#' the second GRanges object. Usually,one GRanges object corresponds to the gene 
#' locations, preferably TSS location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. Meta-data column name should
#' be "reg".
#' @param benchData a GRanges or GrangesList object with at least one meta-data 
#' column which stores the second GRanges object named "reg". Both regions
#' are used in benchmarking procedure. This object stores 
#' benchmarking informations eg interacting region coordinates from techniques 
#' such as HiC,eQTL studies...
#' @param regCol character (default "reg"), a column name of meta-data object.
#' Indicates the column where the 2nd pair of GRanges is stored,
#' (usually regulatory region locations)
#' @param binary (def:FALSE) how many times reg2Gene interactions is observed in
#' the benchmark dataset(s). If TRUE, reports if overlap with benchmark dataset 
#' is observed at least once).
#' @param reportGR (def:TRUE) Expected output format. If TRUE, GRanges object +
#' column with benchmarking results returned. If FALSE only vector of
#' benchmarking results reported.
#' @param nCores possible to be runned in parallel. Argument for mclapply f();
#' how many cores to use.
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
#' 
#' 
#' @return GRanges object or vector/dataframe of benchmarked results. When 
#' reportGR==T, GRanges object corresponding to the reg2Gene input object with 
#' added benchmarked meta-data column(s) is reported. Each column corresponds
#' to one benchmark dataset analyzed. Values can be either 0/1 (not/benchmarked)
#' or 0-n (how many times each gene-enhancer pair is benchmared). When 
#' reportGR==F a vector (one benchmark dataset was used as an input) or a 
#' dataframe (>1 benchmark dataset was used as an input) is reported.
#' 
#' 
#' @details GRanges objects - an output of \code{associateReg2Gene} and 
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
#' benchmark(filterPreBench(GRReg2_toy,GRReg1_toy),GRReg1_toy)
#' @export
benchmark <- function(reg2Gene,
                      benchData,
                      regCol="reg",
                      binary=TRUE,
                      nCores=1,
                      reportGR=FALSE,
                      strandness=FALSE,
                      ...){
  
  
  # checking format of benchmark  and tested input (correct columns,GrangesList)
  tmp <- checkGR(GRangesO=benchData,reg.definition=regCol)
  if (!any(colnames(mcols(reg2Gene)) %in% regCol)){stop("No regCol detected")}
  
  
  
  # run associateReg2Gene for gene chunks and save them separately
  
  
  if (class(benchData)=="GRangesList"){
    
    
    BenchmarkingResults <- benchmarkList(reg2Gene=reg2Gene,
                                         benchList=benchData,
                                         regCol=regCol,
                                         binary=binary,
                                         strandness=strandness)
    
    if (reportGR==TRUE) {mcols(reg2Gene) <- c(mcols(reg2Gene),
                                              BenchmarkingResults)
    return(reg2Gene)}
    
    if (reportGR==F) {return(BenchmarkingResults)}  
    
  }  
  
  if (class(benchData)=="GRanges"){
    
    
    BenchmarkingResults <- benchmarkSimple(reg2Gene=reg2Gene,
                                           benchReg=benchData,
                                           regCol=regCol,
                                           binary=binary,
                                           strandness=strandness)
    
    
    
    if (reportGR==TRUE) {reg2Gene$Bench <- BenchmarkingResults
    return(reg2Gene)}
    if (reportGR==F) {return(BenchmarkingResults)}  
  }
} 


#' Function that performs benchmarking procedure when benchData is GRanges obj
#' 
#' 
#' @param reg2Gene a GRanges object output (from \code{associatereg2Gene}).
#' GRanges object should have at least one meta-data column which stores
#' the second GRanges object.
#' @param benchData a GRanges object with at least one meta-data 
#' column which stores the second GRanges object named "reg". Both regions
#' are used in benchmarking procedure.
#' @param regCol character (default "reg"), a column name of meta-data object.
#' Indicates the column where the 2nd pair of GRanges is stored,
#' @param binary (def:FALSE) how many times reg2Gene interactions is observed in
#' the benchmark dataset(s). If TRUE, reports if overlap with benchmark dataset 
#' is observed at least once). 
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
#' @author Inga Patarcic
#' @keywords internal
benchmarkSimple <- function(reg2Gene,
                            benchReg,
                            regCol="reg",
                            binary=FALSE,
                            strandness=FALSE){
  
  ############################### 
  # benchmarking           
  
  
  # Complexbenchmarking   
  reg2GeneBenchOverlap <- complexOverlaps(reg2Gene,
                                          benchReg,
                                          regCol,
                                          strandness)
  
  # counts or binary reported   
  OverlapVector <- rep(0,length(reg2Gene))
  Duplicatedbenchmarks <- table(reg2GeneBenchOverlap$Coor1Coord2PAIR)
  OverlapVector[as.integer(names(Duplicatedbenchmarks))] <- Duplicatedbenchmarks
  
  # to report vector of Benchmarking results   
  if (binary==T) {return(as.logical(OverlapVector))}
  
  if (binary==F) {return(OverlapVector)}
  
  
}




#' Function that performs benchmarking procedure when benchData is GRangesList
#' 
#' 
#' @param reg2Gene a GRanges object output (from \code{associatereg2Gene}).
#' GRanges object should have at least one meta-data column which stores
#' the second GRanges object.
#' @param benchList a GrangesList object with at least one meta-data 
#' column which stores the second GRanges object named "reg". Both regions
#' are used in benchmarking procedure.
#' @param regCol character (default "reg"), a column name of meta-data object.
#' Indicates the column where the 2nd pair of GRanges is stored,
#' @param binary (def:FALSE) how many times reg2Gene interactions is observed in
#' the benchmark dataset(s). If TRUE, reports if overlap with benchmark dataset 
#' is observed at least once).
#' @param nCores possible to be runned in parallel
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
#' @author Inga Patarcic
#' @keywords internal
benchmarkList <- function(reg2Gene,
                          benchList,
                          regCol="reg",
                          binary=FALSE,
                          nCores=1,
                          strandness=FALSE,
                          ...){
  
  # running benchmarking per benchmark dataset
  bench_list <- mclapply(benchList, function(x) {
    
    benchmarkSimple(reg2Gene=reg2Gene,
                    benchReg=x,
                    regCol = regCol,
                    binary = binary,
                    strandness=strandness)},
    mc.cores = nCores)
  
  bench_list <- do.call("cbind.data.frame",bench_list)
  
  
  
  return(bench_list)
  
  
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
#' @param benchCol character (default NULL), results of benchmarking 
#' procedure; eg a column name of the metadata which indicates the column where 
#' the result of benchmarking procedure is stored
#' 
#' @param prefilterCol character (default NULL), results of prefiltering 
#' procedure; eg a column name of the metadata which indicates the column where 
#' the result of prefiltering procedure is stored (as vector of TRUE and FALSE 
#' entries).
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
#' @examples confusionMatrix(BenchMarkedReg2Gene_toy,thresholdID = "PValue",
#' thresholdValue = 0.05,benchCol = "BenchmarkO",statistics = "ConfusionMatrix") 
#' confusionMatrix(BenchMarkedReg2Gene_toy,thresholdID = "PValue",
#'   thresholdValue = 0.5,benchCol = "BenchmarkO",statistics = "PPV")
#'   # confusion matrix statistics with prefiltering step
#'   confusionMatrix(filterPreBench(BenchMarkedReg2Gene_toy,GRReg2_toy),
#'   thresholdID = "PValue",
#'   thresholdValue = 0.05,
#'   benchCol = "BenchmarkO",
#'   statistics = "ConfusionMatrix")
#' BenchMarkedReg2Gene_toy$prefilter <- c(rep(1,4),0,0)
#' confusionMatrix(reg2GeneBench=BenchMarkedReg2Gene_toy,
#'                thresholdID = "PValue", 
#'                thresholdValue = 0.05,
#'                benchCol = "BenchmarkO",
#'                prefilterCol = "prefilter",
#'                statistics = "ConfusionMatrix")
#'      
#' @export
confusionMatrix <- function(reg2GeneBench,
                            benchCol=NULL,
                            prefilterCol=NULL,
                            thresholdID=NULL,
                            thresholdValue=0.05,
                            statistics="PPV") {
  
   
  if (is.null(thresholdID)){stop("Cannot run thresholding, no column defined")}
  if (!any(stringr::str_detect(colnames(mcols(reg2GeneBench)),thresholdID))){
    stop("thresholdID column not identified")}  
  
  if (!any(stringr::str_detect(colnames(mcols(reg2GeneBench)),benchCol))){
    stop("benchCol column not detected")}  

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




#' Eliminates all regulatory regions and TSSs that do not overlap with benchmark
#' dataset of choice (eliminates high number of true negatives)
#'
#' Function the eliminates all enhancers and TSSs that do not overlap with 
#' benchmark dataset, eg it selects reg2gene regions only when regulatory region 
#' and TSS overlap with at least one benchmarking region.  
#' This is useful to improve the confusion matrix statistics reported by 
#' \code{confusionMatrix} by eliminating the high number of true negatives. TN 
#' are very abundant in the reg2gene dataset since benchmark dataset usually 
#' covers much smaller regions of the genome (method limitations)
#'
#' @param reg2Gene a GRanges object output (from \code{associatereg2Gene}).
#' GRanges object should have at least one meta-data column which stores
#' the second GRanges object. Usually,one GRanges object corresponds to the gene 
#' locations, preferably TSS location, whereas the other GRanges object 
#' corresponds to the regulatory region locations. Meta-data column name should
#' be "reg".
#' @param benchData a GRanges or GrangesList object with at least one meta-data 
#' column which stores the second GRanges object named "reg". Both regions
#' are used in benchmarking procedure. This object stores 
#' benchmarking informations eg interacting region coordinates from techniques 
#' such as HiC,eQTL studies, GWAS
#' @param regCol character (default "reg"), a column name of meta-data object.
#' Indicates the column where the 2nd pair of GRanges is stored,
#' (usually regulatory region locations)
#' @param reportGR (def:TRUE) Expected output format. If TRUE, GRanges object +
#' column with benchmarking results returned. If FALSE only vector of
#' benchmarking results reported.
#' @param nCores possible to be runned in parallel. Argument for mclapply f();
#' how many cores to use.
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
#' 
#' @details Each reg2gene pair is overlapped with the benchmark dataset 
#' (both regions of the benchmark dataset). If both: 1) the regulatory region 
#' and 2) TSS from tested reg2gene pair overlap with at least one benchmarking 
#' region, then this pair is kept for the benchmarking analyses.
#' @examples filterPreBench(GRReg2_toy,GRReg1_toy)
#' GR.list <- GRangesList(GRReg1_toy,GRReg2_toy)
#' names(GR.list) <- c("GRReg1_toy","GRReg2_toy")
#' filterPreBench(GRReg2_toy,GR.list)
#' 
#' @export
filterPreBench <- function(reg2Gene,
                           benchData,
                           regCol="reg",
                           nCores=1,
                           reportGR=FALSE,
                           strandness=FALSE,
                           ...){
  
  
  # checking format of benchmark  and tested input (correct columns,GrangesList)
  tmp <- checkGR(GRangesO=benchData,reg.definition=regCol)
  if (!any(colnames(mcols(reg2Gene)) %in% regCol)){stop("No regCol detected")}
  
  
  
  # run associateReg2Gene for gene chunks and save them separately
  
  
  if (class(benchData)=="GRangesList"){
    
    
    BenchmarkingResults <- filterPreBenchList(reg2Gene=reg2Gene,
                                              benchList=benchData,
                                              regCol=regCol,
                                              strandness=strandness)
    
    if (reportGR==TRUE) {mcols(reg2Gene) <- c(mcols(reg2Gene),
                                              BenchmarkingResults)
    return(reg2Gene)}
    
    if (reportGR==F) {return(BenchmarkingResults)}  
    
  }  
  
  if (class(benchData)=="GRanges"){
    
    
    BenchmarkingResults <- filterPreBenchSimple(reg2geneAssoc=reg2Gene,
                                                benchmarkData=benchData,
                                                regCol=regCol,
                                                strandness=strandness)
    
    
    
    if (reportGR==TRUE) {reg2Gene$Bench <- BenchmarkingResults
    return(reg2Gene)}
    if (reportGR==F) {return(BenchmarkingResults)}  
  }
} 



#' Function that prefilters GRanges object prior benchmarking procedure 
#' when benchmark is GRanges object 
#' 
#' @param reg2geneAssoc a GRanges object output (from \code{associatereg2Gene}).
#' GRanges object should have at least one meta-data column which stores
#' the second GRanges object.
#' @param benchmarkData a Granges object with at least one meta-data 
#' column which stores the second GRanges object named "reg". Both regions
#' are used in benchmarking procedure.
#' @param regCol character (default "reg"), a column name of meta-data object.
#' Indicates the column where the 2nd pair of GRanges is stored,
#' @param nCores possible to be runned in parallel
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
#' @author Inga Patarcic
#' @keywords internal
filterPreBenchSimple <- function(reg2geneAssoc,
                                 benchmarkData,
                                 regCol="reg",
                                 strandness=FALSE){
  
  
  # setting both Benchmark coordinates on top of one another
  benchmark2reg <- c(benchmarkData,switchReg(benchmarkData,regCol))
  
  # checking if region1 or region2 overlap with any region from the
  # benchmark dataset (granges pooled together)
  reg2geneFilterC1 <- findOverlaps(reg2geneAssoc,benchmark2reg,
                                   ignore.strand=strandness)
  reg2geneFilterC2 <- findOverlaps(switchReg(reg2geneAssoc,regCol),
                                   benchmark2reg,
                                   ignore.strand=strandness)
  
  
  # identifying rows of reg2gene object which are present in both the 
  # overlap with benchmarking region1 and benchmarking region2 but not 
  # necessarily they overlap to the same benchmark, only prerequisite is
  # that rows are detected in the both gene coordinates
  reg2geneFilterRow1 <- unique(data.frame(reg2geneFilterC1)$queryHits)
  reg2geneFilterRow2 <- unique(data.frame(reg2geneFilterC2)$queryHits)
  
  # identify rows present in both overlapping procedure
  reg2geneFR <- reg2geneFilterRow1[reg2geneFilterRow1%in%reg2geneFilterRow2]
  
  # reporting vector of zeros&ones
  
  reg2geneFiltered <- rep(0,length(reg2geneAssoc))
  reg2geneFiltered[reg2geneFR] <- 1
  
  return(reg2geneFiltered)
  
  
}

#' Function that prefilters GRanges object prior benchmarking procedure 
#' when benchmark is GRangesList
#' 
#' 
#' @param reg2Gene a GRanges object output (from \code{associatereg2Gene}).
#' GRanges object should have at least one meta-data column which stores
#' the second GRanges object.
#' @param benchList a GrangesList object with at least one meta-data 
#' column which stores the second GRanges object named "reg". Both regions
#' are used in benchmarking procedure.
#' @param regCol character (default "reg"), a column name of meta-data object.
#' Indicates the column where the 2nd pair of GRanges is stored,
#' @param nCores possible to be runned in parallel
#' @param strandness When set to TRUE, the strand information is ignored in 
#' the overlap calculations. ignore.strand argument in 
#' \code{\link[GenomicRanges]{findOverlaps}}
#' @author Inga Patarcic
#' @keywords internal
filterPreBenchList <- function(reg2Gene,
                               benchList,
                               regCol="reg",
                               nCores=1,
                               strandness=FALSE,
                               ...){
  
  # running benchmarking per benchmark dataset
  bench_list <- mclapply(benchList, function(x) {
    
    filterPreBenchSimple(reg2geneAssoc=reg2Gene,
                         benchmarkData=x,
                         regCol = regCol,
                         strandness=strandness)},
    mc.cores = nCores)
  
  bench_list <- do.call("cbind.data.frame",bench_list)
  
  return(bench_list)
  
  
}


#' the function checks if the input GRanges in the GRangesList have
#' the right columns, it should have at least a reg.definition column as a GRanges 
#' object returns TRUE if GRanges has a reg column pointing to a GRanges and
#' if associations is a GRangesList object
#' @param GRangesO a named \code{GRangesList} where each element is a set of
#' associations returned by \code{\link{associateReg2Gene}} or
#' \code{\link{metaAssociations}}.
#' @param reg.definition column where 2nd GRanges object is stored
#' @keywords internal
#' @author Inga P
checkGR<-function(GRangesO,reg.definition){
  
  lapply(GRangesO,function(x){
    
    is.reg=which(colnames(mcols(x)) %in% reg.definition)
    
    if (!length(is.reg)) {stop("No regCol detected")}
    
    if (class(mcols(x)[is.reg][,1]) != "GRanges") {stop("regCol not GRanges")}
  })
  
}















