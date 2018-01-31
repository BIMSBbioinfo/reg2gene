#' Help f() to extract info about cohort/method/algorithm
#' 
#' @author Inga Patarcic
#' 
#' @keywords internal
extractInfo <- function(path,
                        typeOptions=method){
  
  return(str_extract(path,
                     paste0(typeOptions,collapse = "|")))
  
  
} # "algorithm", "cohort"



#' Creates Index table of files produced by reg2gene data
#' integration and modelling
#' 
#' The function creates a data frame with the following columns:
#' 1) path to .rds files - results produced by modelling or data integration 
#' re2gene f() and corresponding info:
#' 2) methods - which method type is used for this analysis: eg. H3K4me1,...
#' 3) algorithms - which algorithm type is used for this analysis: eg. dcor,...
#' 4) cohorts - which cohort type is used for this analysis: eg. Roadmap,...
#' 5) type - define whether results are produced by voting, meta-analysis,from 
#' individual modelling  or a result of regActivityAroundTSS() 
#' 
#' @param pathModels (character) a path to folder where .rds files are stored or 
#' an object with paths to .rds files.  
#' # IMPORTANT! 
#' cohort/method/algorithm names should be indicated in the path, 
#' or .rds files with this name should be present in the defined directory 
#' @param type "ind" (default; character)
#' OPTIONS: ("votingAlgoritm","votingCohorts","votingMethods","metaA","regAct",
#' "ind").
#' Define whether results are produced by voting by("votingAlgoritm",
#' "votingCohorts","votingMethods"), meta-analysis("metaA"), or from individual 
#' modelling ("ind") or a result of regActivityAroundTSS() ("regAct"). 
#' Important for analyses when all these info is gathered together
#' @param method "H3K4me1" (default; character). A character vector of all 
#' methods used for this analysis. 
#' Written for the option: method=c("H3K4me1","H3K27ac","Methylation","DNase")
#' @param algorithm "pearson" (default; character). A character vector of all 
#' methods used for this analysis. 
#' Written for the option: algorithm=c("pearson","spearman","elasticnet","dcor",
#' "randomForest")
#' @param cohort "Roadmap" (default; character). A character vector of all 
#' cohorts used for this analysis. 
#' Written for the option: cohort= c("Roadmap","Blueprint","CEMT","McGill")
#' @param ... maybe to add in the future
#' 
#' @return
#' A data frame with the following columns (n=5):
#' 1) path to .rds files - results produced by modelling or data integration 
#' re2gene f() and corresponding info:
#' 2) methods - which method type is used for this analysis: eg. H3K4me1,...
#' 3) algorithms - which algorithm type is used for this analysis: eg. dcor,...
#' 4) cohorts - which cohort type is used for this analysis: eg. Roadmap,...
#' 5) type - define whether results are produced by voting("votingAlgoritm",
#' "votingCohorts","votingMethods"), meta-analysis("metaA"),  from individual 
#' modelling ("ind") or a result of regActivityAroundTSS() ("regAct")
#' 
#' @details  An IndexTable which is a result of makeIndexTable() is used as
#' an input for following functions: selectEP(),votedEP(),plotGEEA() to
#' gathers paths to different reg2gene data integration and modelling results,
#' and allows easier downstream analysis. This f() should be separately runned
#' for voting("votingAlgoritm","votingCohorts","votingMethods"), 
#' meta-analysis and individual modelling, and regActivityAroundTSS() results, 
#' and individual IndexTables should be pooled together, such that plotGEEA()
#' function can be runned (it requires paths to "regAct" .rds files - results of 
#' regActivityAroundTSS() and for example modelling paths.
#' 
#' @author Inga Patarcic
#' 
#' @examples require(stringr)
#' 
#' makeIndexTable(pathModels="~/H3K4me1pearsonRoadmap.rds",
#' type="ind",method="H3K4me1",algorithm="pearson",cohort="Roadmap")
#' 
#' # cohort/method/algorithm names should be indicated in the path, or .rds 
#' # files with this name should be present in the defined directory
#' 
#' makeIndexTable(pathModels="~/test.rds",
#' type="ind",method="H3K4me1",algorithm="pearson",cohort="Roadmap")
#' 
#' 
#' # HOW-TO-USE-IT EXAMPLE
#' 
#' IndexTable1 <- makeIndexTable(pathModels="~/H3K4me1pearsonRoadmap.rds",
#' type="ind",method="H3K4me1",algorithm="pearson",cohort="Roadmap")
#' 
#' IndexTable2 <- makeIndexTable(pathModels="~/ra_H3K4me1pearsonRoadmap.rds",
#' type="regAct",method="H3K4me1",algorithm="pearson",cohort="Roadmap")
#' 
#' IndexTable <- rbind(IndexTable1,IndexTable2) # IndexTable use downstream  
#' 
#' \dontrun{
#' algorithm <- c("pearson","spearman","elasticnet","dcor","randomForest")
#' cohort=c("Roadmap","Blueprint","CEMT","McGill")
#' method=c("H3K4me1","H3K27ac","Methylation","DNase")
#' 
#' votedPath="/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD/"
#' metaPath="/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/MetaA/"
#' pathRegActTSS="/data/akalin/Projects/AAkalin_reg2gene/Results/regActivityAroundTSS/FishilevichChromHMMoneTOone/"
#' pathModels="/data/akalin/Projects/AAkalin_Catalog_RI/Data/ValidationDataset/Fishillevich/ModellingResults/"
#' votedAlg="/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD/VotingAlg/"
#' votedC="/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD/VotingCohort/"
#' votedM="/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD/VotingMethod/"
#' 
#' IndexTable <- data.frame(rbind(makeIndexTable(metaPath,"metaA"),
#'                                makeIndexTable(votedAlg,"votingAlgorithm"),
#'                                makeIndexTable(votedC,"votingCohorts"),
#'                                makeIndexTable(votedM,"votingMethods"),
#'                                makeIndexTable(pathRegActTSS,"regAct"),
#'                                makeIndexTable(pathModels,"ind")),stringsAsFactors = F)
#' saveRDS(IndexTable,"/data/akalin/Projects/AAkalin_Catalog_RI/Data/ValidationDataset/Fishillevich/IndexTable.rds")
#' 
#' }
#' 
#' @export
makeIndexTable <- function(pathModels,
                           type="ind",
                           method="H3K4me1",
                           algorithm="pearson",
                           cohort="Roadmap"){
  
  
  #1.  OBTAIN modelling Results
  
  if (!str_detect(pathModels,"rds")){
    
    pathModels <- list.files(pathModels, full.names = T)
    
  }
  
  # extracting info about algorithm/method/cohort
  
  return(cbind(path=pathModels,
               methods=extractInfo(path=pathModels,typeOptions=method) ,
               algorithms=extractInfo(path=pathModels,typeOptions=algorithm),
               cohorts=extractInfo(path=pathModels,typeOptions=cohort),
               type=type))
  
}




#' Returns N TOP/BOTTOM/RANDOM associations per 
#' cohort/method/algorithm/voting/metaA
#' 
#' The function returns GInteractions object with selected N TOP/BOTTOM/RANDOM 
#' associations per cohort/method/algorithm/voting/metaA. It can return
#' N random statistically not significant associations as well.
#' 
#' 
#' @param indexTable an output of the makeIndexTable() for given combination of
#' methods/algorithms/cohorts; a table should contain the following columns: 
#' path - paths to .rds files; methods - a method type used for the analysis: 
#' eg. H3K4me1; 
#' algorithms -  algorithm used for the analysis: eg. dcor, cohorts - cohort
#' used for the analysis: eg. Roadmap; type - define whether results are 
#' produced by voting, meta-analysis,from individual modelling or a result of
#' regActivityAroundTSS()[OPTIONAL]
#' @param topN "all" (default) or integer. Defines how many EP associations to 
#' report. If topN=="all" then all associations with predefined statistics
#' [typeSuccess<thresholdSuccess],  eg pval<0.05 are reported. If topN is an
#' integer then corresponding number of EP associations 
#' (as defined by select and success arguments) is reported. 
#' @param select "top" (default). Other options "bottom" and "random". Indicates
#' whether to select top N genes [topN], bottom N genes [topN] or randomly 
#' select N genes [topN]
#' @param success TRUE (default). Whether to select statistically significant EP 
#' associations (succes==TRUE), or those that were found to be statistically
#' not significant EP associations (succes==FALSE). 
#' @param typeSuccess "pval" (default) Other options "qval". Which statistics to
#' threshold to assess a statistical significance.
#' @param thresholdSuccess 0.05 (default, numeric). A threshold useD to 
#' assess a statistical significance.
#' @param method "H3K4me1" (default; character). For which method to 
#' extract Enh~Promoter pairs. Common options: "H3K4me1","H3K27ac",
#' "Methylation","DNase"
#' @param algorithm "pearson" (default; character). For which algorithm to
#' extract Enh~Promoter pairs. Common options:"pearson","spearman",
#' "elasticnet","dcor","randomForest"
#' @param cohort "Roadmap" (default; character). For which cohort  to extract
#' Enh~Promoter pairs. Common options:"Roadmap","Blueprint",
#' "CEMT","McGill"
#' @param metaA FALSE (default). Whether to perform this analysis on results of
#' meta-analysis or not. This argument overwrites cohort argument.
#' 
#' @return  
#' returns a \code{\link[InteractionSet]{GInteractions}} object filtered 
#' for N (or retained all) enhancer-promoter interactions. Corresponding 
#' statistics from either modelling analysis \code{\link{associateReg2Gene}} or 
#' meta-analysis \code{\link{metaAssociations}} is retained. 
#' 
#' @details  
#' This function allows easy export of TOPN gene~enhancer pairs per method,
#' cohort and algorithm combination of modelling results or meta-analysis.
#' It can return all statistically significant/unsignificant genes or 
#' select top, bottom or random N of genes. Which statistics to filter can be
#' choosen by user (typeSuccess), as well as a threshold value for it 
#' (thresholdSuccess).
#' 
#' @author Inga Patarcic
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' require(stringr)
#' library(InteractionSet)
#' 
#' # ONE .rds file will need to be added into package!!!
#' 
#' IndexTable <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD//CohortVoting_McGillH3K4me1.rds")
#' 
#' 
#' TOPgenes <- selectEP(IndexTable,
#' topN=16, #change argument "all" not all
#' select="top", #"bottom","random",
#' success=TRUE,
#' typeSuccess="pval",
#' thresholdSuccess=0.05,
#' method="H3K4me1",
#' algorithm="dcor",
#' cohort="Roadmap",
#' metaA=FALSE)
#' 
#' saveRDS(TOPgenes,"/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevich/TOPgenesH3K4me1dcorRoadmap.rds")
#' }
#' 
#' 
#' 
#' @export
selectEP <- function(indexTable,
                  topN="all", #change argument "all" not all
                  select="top", #"bottom","random",
                  success=TRUE,
                  typeSuccess="pval",
                  thresholdSuccess=0.05,
                  method="H3K4me1",
                  algorithm="pearson",
                  cohort="Roadmap",
                  metaA=FALSE){ # "votingAlgoritm","votingCohorts","votingMethod","metaA"
  
  
  # select Individual modelling results if metaA=T
      if (metaA==FALSE){
       
       # testing if selected method,algorithm or cohort present in datasets
            if (!all(method%in%indexTable$methods,
                algorithm%in%indexTable$algorithms,
                cohort%in%indexTable$cohorts)){
              
              stop("Selected method,algorithm or cohort not in the indexTable")
            
              }
      
      # select appropriate .rds file  
       
        SelectedRDS <-  indexTable$path[which((indexTable$methods==method)&
                                            (indexTable$algorithms==algorithm)&
                                            (indexTable$cohorts==cohort))]
        
        ModellingResults <- readRDS(SelectedRDS)  
        
      }
  
  # select Individual modelling results if metaA=T
      if (metaA==TRUE){
      
          if (!all(method%in%indexTable$methods,
                  cohort%in%indexTable$cohorts,
                  "metaA"%in%indexTable$type)){
            
        stop("Selected method,cohort or meta-analysis is not in the indexTable")
            
          }
     
      
      SelectedRDS <-  indexTable$path[which((indexTable$methods==method)&
                                           (indexTable$type=="metaA")&
                                           (indexTable$algorithms==algorithm))]

      
      ModellingResults <- readRDS(SelectedRDS)  
      
    
      
      }
  
  # Filtering & orginizing results
       
        # identifying + extracting column to filter 
         typeSuccessPos <- which(colnames(mcols(ModellingResults))==typeSuccess)
         StatisticsData <- unlist(mcols(ModellingResults)[typeSuccessPos])
        
        # ordering by "pvalue" or "qvalue"
         ModellingResults <- ModellingResults[order(StatisticsData)]
        # filtering by "pvalue" or "qvalue"
         SuccfRows <- StatisticsData[order(StatisticsData)]<thresholdSuccess
   
         
         
        # select based on topN, select and success arguments
           if (success==TRUE) { SuccfRows <- which(SuccfRows)}
           if (success==FALSE) { SuccfRows <- which(!SuccfRows)}  
             
         ModellingResults <- ModellingResults[SuccfRows]
         
         # filtering of EP when topN!="all"
         if (topN!="all"){ 
         
           if (select=="top") { ModellingResults <- ModellingResults[1:topN]}
           
           if (select=="bottom") { 
             ModellingResults <- ModellingResults[tail(length(ModellingResults),
                                                       topN)]}
           
           if (select=="random") { 
             ModellingResults <- ModellingResults[sample(length(ModellingResults
                                                                ),topN)]} 
         }
         
         return(ModellingResults)
        
         
       }






         
#' Returns voted EP pairs for the combination of cohort/method/
#' algorithm and the voting procedure.
#' 
#' This function allows easy import of results of different voting procedures
#' for given combination of cohort/method/algorithm.
#' 
#' @param indexTable an output of the makeIndexTable() for given combination of
#' methods/algorithms/cohorts; a table should contain the following columns: 
#' path - paths to .rds files; methods - a method type used for the analysis: 
#' eg. H3K4me1; 
#' algorithms -  algorithm used for the analysis: eg. dcor, 
#' cohorts - cohort used for the analysis: eg. Roadmap; 
#' type - define whether results are produced by voting, meta-analysis, from
#' individual modelling or a result of regActivityAroundTSS()
#' @param votingType "votingAlgoritm" (default; character) Describes which
#' voting procedure was used to produce voting results.
#' OPTIONS: ("votingAlgoritm","votingCohorts","votingMethods"). 
#' @param method "H3K4me1" (default; character). For which method to 
#' extract Enh~Promoter pairs. Common options: "H3K4me1","H3K27ac",
#' "Methylation","DNase"
#' @param algorithm "pearson" (default; character). For which algorithm to
#' extract Enh~Promoter pairs. Common options:"pearson","spearman",
#' "elasticnet","dcor","randomForest"
#' @param cohort "Roadmap" (default; character). For which cohort  to extract
#' Enh~Promoter pairs. Common options:"Roadmap","Blueprint",
#' "CEMT","McGill"
#'    
#' @details This function allows easy import of results of different voting
#' procedures for given combination of cohort/method/algorithm. For example, 
#' "votingAlgoritm" correspond to voting based on the results of pearson,
#' spearman,elasticnet,dcor,randomForest modelling procedures for a combination
#' of cohort and method (eg H3K4me1+Roadmap). In this case algorithm 
#' argument is ignored while importing "votingAlgoritm" results 
#' for H3K4me1+Roadmap
#' 
#' @examples
#' 
#' \dontrun{
#' library(InteractionSet)
#' library(stringr)
#' IndexTable <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD//CohortVoting_McGillH3K4me1.rds")
#' 
#' votedEP(indexTable,
#' votingType="votingAlgorithm",
#' method="H3K4me1",
#' algorithm="pearson",
#' cohort="Roadmap")
#'  
#' }
#' @author Inga Patarcic
#' @export
votedEP <- function(indexTable,
                     votingType="votingAlgoritm",
                     method="H3K4me1",
                     algorithm="pearson",
                     cohort="Roadmap"){ # "votingAlgoritm","votingCohorts","votingMethod","metaA"
  
  # select Individual modelling results if metaA=T
 
    # testing if selected method,algorithm or cohort present in datasets, 
    # min 3 overlap expected for voting methods
    if (sum(method%in%indexTable$methods,
             cohort%in%indexTable$cohorts,
             algorithm%in%indexTable$algorithm,
             votingType%in%indexTable$type)<3){
    
      stop("Selected method,cohort or votingType is not in the IndexTable")
      
    }
    
    # select appropriate .rds file  
      # 1. select voted method
  indexTable <-  indexTable[which(indexTable$type==votingType),]      
      
      # 2. identify entries for which 2 overlap from meth/alg/cohort
      IndexALLC <- (c(which(indexTable$methods==method),
                    which(indexTable$algorithms==algorithm),
                    which(indexTable$cohorts==cohort)))
      
      # select identified .rds file
    SelectedRDS <-  indexTable$path[IndexALLC[duplicated(IndexALLC)]]
    
    ModellingResults <- readRDS(SelectedRDS)  
    
  return(ModellingResults)
  
  
}

#' Scatterplots of enhancer activity~gene expression for max 16 EP pairs 
#' 
#' For GInteractions object that contains associations between a regulatory 
#' regions and TSSs (EP pairs),estimated association statistic,and corresponding
#' p-values this function read files that contain quantified (and normalized) 
#' enhancer activity and gene expression levels for each data entry (cell type).
#' Then this info is plotted as a scatterplot of enhancer activity~gene 
#' expression for max 16 EP pairs. 
#' 
#' 
#' @param enhPromPairs A GInteractions object of length 16 with the 
#' corresponding statistics as a results of either reg2gene modelling 
#' \code{\link{associateReg2Gene}}, meta-analysis \code{\link{metaAssociations}}
#' or their pre-filtered versions \code{\link{selectEP}}. 
#' It requires information about statistics 
#' ("coef") and corresonding testing results ("pval"). In the future it could 
#' be adjusted for ("qval" or that it does not require any statistics).
#' @param indexTable an output of the makeIndexTable() for given combination of
#' methods/algorithms/cohorts; a table should contain the following columns: 
#' path - paths to .rds files; methods - a method type used for the analysis: 
#' eg. H3K4me1; 
#' algorithms -  algorithm used for the analysis: eg. dcor, cohorts - cohort
#' used for the analysis: eg. Roadmap; type - define whether results are 
#' produced by meta-analysis, from individual modelling or a result of
#' regActivityAroundTSS() [NECESSARY]
#' @param method "H3K4me1" (default; character). For which method to 
#' extract Enh~Promoter pairs. Common options: "H3K4me1","H3K27ac",
#' "Methylation","DNase"
#' @param algorithm "pearson" (default; character). For which algorithm to
#' extract Enh~Promoter pairs. Common options:"pearson","spearman",
#' "elasticnet","dcor","randomForest"
#' @param cohort "Roadmap" (default; character). For which cohort  to extract
#' Enh~Promoter pairs. Common options:"Roadmap","Blueprint",
#' "CEMT","McGill"
#' 
#' @details This function allows quick visualization of selected EP pairs as a
#' scatterplot of corresponding levels of enhancer activity and gene expression 
#' for all data entries (cell types used in the modelling procedure:
#' \code{\link{associateReg2Gene}}). This info is plotted for max 16 EP pairs. 
#' It is useful in the case when one wants to visualize the background levels of 
#' enhancer activity and gene expression for TOP EP associations. Since, for 
#' each combination of cohort/methods/algorithms separated modelling is 
#' performed one needs to define a combination.
#' 
#' @author Inga Patarcic
#' 
#' @examples 
#' 
#' require(ggplot2)
#' require(stringr)
#' 
#' \dontrun{
#' 
#' TOPgenes <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevich/TOPgenesH3K4me1dcorRoadmap.rds")
#' IndexTable <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD//CohortVoting_McGillH3K4me1.rds")
#' 
#' plotGEEA(TOPgenes,
#'        indexTable=IndexTable,
#'        method="H3K4me1",
#'        algorithm="dcor",
#'        cohort="Roadmap")
#' 
#' }
#' @export
plotGEEA <- function(enhPromPairs, # output of selectEP
                        indexTable,
                        method="H3K4me1",
                        algorithm="dcor",
                        cohort="Roadmap") {
  
  
    
  regActRes <- readRDS(indexTable$path[which(((indexTable$methods==method)
                                   &(indexTable$type=="regAct")
                                   &(indexTable$cohort==cohort)))])
  
  
       # creating per gene~enhancer object to plot
       
       #enhPromPairs <- TOPgenes
       
       EP <- lapply(1:length(enhPromPairs),function(x){
        # x=10
         gEeA <- enhPromPairs[x]
       
              regActRes <- unlist(regActRes[names(regActRes)==gEeA$name])
       
           # obtaining info about gene expression
           geneExpression <- regActRes[1]
           geneExp <- unlist(mcols(geneExpression)[!names(mcols(geneExpression))
                                         %in%c("featureType","name","name2")])
             
           # obtaining info about enhancer activity
           enhActivity <- findOverlaps(first(gEeA),
                                       regActRes, 
                                       type="equal",
                                       select="first")
           enhActivity <- regActRes[enhActivity]
           
           enhAct <- unlist(mcols(enhActivity)[!names(mcols(enhActivity))
                                               %in%c("featureType","name",
                                                     "name2")])
           
           
           ModResult <- paste0(c("pval=", signif(gEeA$pval,digits=2),
                               " coef=", signif(gEeA$coefs,digits=2)),
                               collapse = "")
           
          EPpairs <- paste(geneExpression$name,enhActivity)
           
          return(cbind(geneExp,enhAct,ModResult,EPpairs))
       })
  
       EP  <- do.call("rbind.data.frame",EP)
     
      
   # resolve plotting problems
    
      ggplot(EP, aes(as.numeric(as.character(geneExp)), 
                           (as.numeric(as.character(enhAct))))) +
          geom_point() +
          facet_wrap(~ EPpairs + ModResult , 
                     nrow = 4,ncol=4, scales = "free")  + 
          coord_equal() +
          labs(x = "Gene Expression", y = "Enh Activity") +
          ggtitle(paste0(c("Details:", method,algorithm,cohort),collapse = ""))
  
     
}




#' Reports a number of statistically significant enhancer~
#' promoter associations
#' 
#' compareModelStat() reports the number of statistically significant enhancer~
#' promoter associations for voting or meta-analysis and corresponing individual
#' modelling results. IMORTANT NOTE: since it uses results of previously runned
#' voting procedure it is important that statistics and threshold values used in
#' voting procedure are the same as used in this algorithm. However,
#' meta-analysis statistics thresholding is performed as a part of this function
#' ,thus statistics and threshold [typeSuccess&thresholdSuccess] values are used
#' to filter both: meta-analysis results and individual modelling results.
#'  
#' 
#' 
#' @param indexTable an output of the makeIndexTable() for given combination of
#' methods/algorithms/cohorts; a table should contain the following columns: 
#' path - paths to .rds files; methods - a method type used for the analysis: 
#' eg. H3K4me1; 
#' algorithms -  algorithm used for the analysis: eg. dcor, 
#' cohorts - cohort used for the analysis: eg. Roadmap; 
#' type - define whether results are produced by voting or meta-analysis, 
#' AND individual modelling (indexTable$type=="ind"). 
#' It is necessary that individual modelling paths are
#' present in the index table (indexTable$type=="ind"), because results of these
#' models are compared to the results of voting or meta-analysis procedure.
#' 
#' @param type "votingAlgorithm" (default; character). Describes which voting
#' procedure is used to produce voting results (common options:"votingAlgoritm",
#' "votingCohorts","votingMethods") or if meta-analysis is used (type="metaA").
#' This argument corresponds to the type argument in 
#' \code{\link{makeIndexTable}}, thus if some other type character object is
#' used with \code{\link{makeIndexTable}} then the same argument should be used
#' @param typeSuccess "pval" (default) Other options "qval". Which statistics to
#' threshold to assess a statistical significance.
#' @param thresholdSuccess 0.05 (default, numeric). A threshold useD to 
#' assess a statistical significance.
#' @param method "H3K4me1" (default; character). For which method to 
#' extract Enh~Promoter pairs. Common options: "H3K4me1","H3K27ac",
#' "Methylation","DNase"
#' @param algorithm "pearson" (default; character). For which algorithm to
#' extract Enh~Promoter pairs. Common options:"pearson","spearman",
#' "elasticnet","dcor","randomForest"
#' @param cohort "Roadmap" (default; character). For which cohort  to extract
#' Enh~Promoter pairs. Common options:"Roadmap","Blueprint",
#' "CEMT","McGill"
#' 
#' 
#' @return  
#' returns a dataframe with information about the number of identified 
#' statistically enhancer~promoter associations and info about corresponding
#' input modelling datasets used - cohort,method,algorithm combination.
#' The same statistics is returned for the requested voting method or 
#' meta-analysis. 
#' 
#' 
#' @details  
#' This function allows an easy assesment of the success of voting or 
#' meta-analysis procedure compared to individual modelling procedures, through 
#' returning the number of statically significant cases reported in voting or 
#' meta-analysis and corresponding individual modelling procedures.
#' 
#' @author Inga Patarcic
#' 
#' @examples
#' library(stringr)
#' library(InteractionSet)
#' \dontrun{
#' 
#' IndexTable <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Validation/Fishillevic/VoteD//CohortVoting_McGillH3K4me1.rds")
#' 
#' compareModelStat(indexTable=IndexTable, 
#' type="votingAlgorithm",
#' method="H3K4me1", 
#' algorithm="dcor",
#' cohort="Blueprint")
#' 
#' }
#' 
#' @export
compareModelStat <- function(indexTable=IndexTable,
                    type="votingAlgorithm",
                    typeSuccess="pval",
                    thresholdSuccess=0.1,
                    method="H3K4me1",
                    algorithm="dcor",
                    cohort="Roadmap"){
  
  # extracting info about path for the type of interest: voting or metaA
  
      indexTable.ss <- indexTable[is.na((indexTable$type==type)&
                           (indexTable$method==method)&
                           (indexTable$algorithm==algorithm)&
                           (indexTable$cohort==cohort)),]
  
  # 1. obtaining statistics for type of interest: voting or metaA
  if (type!="metaA"){
      
    # obtaining voting statistics
          positivesN <- length(votedEP(indexTable.ss,
                         votingType=type,
                         method=method,
                         algorithm=algorithm,
                         cohort=cohort))
  
                        }
  if (type=="metaA"){
    
     # obtaining metaA statistics
        positivesN <- length(selectEP(indexTable.ss,
                          topN="all", 
                          success=TRUE,
                          typeSuccess=typeSuccess,
                          thresholdSuccess=thresholdSuccess,
                          method=indexTable.ss$methods,
                          algorithm=indexTable.ss$algorithms,
                          cohort=indexTable.ss$cohorts,
                          metaA=TRUE))
                        }
      
      
      Pooled <- cbind(indexTable.ss[c('methods',"algorithms","cohorts","type")],
                      positivesN)
      
  # 2. obtaining Info from individual modelling procedures
      
    # get info about which cohort, method, algorithm to used in ind analyses 
    # based on which voting has been done
      
        MAC <- unlist(indexTable.ss[c('methods',"algorithms","cohorts")])
        MAC <- MAC[complete.cases(MAC)]
      
    # extract info about individual datasets
    
        indexTable.IND <- indexTable[(indexTable$type=="ind")&
                                    unlist(indexTable[names(MAC[1])]==MAC[1])&
                                    unlist(indexTable[names(MAC[2])]==MAC[2]),]
        
    # obtaining individual model - statistics
    
        positivesN <-  sapply(1:nrow(indexTable.IND), function(x){
            
                          indexTable.IND.ss <- indexTable.IND[x,]
            
                                  return(length(selectEP(indexTable.IND.ss,
                                           topN="all", 
                                           success=TRUE,
                                           typeSuccess=typeSuccess,
                                           thresholdSuccess=thresholdSuccess,
                                           method=indexTable.IND.ss$methods,
                                           algorithm=indexTable.IND.ss$algorithms,
                                           cohort=indexTable.IND.ss$cohorts)))
                      })
        
         
        Ind <- cbind(indexTable.IND[c('methods',"algorithms","cohorts","type")],
                   positivesN)
        
      
    # pool info from voting/metaA and individual models together    
      return(rbind(Ind,Pooled))
  
    }
    
