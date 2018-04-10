#' Function that rearranges list of scores calculated per cell type
#' 
#' @param scores.per.cell.type list of result score matrices for each biwgig
#'  file
#' @param sampleIDs  A vector of unique sample ids/names of bigwig files 
#' @return a dataframe of quatified scores per genome window and bigwig files.
#' Each row corresponds to ranges in the genome, and columns correspond to the 
#' bigwig files. 
#' @keywords internal
#' @author Inga Patarcic
scoresAsMcols <- function(scores.per.cell.type,sampleIDs){
  
  
          mcols.per.cell.type <- plyr::rbind.fill(lapply(scores.per.cell.type,
                      function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))
     
        # adjust column names - remove bw extension
        rownames(mcols.per.cell.type) <- sampleIDs
       
     return(t(mcols.per.cell.type))
}



#' Function that runs normalization for the quantified activity/expression level
#' 
#' @param scoresDF a dataframe of quatified scores per genome window and bigwig 
#'  file
#' @param normalizationMet  charater. which normalization procedure is wanted.
#' @return a normalized dataframe of quatified scores per genome window and 
#' bigwig files. 
#' @keywords internal
#' @author Inga Patarcic
normalizeScores <- function(scoresDF,normalizationMet){

    if (normalizationMet=="ratio") {   
            # estimate size factors and multiply gene expression using these 
            # factors
    sizeFactors <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(scoresDF))
    scoresDF.normalized <- as.matrix(scoresDF)*rep(sizeFactors,each=
                                                            nrow(scoresDF))}
                
    if (normalizationMet=="quantile") {  
 scoresDF.normalized <- preprocessCore::normalize.quantiles(as.matrix(scoresDF))
          rownames(scoresDF.normalized) <- rownames(scoresDF)
          colnames(scoresDF.normalized) <- colnames(scoresDF)}
           
 
  return(scoresDF.normalized)

  
  }
  



#' Function that quantifies exon expression based on strand specificity
#' 
#'  Function that separates stranded and unstranded libraries and quantifies 
#'  exon expression separately for stranded libraries eg expression of exons 
#'  located on the + strand is quantified using forward libraries if libraries
#'  are unstranded then exon orientations is not taken into account
#'  
#' @param exons GRanges object with exon coordinates and metadata about genes
#' @param exons.splitted GRanges object with exon coordinates splitted based on
#' the strand.
#' @param geneActSignals character, paths to bigwig files
#' @param libStrand character, library strandness ("+","-","*")
#' @param sampleIDs character, name of samples/bigwig files
#' @param summaryOperation argument for ScoreMatrixBin() how to quantify scores
#' @param mc.cores integer, number of cores for mclapply
#' @return Granges object with metadata which corresponds to quantified scores 
#' per input exon regions and bigwig files
#' @keywords internal
#' @author Inga Patarcic
exonExpressionStrandAdjusted <- function(exonsForQuant,
                                          exonsSplitted,
                                          geneActSignals,
                                          sampleIDs,
                                          libStrand,
                                          mc.cores,
                                          summaryOperation) {
  
  #exonsForQuant <- exons
  #exonsSplitted <- exons.splitted

  
  if (length(exonsSplitted$`*`)!=0){stop("Strand of exons can be only + or -")}
  
        # detect stranded tracks and unstranded tracks    
        ForwardLibraries <- stringr::str_detect(libStrand,"\\+") 
        Unstranded.libraries <- stringr::str_detect(libStrand,"\\*")
        
  if (sum(ForwardLibraries)!=(sum(!(ForwardLibraries|Unstranded.libraries)))){
    stop("Number of forward and reverse libraries do not match")}
  
  
  # separate analysis for unstranded & stranded libraries
  
  if (sum(Unstranded.libraries)!=length(libStrand)){
    # test strandness - "mixed"&"stranded"
    # Stranded libraries - positive/forward strand
    
    # forward strand
    exonscoresForwardStrand <- mclapply(geneActSignals[ForwardLibraries], 
                                        function(x) {try(
                                        genomation::ScoreMatrixBin(x,
                                              windows =  exonsSplitted$`+`,
                                              bin.num = 1,
                                              type = "bigWig",
                                              bin.op=summaryOperation, 
                                              is.noCovNA=F),
                                              silent = T)},
                                        mc.cores=mc.cores)
    
    if (max(sapply(exonscoresForwardStrand,class)=="try-error")){
      stop(paste("Error with ScoreMatrixBin and following bigwig files",
                 geneActSignals[ForwardLibraries][sapply(
                   exonscoresForwardStrand,class)=="try-error"]))}                                      
    
    exonscoresForwardStrand.df <- scoresAsMcols(scores.per.cell.type=
                                                  exonscoresForwardStrand,
                                        sampleIDs=sampleIDs[ForwardLibraries])
    
    # reverse strand               
    exonscoresReverseStrand <- mclapply(geneActSignals[!(ForwardLibraries|
                                                        Unstranded.libraries)], 
                                        function(x) {try(
                                          genomation::ScoreMatrixBin(x,
                                               windows =  exonsSplitted$`-`,
                                               bin.num = 1,
                                               type = "bigWig", 
                                               bin.op=summaryOperation,
                                               is.noCovNA=F),
                                               silent = T)},
                                        mc.cores=mc.cores)
    
    if (max(sapply(exonscoresReverseStrand,class)=="try-error")){
      stop(paste("This bigwig files do not overlap exon regions",
                 geneActSignals[!(ForwardLibraries|Unstranded.libraries)]
                 [sapply(exonscoresReverseStrand,class)=="try-error"]))}
    
    exonscoresReverseStrand.df <- abs(scoresAsMcols(scores.per.cell.type=
                                                      exonscoresReverseStrand,
                                        sampleIDs=sampleIDs[ForwardLibraries]))
    
    
    # test if exon length is equal to the length of splitted exons
    if (length(exonsForQuant)!=sum(c(nrow(exonscoresForwardStrand.df),
                             nrow(exonscoresReverseStrand.df)))) {stop("Exons
                                    and splitted exon object are not equal")}
    
    
    # pool reverse and forward strand together    
    
    if (any(c(exonsSplitted$`-`$name,exonsSplitted$`+`$name)!=
            exonsForQuant$name)){break("Exons strand after quantification not correct")}
    
    
    DF_exons <- rbind(exonscoresReverseStrand.df,
                      exonscoresForwardStrand.df)
        row.names(DF_exons) <- NULL # adjusting for problem of equal row.names
        
        values(exonsForQuant) <- cbind(rbind(mcols(exonsSplitted$`-`),
                                             mcols(exonsSplitted$`+`)),
                                   DataFrame(DF_exons))
    
        
        
    }
  
  if (sum(Unstranded.libraries)!=0){
    
    
    # Unstranded
    exonscoresUnstranded <- mclapply(geneActSignals[Unstranded.libraries], 
                                     function(x){try(
                                       genomation::ScoreMatrixBin(x,
                                                      windows = exonsForQuant,
                                                      bin.num = 1,
                                                      type = "bigWig",
                                                      bin.op=summaryOperation,
                                                      is.noCovNA=F),
                                                    silent = T)},
                                     mc.cores=mc.cores)
    
    if (max(sapply(exonscoresUnstranded,class)=="try-error")){
      stop(paste("This bigwig files do not overlap exon regions",
                 geneActSignals[Unstranded.libraries]
                 [sapply(exonscoresUnstranded,class)=="try-error"]))}                                      
    
    
    exonscoresUnstranded.df <- scoresAsMcols(scores.per.cell.type=
                                               exonscoresUnstranded,
                                    sampleIDs=sampleIDs[Unstranded.libraries])                  
    
    values(exonsForQuant) <- cbind(mcols(exonsForQuant),
                                  DataFrame(exonscoresUnstranded.df))
    
    
  }
  
  return(exonsForQuant)                
}



#' Function that quantifies gene expression based on prequantified exon exp
#' 
#' Function that calculates per gene expression based on exon expression results
#' sum across exons (mean value across exon*exon width)/
#' gene_length(sum_of_all_exon_lengths)
#'  
#' @param expressionPerGene GRanges object with exon coordinates and 
#' metadata info about corresponding genes and quantifed exon expression
#' @return Granges object of genes with metadata which corresponds to quantified
#' gene expression
#' @keywords internal
#' @author Inga Patarcic
quantifyGeneExpression <- function(expressionPerGene,exonsForGene){
  
  ExonLength <- width(expressionPerGene)
  GeneLength <- sum(ExonLength,na.rm=T)
  ExpressionPerSample=as.data.frame(mcols(expressionPerGene)[-(1:exonsForGene)])
        
    GeneExpression <- apply(ExpressionPerSample*ExonLength,
                            2,sum,na.rm=T)/GeneLength
  
      return(GeneExpression)
}




#' Calculates regulatory activity over pre-defined regions
#' 
#' The function calculates regulatory activity from histone
#' modification, DNAse or methylation signals for pre-defined regulatory
#' regions and returns a GRanges object with regulatory region locations
#' and their activity over a set of samples.
#' 
#' 
#' @param regRegions a GRanges object that contains regulatory regions
#' over which the regulatory activity will be calculated. 
#' It is strongly suggested to adjust 
#' seqlengths of this object to be equal to the seqlenghts Hsapiens from the 
#' appropriate package (BSgenome.Hsapiens.UCSC.hg19 or whatever needed version).
#' @param activitySignals a named list of BigWig files. Names correspond to 
#'        unique sample ids/names.
#' @param sampleIDs NULL (default). A vector of unique sample
#' ids/names(.bw files), ordered as the bigwig files are ordered. When NULL
#' basenames of .bw files is used as a unique sample ids/names.
#' @param isCovNA (def:FALSE), if this is set to TRUE, uncovered
#' bases are set to NA, this is important when dealing with methylation
#' data, where uncovered bases on a bigWig file do not mean zero methylation.
#' @param weightCol a numeric column in meta data part used as weights. Useful
#' when genomic regions have scores other than their coverage values, 
#' such as percent methylation, conservation scores, GC content, etc.
#' @param normalize NULL(default). Optional "quantile" and "ratio"
#' If set to "quantile" activity measures are quantile normalized as 
#' implemented in \code{\link[preprocessCore]{normalize.quantiles}} and
#' returned ; if set to "ratio" then "median ratio method" implemented
#' as \code{\link[DESeq2]{estimateSizeFactorsForMatrix}} is used to normalize 
#' activity measures.
#' @param summaryOperation "mean"(default).  An argument for 
#' \code{\link[genomation]{ScoreMatrixBin}} that is in the nutshell of 
#' quantifying enhancer activity across pre-defined enhancer regions.
#' This designates which summary operation should be used over the regions.
#' Currently, only mean is available, but "median" or "sum" will be implemented
#' in the future.
#' @param mc.cores (def:1) Define the number of cores to use;
#' at most how many child processes will be run simultaneously 
#' using mclapply from parallel package. Parallelization requires at 
#' least two cores.
#' 
#' @return a GRanges object where its meta-columns correspond
#'         to calculated acvitity measures and column names 
#'         correspond to provided sample ids or names.
#' 
#' 
#' @details regulatory activity is measured by averaging logFC for
#' histone modification ChIP-seq profiles, or DNAse signal, or methylation
#' per base.Currently, relevant bigWig files are required to calculate activity       
#' activity. This function might be extended to work with BAM files
#' in the future. 
#' 
#' @examples # INPUT1: defining .bw files:
#' 
#' test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
#' test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
#' 
#' # INPUT2: defining regulatory regions:regRegions
#' 
#' regRegions <- GRanges(c(rep("chr1",4),rep("chr2",2)),
#'                       IRanges(c(1,7,9,15,1,15),c(4,8,14,20,4,20)),
#'                                             c(rep("+",3),rep("-",3)))
#' regRegions$reg <-  regRegions[c(1,1,3:6)]
#' regRegions$name2 <- regRegions$name <- paste0("TEST_Reg",
#'                                         c(1,1,3:length(regRegions)))
#'                         
#'  # OUTPUT  regActivity():                                                      
#' regActivity(regRegions,c(test.bw,test2.bw))   
#' 
#' # Additionaly, sample names can changed with sampleIDs argument
#' 
#' regActivity(regRegions,c(test.bw,test2.bw),sampleIDs=c("Cell1","Cell2"))
#' 
#' # Additionaly, it supports different normalization procedures:
#' 
#' regActivity(regRegions,c(test.bw,test2.bw),normalize ="ratio")
#' regActivity(regRegions,c(test.bw,test2.bw),normalize ="quantile")
#' 
#' @export
regActivity <- function(regRegions,
                        activitySignals,
                        sampleIDs=NULL,
                        isCovNA=FALSE,
                        weightCol=NULL,
                        summaryOperation="mean",
                        normalize=NULL,
                        mc.cores=1){
  
  
       # test input - ranges
  if (!exists("regRegions")) {stop("regRegions object missing")}
      # test input - bigwig files
  if (!min(sapply(activitySignals,file.exists))) 
          {stop("activitySignals object missing")}
  if (is.null(sampleIDs)) {
      bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
      sampleIDs <- stringr::str_replace(basename(activitySignals),
                               paste(bw.exts,collapse="|"),"")}
          
          
  # calculating coverage
         
          scores.per.cell.type <- parallel::mclapply(activitySignals,
                                    genomation::ScoreMatrixBin,
                                           windows = regRegions,
                                           bin.num = 1,
                                           type = "bigWig", 
                                           is.noCovNA=isCovNA,
                                           weight.col=weightCol,
                                           bin.op=summaryOperation,
                                           mc.cores = mc.cores)
          
         
           # adding scores as mcols  
        scores.per.cell.type.df <- scoresAsMcols(scores.per.cell.type=
                            scores.per.cell.type,sampleIDs=sampleIDs)
                                                   
        
          
          # add normalization
      if (!is.null(normalize)){
        scores.per.cell.type.df <- normalizeScores(scores.per.cell.type.df,
                                                   normalize)}

          mcols(regRegions) <- scores.per.cell.type.df
           
          return(regRegions)
          
  
  
}

#' Identify regulatory regions around provided TSSes
#' 
#' The function identifies the regulatory regions around provided
#' TSSes over a pre-defined window. The function needs a GRanges
#' object for TSSes with meta-columns corresponding to expression
#' levels in different cell types or conditions. 
#' 
#' @param regActivity a GRanges object output from \code{regActivity}
#'        function
#' @param tss a GRanges object output from \code{bwToGeneExp} that
#' contains the TSS location and associated gene expression values 
#' per cell type or condition as meta data. Each row should have a
#' "name" and "name2" columns for unique id or name/symbol for the
#' gene which the TSS is associated with. One is Ensembl id 
#' and the other could be used for the gene symbol. Other metadata
#' column names should represent sample names/ids and should
#' match the GRanges object provided via regActivity argument.
#' @param upstream number of basepairs upstream from TSS to look for
#' regulatory regions. default 500kb
#' @param downstream number of basepairs downstream from TSS to look for 
#' regulatory regions. default 500kb
#' @param mc.cores (def:1) Define the number of cores to use;
#' at most how many child processes will be run simultaneously 
#' using mclapply from parallel package. Parallelization requires at 
#' least two cores.
#' @param force (DEFAULT: FALSE) This argument allows to test input enhancers
#' and genes in 1-to-1 manner. Meaning,a gene expression of gene1 is modelled 
#' using enhancer1 enhancer activity, gene2 gene expression is modelled with
#' enhancer2 enhancer activity, etc. Thus, input gene expression GRanges object
#' (tss) needs to have equal length as enhancer activity GRanges object 
#' (regActivity). It overwrites upstream and downstream arguments since 1-to-1
#' relationship is tested. 
#' 
#' @return An output of regActivityAroundTSS() is a GRangesList object that 
#' contains per gene GRanges with location of the corresponding TSS and
#' regulatory regions identified around that gene. Names for the GRangesList 
#' are unique gene ids/names. Metadata for each GRanges object in the 
#' GRangesList represents regulatory activity and gene expression quantified 
#' across a number of samples (.bw files) that have matched IDs in RNA-Seq and
#' CHiP-Seq,DNase-Seq or bisulfite sequencing experiment. 
#' The GRanges objects have the following metadata columns:
#'  1. featureType: either "gene" or "regulatory"
#'  2. name: name/id for gene and enhancers. Gene name could be id from a 
#'  database enhancer name should be in the format as follows "chr:start-end"
#'  3. name2: a secondary name for the feature, such as gene symbol "PAX6" etc. 
#'  not necessary for enhancers could be NA
#'  4. other columns: numeric values for gene expression or regulatory actvity.
#'    Column names represent sample names/ids.
#' 
#' @examples  # INPUT1: CREATING Toy example for quantified gene expression:
#' 
#' geneExpression <- GRReg1_toy
#'   geneExpression$bw1 <- rep(1,length(GRReg1_toy))
#'   geneExpression$bw2 <- rep(2,length(GRReg1_toy))
#'   geneExpression$bw3 <- rep(3,length(GRReg1_toy))
#'   
#'   # INPUT2: CREATING TOY example for quantified enhancer activity:
#'   
#' regRegion <- GRReg2_toy
#'    regRegion$bw1 <- rep(3,length(regRegion))
#'    regRegion$bw2 <- rep(4,length(regRegion))
#' 
#' # Overlapping quantified enhancer activity and gene expression:
#' # OUTPUT  regActivityAroundTSS():
#' regActivityAroundTSS(regActivity=regRegion, tss=geneExpression,
#' upstream=1,downstream=1)
#' 
#' # when different upstream/downstream argument is used:
#' 
#' regActivityAroundTSS(regRegion,geneExpression,upstream=5,downstream=5)
#' 
#' # force=T, forcing 1-to-1 relationship
#'  regActivityAroundTSS(regRegion[1:length(geneExpression)],
#'  geneExpression,upstream=5,downstream=5,force=TRUE)
#' 
#' @details only enhancers located within (+/-)upstream/downstream of TSS 
#' are identified,extracted and reported in output (together with info
#' about gene expression). Sample id's (corresponding to the cell types 
#' or conditions) are included only if both, 1) gene expression values and
#' 2) quantified regulatory activity are available in TSS and 
#' regActivity objects. Non-overlapping cell types are excluded.
#'  
#' 
#' 
#' @export
regActivityAroundTSS <- function(regActivity,
                                 tss,
                                 upstream=500000,
                                 downstream=500000,
                                 mc.cores=1,
                                 force=FALSE){
  
  
  if (is.null(tss$name)) {stop("TSS object does not contain info about gene 
                               name. TSS GRanges object should be 2nd arg")}
  
  
  # getting modelling data for testing 1-to-1 relationship
  
    if (force==TRUE){
    
    if (length(regActivity)!=length(tss)){stop(
          "Force=T, but data doesn't have an equal length")
      
      }else if (length(regActivity)==length(tss)){
    
        tssActivity <- parallel::mclapply(1:length(tss),function(x){
    
      GeneInfo <- tss[x]
      EnhancerInfo <- regActivity[x]
    
      #keep names, because GeneInfo is stripped off this info afterwards
      names <- c(GeneInfo$name,GeneInfo$name2)
      
      # getting shared info between gene expression and enh activity datasets
      CommonData <- colnames(mcols(GeneInfo))[colnames(mcols(GeneInfo))%in%
                                                colnames(mcols(EnhancerInfo))]
      
      # adjusting for common cell types
          mcols(GeneInfo) <- mcols(GeneInfo)[CommonData]
          mcols(EnhancerInfo) <- mcols(EnhancerInfo)[CommonData]
          
      EnhGene <- c(GeneInfo,EnhancerInfo)
            # adding meta-data   
      EnhGene$featureType <- c("gene","regulatory")
      EnhGene$name <- c(names[1],as.character(EnhancerInfo))
      EnhGene$name2 <- c(names[2],as.character(EnhancerInfo))
      
      #ordering meta-data
      FirstMetaData <- c("featureType","name","name2")
      BigWigMetaData <- colnames(mcols(EnhGene))[!colnames(mcols(EnhGene))
                                                 %in%FirstMetaData]
        mcols(EnhGene) <- mcols(EnhGene)[c(FirstMetaData,BigWigMetaData)]
      
          return(EnhGene)
            })
  
         names(tssActivity) <- tss$name
  
      }
  }
  
  # testing many-to-many relationship
  
    if (force==FALSE){
      
      # per gene analysis
      tss.extended <- split(tss,
                            as.character(tss$name))
    
        tssActivity <- parallel::mclapply(tss.extended,function(x){
          
          # overlap extended TSS (+/-downstream and upstream) & regRegion
        tss.regAct.overlap <- as.data.frame(findOverlaps(
                trim(suppressWarnings(promoters(x,upstream,downstream))),
                                                         regActivity))
          
          # if there is any overlap - GO!
          if (nrow(tss.regAct.overlap)!=0) {
            regActivity <- regActivity[tss.regAct.overlap$subjectHits]
            
            # adjusting mcols for reg region
                  name <- as.character(regActivity)
                  featureType <- rep("regulatory",length(name))
                  name2 <- rep(NA,length(name))
                  
                  mcols(regActivity) <- cbind(featureType,
                                              name,
                                              name2,
                                              data.frame(mcols(regActivity)))
                  
            # adjusting mcols for gene expression object
                  featureType <- "gene"
                  mcols(x) <- cbind(featureType,
                                    data.frame(mcols(x)))
                  
                  tss.colnames <- colnames(mcols(x))[colnames(mcols(x))%in%
                                                     colnames(mcols(regActivity))]
                  mcols(x) <-  mcols(x)[tss.colnames]
                  mcols(regActivity) <-  mcols(regActivity)[tss.colnames]
                  
                  geneReg <- c(x,regActivity)
                  geneReg$featureType  <- as.character(geneReg$featureType)
                  geneReg$name   <- as.character(geneReg$name)
                  geneReg$name2   <- as.character(geneReg$name2)
                  
            return(geneReg)
            
          }
        },
        mc.cores = mc.cores)
  
    }
  
  
  # removing NULL entries
  tssActivity <- tssActivity[!sapply(tssActivity,is.null)]
  # returning GRangesList
  tssActivity <- GRangesList(tssActivity)
  
  return(tssActivity)
  
}





#' Quantifies gene expression measured by RNA-Seq  
#' 
#' 
#' The function quantifies gene expression as a sum of exon 
#' expressions quantified over pre-defined exon regions using 
#' signal from RNA-Seq tracks (bigwig files) and returns a
#' GRanges object with TSS location and corresponding gene
#' expression levels quantified over a set of samples. 
#'  
#' 
#' 
#' @param exons A \code{\link[InteractionSet]{GInteractions}} object with 
#' Anchor1 representing the exon locations, and Anchor2 represents location of 
#' the corresponding gene. A meta-columns with following info should be present:
#' 1) name (character)- ENSEMBL or other gene identifier; 2) name2 (character,
#' optional) - 2nd gene identifier. 
#' Optionally, GRanges object can be used as input object as well. It should
#' contain exon regions over which the expression will be calculated.
#' For this object, a meta-columns with following info are necessary: 
#' 1) reg - GRanges object - gene location, or TSS location, 2) name (character)
#' - ENSEMBL or other gene identifier; 3) name2(character,optional) - 
#' 2nd gene identifier.  
#' It is strongly suggested to adjust seqlengths of this object to be equal to 
#' the seqlenghts Hsapiens from the appropriate package 
#' (BSgenome.Hsapiens.UCSC.hg19 or whatever needed version).
#' @param geneActSignals a named list of RNA-Seq BigWig files. Names correspond 
#' to the unique sample ids/names. Stranded and unstranded libraries allowed.
#'  BUT!!! It is crucial that forward and reverse RNA-Seq libraries are listed
#' in a row (eg one on top of each other)
#' @param sampleIDs NULL (default). A vector of unique sample
#' ids/names(.bw files), ordered as the bigwig files are ordered. When NULL
#' basenames of .bw files is used as a unique sample ids/names.
#' @param libStrand a vector of "*","+,"-" (default NULL) which needs to be 
#' entered as an argument.This vector provides info about the order of RNA-Seq 
#' libraries based on their strandness: "+" corresponds to forward/positive 
#' RNA-Seq bigwig files; "-" corresponds to reverse/negative RNA-Seq bigwig
#' files and "*" is unstranded library. When all libraries are unstranded then
#' the vector should contain a list of "*" with the lenght equal to the number
#' of analyzed RNA-Seq libraries (eg bigwig files). 
#' If libStrand=NULL than function will do it automatically, eg
#' create a vector of "*".It is crucial that stranded RNA-Seq libraries are
#' listed in a row (eg one on top of each other)
#' @param summaryOperation "mean"(default). An argument for 
#' \code{\link[genomation]{ScoreMatrixBin}} that is in the nutshell of 
#' quantifying exon expression across pre-defined exon regions.
#' This designates which summary operation should be used over the regions.
#' Currently, only mean is available, but "median" or "sum" will be implemented
#' in the future.
#' @param mc.cores (def:1) Define the number of cores to use;
#' at most how many child processes will be run simultaneously 
#' using mclapply from parallel package. Parallelization requires at 
#' least two cores.
#' @param normalize NULL(default). Optional "quantile" and "ratio"
#' If set to "quantile" activity measures are quantile normalized as 
#' implemented in \code{\link[preprocessCore]{normalize.quantiles}} and
#' returned ; if set to "ratio" then "median ratio method" implemented
#' as \code{\link[DESeq2]{estimateSizeFactorsForMatrix}} is used to normalize 
#' activity measures.
#' 
#'         
#' @return a GRanges object where its meta-columns correspond to quantified
#' gene expression across cell type or condition. GRanges correspond
#' to the TSS location of the analyzed gene. Additionally, each TSS 
#' contains following metadata:a "name" and "name2" columns for 
#' unique id or name/symbol for the gene which the TSS is associated with. 
#' One is Ensembl id and the other could be used for the gene symbol. 
#' Other metadata column names should represent sample names/ids and should
#' match the GRanges object provided via regActivity argument.
#' 
#' @details For a gene of interest and across all samples (.bw files)
#' individaully,the expression is firstly calculated for all exon 
#' regions which correspond to the gene of interest. Then, per gene
#' exon expression scores are summed together and divided by a total exon length
#' as followed:
# #'  \deqn{Gene Expression=\sum{n=1}^{Exons{n}}
# #'  (MeanExonExpression{n}*ExonLength{n})/(GeneLength)}
#' Normalization is runned on the level of gene expression. Currently, the 
#' relevant bigWig files are required to calculate activity. This function might 
#' be extended to work with BAM files in the future. RNA-Seq .bw files can 
#' originate from stranded,unstranded or mixed libraries.
#' 
#' 
#' 
#' @examples #INPUT1 DEFINING .BW FILES:
#' 
#' 
#' test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
#' test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
#' 
#' #INPUT2 DEFINING EXONS:
#' exons <- GRanges(c(rep("chr1",2),"chr2",rep("chr1",3)),
#'                       IRanges(c(1,7,9,15,1,21),c(4,8,14,20,4,25)),
#'                                             c(rep("+",3),rep("-",3)))
#'  exons$reg <-  exons[c(1,1,3,5,5,5)]
#'  exons$name2 <- exons$name <- paste0("TEST_Reg",c(1,1,3,5,5,5))
#'  bwToGeneExp(exons = exons,geneActSignals = c(test.bw,test2.bw),
#'          sampleIDs=c("CellType1","CellType2"))
#'  
#'  # OUTPUT bwToGeneExp():                                                                                                                       
#' bwToGeneExp(exons = exons,geneActSignals = c(test.bw,test2.bw))
#' 
#' # adding different sample IDs
#' 
#' bwToGeneExp(exons = exons,geneActSignals = c(test.bw,test2.bw),
#'             sampleIDs=c("CellType1","CellType2"))
#' 
#' # if exons input is GInteractions class object,the same output is obtained
#' 
#' 
#' exonsGI= GInteractions(exons,exons$reg)
#'    exonsGI$name=exons$name
#'    exonsGI$name2=exons$name2
#'    
#' bwToGeneExp(exons = exonsGI,geneActSignals = c(test.bw,test2.bw))
#' 
#' @export
bwToGeneExp <- function(exons,
                        geneActSignals,
                        sampleIDs=NULL,
                        libStrand=NULL,
                        summaryOperation="mean",
                        mc.cores=1,
                        normalize=NULL){
  
  
 
  #test if there is any info about genes in exon granges object
  if (any(stringr::str_detect(colnames(mcols(exons)),"name"))==0){
    stop("No info about the gene name, e.g. name column is missing")}
  
  if ((class(exons)=="GRanges")&
      (any(stringr::str_detect(colnames(mcols(exons)),"reg"))==0)){
    stop("No info about the gene location, e.g. reg column is missing ")}
  
  if ((class(exons)=="GInteractions")&
      (any(stringr::str_detect(colnames(mcols(exons)),"name2"))==0)){
    stop("No info about gene name2,e.g. name2 column needed for GInteractions")}
  
  
  
  if (is.null(sampleIDs)) {
    bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
    sampleIDs <- stringr::str_replace(basename(geneActSignals),
                                      paste(bw.exts,collapse="|"),"")
  }
  
  # adjusting for GInteractions object input - rearranging to GRanges object
  
  if (class(exons)=="GInteractions"){
    
    require(GenomicInteractions)
    
    exons <- GRanges(anchorOne(exons),
                        reg=anchorTwo(exons),
                        name=exons$name,
                        name2=exons$name2)
    
  }
  
  
  # order exons such that you wouldn't have problems afterwards with orientation
  exons <- exons[order(as.character(strand(exons)))]
  
  # separate genes form forward and reverse strand
  exons.splitted <- split(exons,strand(exons))
  gene.metadata <- length(mcols(exons)) # metadata info about corresponding gene
  
  # arranging a problem of strandness of RNA-Seq libraryies
  # 1. if object libStrand is empty then strands are assumed to be "*"
  if (is.null(libStrand)) {libStrand <- rep("*",length(geneActSignals))}
  
  # quantify expression per exon with adjustment to strandness of RNA-Seq 
  # libraries
  ExonExpression <- exonExpressionStrandAdjusted(exonsForQuant = exons,
                                                 exonsSplitted = exons.splitted,
                                                 geneActSignals = geneActSignals,
                                                 sampleIDs=sampleIDs, 
                                                 libStrand=libStrand,
                                                 mc.cores=mc.cores,
                                                 summaryOperation=summaryOperation)         
  # quantify expression per gene   
  exonExpressionPerGene <- split(ExonExpression,
                                 as.character(ExonExpression$name))
  
  GeneExpression <- list()
     for (i in 1:length(exonExpressionPerGene)){
           gene.id <- names(exonExpressionPerGene)[[i]]
          GeneExpression[[gene.id]] <- quantifyGeneExpression(
                                              exonExpressionPerGene[[gene.id]],
                                              gene.metadata)}
  
  
    GeneExpression.df <- do.call("rbind",GeneExpression)

   #############            
   # normalization
   if (!is.null(normalize)){
     GeneExpression.df <-  normalizeScores(GeneExpression.df,
                                           normalize)}


        # arranging gene coordinates as TSS
    Genes.object <- data.frame(unique(mcols(ExonExpression)[1:gene.metadata]),
                               stringsAsFactors = F)
    
    
        #extract gene location
              # TSS coord of genes
    
                 tss <- promoters(GRanges(Genes.object$reg.seqnames,
                                IRanges(Genes.object$reg.start,
                                        Genes.object$reg.end),
                                strand=Genes.object$reg.strand),
                                1,1)


        # adding meta-data
             name <- Genes.object$name
             name2 <- Genes.object$name2
                       if (length(name2)==0){
                         name2 <-  rep(NA,length(name))}# if no 2nd name
             metadata <- data.frame(name,
                                    name2,
                  GeneExpression.df[match(name,rownames(GeneExpression.df)),],
                            stringsAsFactors = F)
    
             if (length(colnames(metadata))==3) {
               
               colnames(metadata)[3] <- sampleIDs[1]
               
               }
             
          mcols(tss) <- metadata

  return(tss)
  
} 







