.ScoresAsMcols <- function(scores.per.cell.type,activitySignals,sampleIDs){
  
  require(stringr)
  # function that rearranges list of scores calculated per cell type
  
        mcols.per.cell.type <- rbind.data.frame(scores.per.cell.type,stringsAsFactors=F,
                                                make.row.names=T)
        
        # adjust column names - remove bw extension
        
        
        colnames(mcols.per.cell.type) <- sampleIDs
       
        
        return(mcols.per.cell.type)
}


.NormalizeScores <- function(ScoresDF,NormalizationMet=normalize){
  # f() that runs normalization
  
    require(DESeq2)
    require(preprocessCore)
  
          if (NormalizationMet=="ratio") {   
                # estimate size factors and multiply gene expression using these factors
                      sizeFactors <- estimateSizeFactorsForMatrix(as.matrix(ScoresDF))
                      ScoresDF.normalized <- as.matrix(ScoresDF)*rep(sizeFactors,
                                                                     each=nrow(ScoresDF))}
                
          if (NormalizationMet=="quantile") {  
            ScoresDF.normalized <- normalize.quantiles(as.matrix(ScoresDF))
                    rownames(ScoresDF.normalized) <- rownames(ScoresDF)
                    colnames(ScoresDF.normalized) <- colnames(ScoresDF)}
           
 
  return(ScoresDF.normalized)

  
  }
  

.ExonExpressionStrandAdjusted <- function(Exons,GeneExpSignals,sampleIDs,LibStrand,
                                          mc.cores=1,summaryOperation="mean") {
  
  # function that separates stranded and unstranded libraries and quantifies 
  # exon expression separately for stranded libraries
  # eg expression of exons located on the + strand is quantified using forward libraries 
  # if libraries are unstranded then exon orientations is not taken into account
        require(stringr)
  # detect stranded tracks and unstranded tracks    
  ForwardLibraries <- str_detect(LibStrand,"\\+") 
  Unstranded.libraries <- str_detect(LibStrand,"\\*")
  
  if (sum(ForwardLibraries)!=(sum(!(ForwardLibraries|Unstranded.libraries)))){
    stop("Number of forward and reverse libraries do not match")
  }
  
  
  # separate analysis for unstranded & stranded libraries
  
  if (sum(Unstranded.libraries)!=length(LibStrand)){
    # test strandness - "mixed"&"stranded"
    # Stranded libraries - positive/forward strand
    
    # forward strand
    ExonScoresForwardStrand <- mclapply(GeneExpSignals[ForwardLibraries], 
                                        function(x) {try(ScoreMatrixBin(x,
                                                                        windows =  Exons$`+`,
                                                                        bin.num = 1,
                                                                        type = "bigWig",
                                                                        bin.op=summaryOperation, 
                                                                        is.noCovNA=T),
                                                         silent = T)},mc.cores=mc.cores)
    
          if (max(sapply(ExonScoresForwardStrand,class)=="try-error")){
                stop(paste("This bigwig files do not overlap exon regions",
                            GeneExpSignals[ForwardLibraries][sapply(ExonScoresForwardStrand,
                            class)=="try-error"]))}                                      
    
    ExonScoresForwardStrand.df <- .ScoresAsMcols(ExonScoresForwardStrand,
                                                 GeneExpSignals[ForwardLibraries],
                                                 sampleIDs[ForwardLibraries])
    
    # reverse strand               
    ExonScoresReverseStrand <- mclapply(GeneExpSignals[!(ForwardLibraries|Unstranded.libraries)], 
                                        function(x) {try(ScoreMatrixBin(x,windows =  Exons$`-`,
                                                                        bin.num = 1,
                                                                        type = "bigWig", 
                                                                        bin.op=summaryOperation,
                                                                        is.noCovNA=T),
                                                         silent = T)},mc.cores=mc.cores)
    
          if (max(sapply(ExonScoresReverseStrand,class)=="try-error")){
            stop(paste("This bigwig files do not overlap exon regions",
                          GeneExpSignals[!(ForwardLibraries|Unstranded.libraries)][sapply(ExonScoresReverseStrand,
                          class)=="try-error"]))}
    
    ExonScoresReverseStrand.df <- abs(.ScoresAsMcols(ExonScoresReverseStrand,
                                                     GeneExpSignals[ForwardLibraries],
                                                     sampleIDs[ForwardLibraries]))
    # pool reverse and forward strand together        
    Exons <- unlist(Exons)
    values(Exons) <- (cbind(mcols(Exons),rbind(ExonScoresForwardStrand.df,ExonScoresReverseStrand.df)))
    
  }
  
  if (sum(Unstranded.libraries)!=0){
    
    
    # Unstranded
    ExonScoresUnstranded <- mclapply(GeneExpSignals[Unstranded.libraries], 
                                     function(x) { try(ScoreMatrixBin(x,windows = Exons,bin.num = 1,type = "bigWig",
                                                                      bin.op=summaryOperation,is.noCovNA=T),
                                                       silent = T)},mc.cores=mc.cores)
    
            if (max(sapply(ExonScoresUnstranded,class)=="try-error")){
              stop(paste("This bigwig files do not overlap exon regions",
                      GeneExpSignals[Unstranded.libraries][sapply(ExonScoresUnstranded,
                      class)=="try-error"]))}                                      
            
            
    ExonScoresUnstranded.df <- .ScoresAsMcols(ExonScoresUnstranded,
                                              GeneExpSignals[Unstranded.libraries],
                                              sampleIDs[Unstranded.libraries])                  
    
    values(Exons) <- (cbind(mcols(Exons),rbind(ExonScoresUnstranded.df)))
    
    
  }
  
  return(Exons)                
}


.QuantifyGeneExpression <- function(ExonExpressionPerGene){
  
  # f() that calculates per gene expression based on exon expression results
  # sum across exons (mean value across exon*exon width)/ 
  # gene_length(sum_of_all_exon_lengths)
  
        ExonLength <- width(ExonExpressionPerGene)
        GeneLength <- sum(ExonLength,na.rm=T)
        ExpressionPerSample=as.data.frame(mcols(ExonExpressionPerGene)[-(1:7)])
        
        GeneExpression <- apply(ExpressionPerSample*ExonLength,2,sum,na.rm=T)/GeneLength
  
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
#' @param activitySignals a named list of BigWig files. Names correspond to 
#'        unique sample ids/names.
#' @param sampleIDs NULL (default). A vector of unique sample
#' ids/names(.bw files), ordered as the bigwig files are ordered. When NULL
#' basenames of .bw files is used as a unique sample ids/names.
#' @param isCovNA (def:FALSE), if this is set to TRUE, uncovered
#' bases are set to NA, this is important when dealing with methylation
#' data, where uncovered bases on a bigWig file do not mean zero methylation.
#' @param normalize NULL(default). Optional "quantile" and "ration"
#' If set to "quantile" activity measures are quantile normalized as 
#' implemented in \code{\link[preprocessCore]{normalize.quantiles}} and
#' returned ; if set to "ration" then "median ration method" implemented
#' as \code{\link[DESeq2]{estimateSizeFactorsForMatrix}} is used to normalize 
#' activity measures.
#' @param summaryOperation "mean"(default).
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
#' @import genomation
#' @import GenomicRanges
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import preprocessCore
#' @import DESeq2
#' @import parallel
#' @import stringr
#' 
#' @details regulatory activity is measured by averaging logFC for
#' histone modification ChIP-seq profiles, or DNAse signal, or methylation
#' per base.Currently, relevant bigWig files are required to calculate activity       
#' activity. This function might be extended to work with BAM files
#' in the future. 
#' 
#' @examples 
#'      library(genomation)
#'      library(GenomicRanges)
#'      load("pkg/inst/extdata/regRegions.RData")
#'      activitySignals <- c("pkg/inst/extdata/E085-H3K27ac.chr10.fc.signal.bigwig",
#'                           "pkg/inst/extdata/E066-H3K27ac.chr10.fc.signal.bigwig")
#'      regActivity(regRegions,activitySignals)
#'      sampleIDs <- c("E085","E066")
#'      regActivity(regRegions,activitySignals,sampleIDs)
#'      
#' 



regActivity <- function(regRegions,activitySignals,sampleIDs=NULL,isCovNA=FALSE,
                      summaryOperation="mean",normalize=NULL,mc.cores=1){
  
  
          # test input - ranges
          if (!exists("regRegions")) {stop("regRegions object missing")}
          # test input - bigwig files
          if (!min(sapply(activitySignals,file.exists))) {stop("activitySignals object missing") }
          
          if (is.null(sampleIDs)) {
            
            bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
            sampleIDs <- str_replace(basename(activitySignals),
                                                         paste(bw.exts,collapse="|"),"")
          }
  
  
  
  
  # adjust chromosome sizes 
  seqlengths(regRegions) <- seqlengths(Hsapiens)[names(seqlengths(regRegions))]        
  
          
          
  # calculating coverage
         
          scores.per.cell.type <- mclapply(activitySignals,ScoreMatrixBin,
                                           windows = regRegions,
                                           bin.num = 1,type = "bigWig", 
                                           is.noCovNA=isCovNA,
                                           bin.op=summaryOperation,
                                           mc.cores = mc.cores)
          
         
           # adding scores as mcols  
          scores.per.cell.type.df <- .ScoresAsMcols(scores.per.cell.type,activitySignals,sampleIDs)
        
          
          # add normalization
          if (!is.null(normalize)) {scores.per.cell.type.df <-  .NormalizeScores(ScoresDF=scores.per.cell.type.df,
                                                                              NormalizationMet=normalize)}
          
          
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
#' @param TSS a GRanges object output from \code{bwToGeneExp} that
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
#' 
#' @return a GRangesList object per gene that contain location of TSS
#' and regulatory regions around that gene. Names for the GRangesList
#' are unique gene ids/names. 
#' Metadata for a GRanges object in the list represents regulatory 
#' activity and gene expression accross the same samples. 
#' The GRanges objects have the following metadata columns:
#'  1. featureType: either "gene" or "regulatory"
#'  2. name: name/id for gene and enhancers. Gene name could be id from a database
#'          enhancer name should be in the format as follows "chr:start-end"
#'  3. name2: a secondary name for the feature, such as gene symbol "PAX6" etc. not
#'    necessary for enhancers could be NA
#'  4. other columns: numeric values for gene expression or regulatory actvity.
#'    Column names represent sample names/ids.
#' 
#' @examples 
#'  load("~/TSS.RData")
#'  load("~/regActivity.RData")
#'  regActivityAroundTSS(regActivity=regActivity,TSS=TSS,upstream=500000,downstream=500000)
#' 
#' 
#'  load(file="~/LargeregRegions.RData")
#'  RoadmapH3K4me1 <- readRDS("/data/akalin/Projects/AAkalin_reg2gene/Results/regActivity/Roadmap_H3K4me1_meanratio.rds")
#'  RoadmapRNASeq <- readRDS("/data/akalin/Projects/AAkalin_reg2gene/Results/bwToGeneExp/Roadmap_RNASeq_meanratio.rds")
#'
#'    regActivityAroundTSS(RegRegionActivity = RoadmapH3K4me1, TSS = RoadmapRNASeq[1:10])
#' 
#' 
#' RegRegionActivity <- readRDS(pathRegActivity)
#' TSS <- readRDS(pathTSS)[1:10]
#' Method <- basename(pathRegActivity)
#' mc.cores=25
#'  regActivityAroundTSS(RegRegionActivity, 
#'                                   TSS,
#'                                   downstream=1000000,
#'                                   upstream=1000000,
#'                                   mc.cores)
#'
#' 
#' 
#' @details only enhancers located within (+/-)upstream/downstream of TSS 
#' are identified,extracted and reported in output (together with info
#' about gene expression). Sample id's (corresponding to the cell types 
#' or conditions) are included only if both, 1) gene expression values and
#' 2) quantified regulatory activity are available in TSS and 
#' regActivity objects. Non-overlapping cell types are excluded.
#'  
#' 
#' @import GenomicRanges
#' 
#' 


regActivityAroundTSS <- function(RegRegionActivity,TSS,upstream=500000,
                                 downstream=500000,mc.cores=10){
  
  
  
  # extend TSS for wanted region & split per gene such that analysis goes per gene
  TSS.extended <- promoters(TSS,upstream=upstream,downstream=downstream)
  TSS.extended <- split(TSS.extended,as.character(TSS.extended$name))
  
  # getting EnhRegions which overlap extended TSS
  
  TSSexprRegRegionActivity <- mclapply(TSS.extended,function(x){
    
    # overlap extended TSS & regRegion
    TSS.regAct.overlap <- as.data.frame(findOverlaps(x,RegRegionActivity))
    
    # if there is any overlap - GO!
    if (nrow(TSS.regAct.overlap)!=0) {
      
      RegRegionActivity <- RegRegionActivity[TSS.regAct.overlap$subjectHits]
      
      # adjusting mcols
      name <- as.character(RegRegionActivity)
      featureType <- rep("regulatory",length(name))
      name2 <- rep(NA,length(name))
      gene.indicator <- TSS.regAct.overlap$queryHits
      
      mcols(RegRegionActivity) <- cbind(gene.indicator,featureType,name,name2,
                                        data.frame(mcols(RegRegionActivity)))
      
      # adjusting mcols for gene expression object
      featureType <- "gene"
      gene.indicator <- unique(TSS.regAct.overlap$queryHits)
      
      mcols(x) <- cbind(gene.indicator,featureType,
                        data.frame(mcols(x)))
      
      TSS.colnames <- colnames(mcols(x))[colnames(mcols(x))%in%colnames(mcols(RegRegionActivity))]
      mcols(x) <-  mcols(x)[TSS.colnames]
      mcols(RegRegionActivity) <-  mcols(RegRegionActivity)[TSS.colnames]
      
      return(c(x,RegRegionActivity))
      
    }
  },mc.cores = mc.cores)
  
  
  return(TSSexprRegRegionActivity)
  
}





#' Quantifies gene expression measured by RNA-Seq  
#' 
#' 
#'  The function quantifies gene expression as a sum of exon 
#'  expressions quantified over pre-defined exon regions using 
#'  signal from RNA-Seq tracks (bigwig files) and returns a
#'  GRanges object with TSS location and corresponding gene
#'  expression levels quantified over a set of samples. 
#'  
#' 
#' 
#'  @param Exons A GRanges object that contains exon regions over which 
#'  the expression will be calculated. A meta-columns with following info are
#'  necessary: need to adjust these
#'  @param GeneExpSignals a named list of RNA-Seq BigWig files. Names correspond to 
#'       the unique sample ids/names. Stranded and unstranded libraries allowed.
#'       BUT!!! It is crucial that forward and reverse RNA-Seq libraries are listed
#'       in a row (eg one on top of each other)
#'  @param sampleIDs NULL (default). A vector of unique sample
#'  ids/names(.bw files), ordered as the bigwig files are ordered. When NULL
#'  basenames of .bw files is used as a unique sample ids/names.
#'  @param LibStrand a vector of "*","+,"-" (default NULL) which needs to be entered
#'  as an argument.This vector provides info about the order of RNA-Seq libraries
#'   based on their strandness:
#'  "+" corresponds to forward/positive RNA-Seq bigwig files, whereas 
#'  "-" corresponds to reverse/negative RNA-Seq bigwig files and "*" is unstranded
#'  library. When all libraries are unstranded then the vector should contain a list
#'  of "*" with the lenght equal to the number of analyzed RNA-Seq libraries
#'  (eg bigwig files). If LibStrand=NULL than function will do it automatically, eg
#'  create a vector of "*".It is crucial that stranded RNA-Seq libraries are
#'   listed in a row (eg one on top of each other)
#'   @param summaryOperation "mean"(default).
#'  This designates which summary operation should be used over the regions.
#'  Currently, only mean is available, but "median" or "sum" will be implemented
#'  in the future.
#'  @param mc.cores (def:1) Define the number of cores to use;
#'  at most how many child processes will be run simultaneously 
#'  using mclapply from parallel package. Parallelization requires at 
#'  least two cores.
#' @param normalize NULL(default). Optional "quantile" and "ration"
#' If set to "quantile" activity measures are quantile normalized as 
#' implemented in \code{\link[preprocessCore]{normalize.quantiles}} and
#' returned ; if set to "ration" then "median ration method" implemented
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
#' exon expression scores are summed together 
#' and divided by a total exon length. Normalization is runned
#' on the level of gene expression. Currently, the relevant bigWig 
#' files are required to calculate activity. This function might 
#' be extended to work with BAM files in the future. 
#' RNA-Seq .bw files can come from stranded,unstranded or mixed 
#' libraries.
#' 
#' @import GenomicRanges
#' @import genomation
#' @import parallel
#' @import DESeq2
#' @import stringr
#' @import preprocessCore
#' @import BSgenome.Hsapiens.UCSC.hg19
#' 
#' 
#' @examples
#' 
#' load(file="~/GeneExpSignals.RData")
#' load(file="~/LibStrand.RData")
#' load(file="~/Exons.RData")
#' load(file="~/sampleIDs.RData")
#' 
#'     bwToGeneExp(Exons,GeneExpSignals,LibStrand, sampleIDs=NULL,mc.cores=1,normalize="ratio",summaryOperation="mean")
#'     bwToGeneExp(Exons,GeneExpSignals,LibStrand, mc.cores=1,normalize="quantile", sampleIDs=NULL,summaryOperation="mean")
#'     bwToGeneExp(Exons,GeneExpSignals,LibStrand, mc.cores=1,normalize=NULL, sampleIDs=NULL,summaryOperation="mean")
#'     bwToGeneExp(Exons=Exons,GeneExpSignals=GeneExpSignals,LibStrand=LibStrand,sampleIDs=sampleIDs, 
#'                   mc.cores=1,normalize="quantile",summaryOperation="mean")
#' 




bwToGeneExp <- function(Exons,
                        GeneExpSignals,
                        sampleIDs=NULL,
                        LibStrand=NULL,
                        summaryOperation="mean",
                        mc.cores=1,
                        normalize="ratio"){
  
  # test if there is any info about genes in exon granges object
  if (sum(str_detect(colnames(mcols(Exons)),"gene.id"))==0)  {stop("No info about the gene name")}
  
  
  if (is.null(sampleIDs)) {
    
        bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
        sampleIDs <- str_replace(basename(GeneExpSignals),
                                 paste(bw.exts,collapse="|"),"")
  }
      
  
  # adjust chromosome sizes 
  seqlengths(Exons) <- seqlengths(Hsapiens)[names(seqlengths(Exons))]
  
  # separate genes form forward and reverse strand
   Exons.splitted <- split(Exons,strand(Exons))


        # arranging a problem of strandness of RNA-Seq libraryies
        # 1. if object LibStrand is empty then strands are assumed to be "*"
             if (is.null(LibStrand)) {LibStrand <- rep("*",length(GeneExpSignals))}
        
                    # quantify expression per exon with adjustment to strandness of RNA-Seq libraries
                       ExonExpression <- .ExonExpressionStrandAdjusted(Exons=Exons.splitted,
                                                                       GeneExpSignals=GeneExpSignals,
                                                                       sampleIDs=sampleIDs, 
                                                                       LibStrand=LibStrand,
                                                                       mc.cores=mc.cores,
                                                                       summaryOperation=summaryOperation)
                       
                    # quantify expression per gene   
                       ExonExpressionPerGene <- split(ExonExpression,as.character(ExonExpression$gene.id))
                            GeneExpression <- lapply(ExonExpressionPerGene,.QuantifyGeneExpression)
                            GeneExpression.df <- do.call("rbind",GeneExpression)
                     
                            
                  # normalization
                       if (!is.null(normalize)) {GeneExpression.df <-  .NormalizeScores(ScoresDF=GeneExpression.df,
                                                                                        NormalizationMet=normalize)}
          
                            
          # arranging gene coordinates as TSS  
               Genes.object <- data.frame(unique(mcols(ExonExpression)[1:7]),stringsAsFactors = F) # extract gene location
                # TSS coord of genes
                   TSS <- promoters(GRanges(Genes.object$seqnames,
                                                   IRanges(Genes.object$start,Genes.object$end),
                                                   strand=Genes.object$strand),1,1)
                  
                     
          # adding meta-data    
               name <- as.character(Genes.object$gene.id)
               name2 <- as.character(Genes.object$gene.name)
                         if (length(name2)==0) { name2 <-  rep(NA,length(name))} # if there is no 2nd name
               metadata <- data.frame(name,name2,GeneExpression.df[match(name,rownames(GeneExpression.df)),],stringsAsFactors = F)
               
            mcols(TSS) <- metadata
              
              
              
                return(TSS) 
         
          
} 
     






