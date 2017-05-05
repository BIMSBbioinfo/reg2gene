.ScoresAsMcols <- function(scores.per.cell.type,activitySignals){
  
  # function that rearranges list of scores calculated per cell type
  mcols.per.cell.type <- rbind.data.frame(scores.per.cell.type,stringsAsFactors =F, make.row.names=T)
  
  # adjust column names - remove bw extension
  bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
  colnames(mcols.per.cell.type) <- str_replace(basename(activitySignals),paste(bw.exts,collapse="|"),"")
  
  return(mcols.per.cell.type)
}


.NormalizeScores <- function(ScoresDF,NormalizationMet=normalize){
  
    require(DESeq2)
    require(preprocessCore)
  
          if (NormalizationMet=="ratio") {   
                # estimate size factors and multiply gene expression using these factors
                      sizeFactors <- estimateSizeFactorsForMatrix(as.matrix(ScoresDF))
                      ScoresDF.normalized <- as.matrix(ScoresDF)*rep(sizeFactors,each=nrow(ScoresDF))}
                
          if (NormalizationMet=="quantile") {  
            ScoresDF.normalized <- normalize.quantiles(as.matrix(ScoresDF))
                    rownames(ScoresDF.normalized) <- rownames(ScoresDF)
                    colnames(ScoresDF.normalized) <- colnames(ScoresDF)}
           
 
  return(ScoresDF.normalized)

  
  }
  


#' Calculates regulatory activity over pre-defined regions
#' 
#' The function calculates regulatory activity from histone
#' modification, DNAse or methylation signals for pre-defined regulatory
#' regions and returns a GRanges object with regulatory region locations
#' and their activity over a set of samples.
#' 
#' @param regRegions a GRanges object that contains regulatory regions
#' over which the regulatory activity will be calculated.
#' @param activitySignals a named list of BigWig files. Names correspond to 
#'        unique sample ids/names.
#' @param isCovNA (def:FALSE), if this is set to TRUE, uncovered
#' bases are set to NA, this is important when dealing with methylation
#' data, where uncovered bases on a bigWig file
#'  do not mean zero methylation.
#' @param summaryOperation "mean"(default),"median" or "sum". This
#' designates which summary operation should be used over the regions
#' 
#' @param normalize NULL(default). Optional "quantile" and "ration"
#' If set to "quantile" returned activity measures are quantile normalized. 
#' If set to "ration" then "median ration method" implemented as 
#' DESeq2::estimateSizeFactorsForMatrix() is used
#' 
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
#'      regActivity <- regActivity(regRegions,activitySignals)
#'      regActivity
#' 

regActivity<-function(regRegions,activitySignals,
                      isCovNA=FALSE,normalize=NULL){
  
  
          # test input - ranges
          if (!exists("regRegions")) {stop("regRegions object missing")}
          # test input - bigwig files
          if (!min(sapply(activitySignals,file.exists))) {stop("activitySignals object missing") }
          
  # adjust chromosome sizes 
  seqlengths(regRegions) <- seqlengths(Hsapiens)[names(seqlengths(regRegions))]        
  
          
          
  # calculating coverage
          scores.per.cell.type <- lapply(activitySignals,ScoreMatrixBin,windows = regRegions,bin.num = 1,type = "bigWig", is.noCovNA=isCovNA,bin.op=summaryOperation)
          #scores.exp=mclapply(activitySignals,ScoreMatrixBin,windows = regRegions,bin.num = 1,type = "bigWig", is.noCovNA=isCovNA),mc.cores=round(length(activitySignals)/5))
          
  # add normalization
          if (!is.null(normalize)) {scores.per.cell.type <-  .NormalizeScores(ScoresDF=scores.per.cell.type,NormalizationMet=normalize)}
          
           # adding scores as mcols  
          mcols(regRegions) <- .ScoresAsMcols(scores.per.cell.type,activitySignals)
         
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
#' @param TSS a GRanges object that contains the TSS location and associated
#' gene expression values per cell type or condition as meta data. Each
#' row should have a "name" and "name2" columns for unique id or name/symbol
#' for the gene which the TSS is associated with. One could be Ensembl id and the
#' other could be used for gene symbol.
#' Other metadata column names should represent sample names/ids and should
#' match the GRanges object provided via regActivity argument.
#' @param upstream number of basepairs upstream from TSS to look for regulatory
#' regions. default 500kb
#' @param downstream number of basepairs downstream from TSS to look for regulatory
#' regions. default 500kb
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


regActivityAroundTSS <- function(regActivity,TSS,upstream=500000,
                                 downstream=500000){
  
  
  # extend TSS for wanted region
  TSS.extended <- promoters(TSS,upstream=upstream,downstream=downstream)
  # getting EnhRegions which overlap extended TSS
        TSS.regAct.overlap <- as.data.frame(findOverlaps(TSS.extended,regActivity))
        regActivity <- regActivity[TSS.regAct.overlap$subjectHits]
        
  # adjusting mcols
        name <- as.character(regActivity)
        featureType <- rep("regulatory",length(name))
        name2 <- rep(NA,length(name))
        gene.indicator <- TSS.regAct.overlap$queryHits
  
    mcols(regActivity) <- cbind(gene.indicator,featureType,name,name2,data.frame(mcols(regActivity)))
     
  # adjusting mcols for gene expression object
        featureType <- "gene"
       # adding gene indicator
            gene.indicator <- unique(TSS.regAct.overlap$queryHits)
            mcols(TSS) <- cbind(gene.indicator,featureType,data.frame(mcols(TSS)))

  
  # adjusting corresponding colnames
  
          TSS.colnames <- colnames(mcols(TSS))[colnames(mcols(TSS))%in%colnames(mcols(regActivity))]
                   mcols(TSS) <-  mcols(TSS)[TSS.colnames]
                   mcols(regActivity) <-  mcols(regActivity)[TSS.colnames]
          
 
  TSSexprRegActivity <- c(TSS,regActivity)
      TSSexprRegActivity$gene.indicator <- TSS$name[TSSexprRegActivity$gene.indicator]
      TSSexprRegActivity <- split(TSSexprRegActivity,TSSexprRegActivity$gene.indicator)
     
  
  
  return(TSSexprRegActivity)
  
}




#' Quantifies gene expression measured by RNA-Seq  
#' 
#' 
#'  The function quantifies gene expression over pre-defined gene
#'  regions using signal from RNA-Seq tracks (bigwig files) and 
#'  returns a GRanges object with TSS location and corresponding gene
#'  expression levels quantified over a set of samples. 
#'  For a given geneExpression is firstly calculated across all exon 
#'  regions. Then, per gene exon expression  scores are summed together
#'  and divided by a total exon length. 
#' 
#' 
#'  @param Exons A GRanges object that contains exon regions over which 
#'  the expression will be calculated. A meta-column with corresponding 
#'  gene names is necessary.
#'    
#'  @param GeneExpSignals a named list of RNA-Seq BigWig files. Names correspond to 
#'       the unique sample ids/names. Stranded and unstranded libraries allowed.
#'        
#'  @param LibStrand a vector of "*","+,"-" (default NULL which is in scripts
#'  tackled as "*" strand for all entered RNA-Seq .bw tracks. Otherwise a vector
#'  of "*","+,"-" is expected
#'   to the + and - strand) are expected to be lsted in a row 
#' bases are set to NA, this is important when dealing with methylation
#' data, where uncovered bases on a bigWig file
#'  do not mean zero methylation. 
#' 
#' @param is.noCovNA Possibly, need to test
#' 
#' @param mc.cores (def:1)
#' 
#' @param normalize NULL(default). Optional "quantile" and "ration"
#' If set to "quantile" returned activity measures are quantile normalized. 
#' If set to "ration" then "median ration method" implemented as 
#' DESeq2::estimateSizeFactorsForMatrix() is used
#' 
#' 
#' @return
#' 
#' @details  regulatory activity is measured by averaging logFC for
#' histone modification ChIP-seq profiles, or DNAse signal, or methylation
#' per base.Currently, relevant bigWig files are required to calculate activity       
#' activity. This function might be extended to work with BAM files
#' in the future. 
#' 
#' @import GenomicRanges
#' @import genomation
#' @import parallel
#' @import DESeq2
#' @import stringr
#' @import preprocessCore
#' 
#' 
#' 
#' @examples
#' 
#' load(file="~/GeneExpSignals.RData")
#' load(file="~/LibStrand.RData")
#' load(file="~/Exons.RData")
#' 
#' 
#' bwToGeneExp(Exons=Exons,GeneExpSignals=GeneExpSignals,LibStrand=LibStrand, mc.cores=1,normalize="ratio")
#' bwToGeneExp(Exons=Exons,GeneExpSignals=GeneExpSignals,LibStrand=LibStrand, mc.cores=1,normalize="quantile")
#' bwToGeneExp(Exons=Exons,GeneExpSignals=GeneExpSignals,LibStrand=LibStrand, mc.cores=1,normalize=NULL)
#' 
#' 

.ExonExpressionStrandAdjusted <- function(Exons,GeneExpSignals,LibStrand=NULL,
                                          mc.cores=mc.cores) {
  
  # function that separates stranded and unstranded libraries and quantifies exon expression separately for stranded libraries
  # eg expression of exons located on the + strand is quantified using forward libraries 
  # if libraries are unstranded then exon orientations is not taken into account
  
  # detect stranded tracks and unstranded tracks    
  ForwardLibraries <- str_detect(LibStrand,"\\+") 
  Unstranded.libraries <- str_detect(LibStrand,"\\*")
  
  # separate analysis for unstranded & stranded libraries
  
  if (sum(Unstranded.libraries)!=length(LibStrand)){
    # test strandness - "mixed"&"stranded"
    # Stranded libraries - positive/forward strand
    
    # forward strand
    ExonScoresForwardStrand <- mclapply(GeneExpSignals[ForwardLibraries], 
                                        function(x) {try(ScoreMatrixBin(x,windows =  Exons$`+`,bin.num = 1,type = "bigWig",bin.op="mean", is.noCovNA=T),
                                                         silent = T)},mc.cores=mc.cores)
    ExonScoresForwardStrand.df <- .ScoresAsMcols(ExonScoresForwardStrand,GeneExpSignals[ForwardLibraries])
    
    # reverse strand               
    ExonScoresReverseStrand <- mclapply(GeneExpSignals[!(ForwardLibraries|Unstranded.libraries)], 
                                        function(x) {try(ScoreMatrixBin(x,windows =  Exons$`-`,bin.num = 1,type = "bigWig", bin.op="mean",is.noCovNA=T),
                                                         silent = T)},mc.cores=mc.cores)
    ExonScoresReverseStrand.df <- abs(.ScoresAsMcols(ExonScoresReverseStrand,GeneExpSignals[ForwardLibraries]))
    # pool reverse and forward strand together        
    Exons <- unlist(Exons)
    values(Exons) <- (cbind(mcols(Exons),rbind(ExonScoresForwardStrand.df,ExonScoresReverseStrand.df)))
    
  }
  
  
  
  if (sum(Unstranded.libraries)!=0){
    
    Exons <- unlist(Exons)
    # Unstranded
    ExonScoresUnstranded <- mclapply(GeneExpSignals[Unstranded.libraries], 
                                     function(x) { try(ScoreMatrixBin(x,windows = Exons,bin.num = 1,type = "bigWig",bin.op=bin.op,is.noCovNA=T),
                                                       silent = T)},mc.cores=mc.cores)
    ExonScoresUnstranded.df <- .ScoresAsMcols(ExonScoresUnstranded,GeneExpSignals[Unstranded.libraries])                  
    
    values(Exons) <- (cbind(mcols(Exons),rbind(ExonScoresUnstranded.df)))
    
    
  }
  
  return(Exons)                
}


.QuantifyGeneExpression <- function(ExonExpressionPerGene){
  
  # sum across exons (mean value across exon*exon width)/ gene_length(sum_of_all_exon_lengths)
  
  ExonLength <- width(ExonExpressionPerGene)
  GeneLength <- sum(ExonLength,na.rm=T)
  ExpressionPerSample=as.data.frame(mcols(ExonExpressionPerGene)[-(1:7)])
  
  GeneExpression <- apply(ExpressionPerSample*ExonLength,2,sum,na.rm=T)/GeneLength
  
  return(GeneExpression)
}


.ScoresAsMcols <- function(scores.per.cell.type,activitySignals){
  
  # function that rearranges list of scores calculated per cell type
  mcols.per.cell.type <- rbind.data.frame(scores.per.cell.type,stringsAsFactors =F, make.row.names=T)
  
  # adjust column names - remove bw extension
  bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
  colnames(mcols.per.cell.type) <- str_replace(basename(activitySignals),paste(bw.exts,collapse="|"),"")
  
  return(mcols.per.cell.type)
}


bwToGeneExp <- function(Exons,GeneExpSignals,LibStrand=NULL,bin.op="mean",
                        mc.cores=1,normalize="ratio"){
  
  # adjust chromosome sizes 
  seqlengths(Exons) <- seqlengths(Hsapiens)[names(seqlengths(Exons))]
  
  # separate genes form forward and reverse strand
   Exons.splitted <- split(Exons,strand(Exons))


        # arranging a problem of strandness of RNA-Seq libraryies
        # 1. if object LibStrand is empty then strands are assumed to be "*"
             if (is.null(LibStrand)) {LibStrand <- rep("*",length(GeneExpSignals))}
        
                    # quantify expression per exon with adjustment to strandness of RNA-Seq libraries
                       ExonExpression <- .ExonExpressionStrandAdjusted(Exons=Exons.splitted,GeneExpSignals=GeneExpSignals,
                                                                       LibStrand=LibStrand,mc.cores=mc.cores)
                       
                    # quantify expression per gene   
                       ExonExpressionPerGene <- split(ExonExpression,as.character(ExonExpression$sample))
                            GeneExpression <- lapply(ExonExpressionPerGene,.QuantifyGeneExpression)
                            GeneExpression.df <- do.call("rbind",GeneExpression)
                     
                  
                  # normalization
                       if (!is.null(normalize)) {GeneExpression.df <-  .NormalizeScores(ScoresDF=GeneExpression.df,NormalizationMet=normalize)}
          
       
          # arranging gene coordinates as TSS  
               Genes.object <- as.data.frame(mcols(ExonExpression)[1:5]) # extract gene location
                # TSS coord of genes
                   Genes.GRanges <- promoters(reduce(GRanges(Genes.object$seqnames,IRanges(Genes.object$start,Genes.object$end),strand=Genes.object$strand)),1,1)
                  
                      
          # adding meta-data    
               Gene.names <- unique(as.character(ExonExpression$sample))
               mcols(Genes.GRanges) <- cbind(Gene.names,GeneExpression.df[match(Gene.names,rownames(GeneExpression.df)),])
            
                            
                return(Genes.GRanges) 
         
          
} 
     





library(genomation)
library(stringr)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg19)

load(file="~/GeneExpSignals.RData")
load(file="~/LibStrand.RData")
load(file="~/Exons.RData")
mc.cores=1
bin.op="mean"
normalize="quantile"


bwToGeneExp(Exons=Exons,GeneExpSignals=GeneExpSignals,LibStrand=LibStrand, mc.cores=1,normalize="ratio")
bwToGeneExp(Exons=Exons,GeneExpSignals=GeneExpSignals,LibStrand=LibStrand, mc.cores=1,normalize="quantile")
bwToGeneExp(Exons=Exons,GeneExpSignals=GeneExpSignals,LibStrand=LibStrand, mc.cores=1,normalize=NULL)

