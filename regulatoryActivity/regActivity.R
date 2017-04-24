.ScoresAsMcols <- function(scores.per.cell.type,activitySignals){
  
  # function that rearranges list of scores calculated per cell type
  mcols.per.cell.type <- rbind.data.frame(scores.per.cell.type,stringsAsFactors =F, make.row.names=T)
  
  # adjust column names - remove bw extension
  bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
  colnames(mcols.per.cell.type) <- str_replace(basename(activitySignals),paste(bw.exts,collapse="|"),"")
  
  return(mcols.per.cell.type)
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
#' @param normalize NULL(default). If set to "quantile" returned activity
#' measures are quantile normalized
#' 
#' @return a GRanges object where its meta-columns correspond
#'         to calculated acvitity measures and column names 
#'         correspond to provided sample ids or names.
#' 
#' @import genomation
#' @import GenomicRanges
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
                      isCovNA=FALSE,summaryOperation="mean",
                      normalize=NULL){
  
  
          # test input - ranges
          if ( !exists("regRegions")) { stop("regRegions object missing") }
          # test input - bigwig files
          if ( !min(sapply(activitySignals,file.exists))) { stop("regRegions object missing") }
          
          
          
          #scores.exp=mclapply(activitySignals,ScoreMatrixBin,windows = regRegions,bin.num = 1,type = "bigWig", is.noCovNA=isCovNA),mc.cores=round(length(activitySignals)/5))
          # calculating coverage
          scores.per.cell.type <- lapply(activitySignals,ScoreMatrixBin,windows = regRegions,bin.num = 1,type = "bigWig", is.noCovNA=isCovNA,bin.op=summaryOperation)
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
#' @param TSS a GRanges object that have the TSS location and associated
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
#' are unique gene ids/names. Metadata for a GRanges
#' object in the list represents regulatory activity and gene expression accross
#' the same samples. The GRanges objects have the following 
#' metadata columns:
#' featureType: either "gene" or "regulatory"
#' name: name/id for gene and enhancers. Gene name could be id from a database
#' enhancer name should be in the format as follows "chr:start-end"
#' name2: a secondary name for the feature, such as gene symbol "PAX6" etc. not
#' necessary for enhancers could be NA
#' other columns: numeric values for gene expression or regulatory actvity.
#' Column names represent sample names/ids.
#' 
#' @examples #NB provide minimal examples that work on beast
#' 
#' @details 
#' 
#' @import # provide a list of required packages and functions
#' 
regActivityAroundTSS<-function(regActivity,TSS,upstream=500000,
                               downstream=500000){
  
}
