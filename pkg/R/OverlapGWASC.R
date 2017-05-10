
#' @param GWASCat - GRanges with of SNP positions,
#'  mcols correspond to rsSNP, LDrsSNP ,associated gene, and disease

#' @param RegCat - GRanges with of SNP positions,
#'  mcols correspond to rsSNP, LDrsSNP ,associated gene, and disease
#'  
#'  
#'  @param orientation - character "GWASCat" or "RegCat" (default GWASCat). 
#'  Whether to report GWASCatalog entries taht overlap RegCatalog, or other
#'  Regulatory regions which overlap GWAS SNPs (eg RegCat). 
#'  
#'  @import GenomicRanges



GWASCat <- readRDS("D:/Projects_Helping/PhD/GWASCat_example_17_05_10.rds")
RegCat <- readRDS("D:/Projects_Helping/PhD/RegCat_17_05_10.rds")



OverlapGWASC <- function(GWASCat,RegCat){
  
  
  require(GenomicRanges)
  # need to add part where input df will be organized as dataframes
  
  
  # Adjustments for GWAS 
  GWASCat.gr <- GRanges(GWASCat[,"chr1"],IRanges(as.integer(GWASCat[,"start1"]),
                                                 as.integer(GWASCat[,"end1"])))
  colnames(GWASCat) <- paste0("GWAS_'",colnames(GWASCat))
  mcols(GWASCat.gr) <- GWASCat
  # Adjustments for RegCatalog 
  #ensuring that all columns are equal      
  RegCat <- apply(RegCat,2,as.character)
  
  # creating GR  
  RegCat.gr <- GRanges(RegCat[,"chr1"],IRanges(as.integer(RegCat[,"start1"]),
                                               as.integer(RegCat[,"end1"])))
  colnames(RegCat) <- paste0("RegCat_'",colnames(RegCat))    
  mcols(RegCat.gr) <- RegCat
  
  GWASCat.RegCat.Overlap <- data.frame(findOverlaps(GWASCat.gr,RegCat.gr))
  
  GWASCatRegCat <- cbind(GWASCat[GWASCat.RegCat.Overlap$queryHits,],
                         RegCat[GWASCat.RegCat.Overlap$subjectHits,])
  
  return(GWASCatRegCat)
  
  
}

OverlapGWASC(GWASCat,RegCat)
