#' Overlap GWAS Catalog and EnhanceRegions  
#' 
#' 
#'  The function overlaps SNPs reported in GWAS Catalog
#'  and enhancer regions and returns a dataframe of overlap
#'  regions with associated metadata. Any predfined region 
#'  in the genome can used as input.
#'  
#'       
#'  @param GWASCat - A dataframe with information about the coordinates of
#'  GWAS Catalog SNPs. Coordinates should be stored in three columns:
#'  "chr1","start1","end1". Additional information, such as
#'  SNP ID, proxy SNP ID, or reported gene can be stored in the additional 
#'  columns.
#'   
#'  @param RegCat - A dataframe with information about the coordinates of
#'  regulatory regions or any other predefined regions in the 
#'  genome. coordinates should be stored in three columns:
#'  "chr1","start1","end1". Additional information, such as
#'   gene name can be stored in the additional columns.
#'  
#'  @return A dataframe of coordinates of all overlapping regions.
#'  Each row corresponds to one-SNP-one-regulatory-region overlap 
#'  with addional metadata. Column names indicate whether coordinates and
#'  additional info come from GWASCat or RegCatalog.
#'  
#'  
#'  
#'  
#'  @details Function performs simple overlap between 2 
#'  dataframes of coordinates stored as  "chr1","start1","end1". First,
#'  it creates GRanges object, and then it runs
#'  
#'  @examples
#'    #load("~/GWASCat.RData"))
#'    #load("~/RegCat.RData"))
#'    #OverlapGWASC(GWASCat,RegCat)
#'  
#GWASCat <- readRDS("D:/Projects_Helping/PhD/GWASCat_example_17_05_10.rds")
#RegCat <- readRDS("D:/Projects_Helping/PhD/RegCat_17_05_10.rds")
OverlapGWASC <- function(GWASCat,RegCat){
  
  
  require(GenomicRanges)
  # need to add part where input df will be organized as dataframes
  
  
  # Adjustments for GWAS 
  GWASCat.gr <- GRanges(GWASCat[,"chr1"],IRanges(as.integer(GWASCat[,"start1"]),
                                                 as.integer(GWASCat[,"end1"])))
  colnames(GWASCat) <- paste0("GWAS_",colnames(GWASCat))
  mcols(GWASCat.gr) <- GWASCat
  # Adjustments for RegCatalog 
  #ensuring that all columns are equal      
    RegCat <- apply(RegCat,2,as.character)
  
  # creating GR  
    RegCat.gr <- GRanges(RegCat[,"chr1"],IRanges(as.integer(RegCat[,"start1"]),
                                               as.integer(RegCat[,"end1"])))
    colnames(RegCat) <- paste0("RegCat_",colnames(RegCat))    
    mcols(RegCat.gr) <- RegCat
  
  GWASCat.RegCat.Overlap <- data.frame(findOverlaps(GWASCat.gr,RegCat.gr))
  
  GWASCatRegCat <- cbind(GWASCat[GWASCat.RegCat.Overlap$queryHits,],
                         RegCat[GWASCat.RegCat.Overlap$subjectHits,])
  
  return(GWASCatRegCat)
  
  
}


