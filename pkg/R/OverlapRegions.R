#' Overlap regions (+/-metadata) or interacting regions in the genome 
#' 
#'  The function that compares two input regions or a combinations of two
#'  interacting regions in the genome and reportes regions present in both
#'  input datasets. If one coordinate set per input dataframe is reported
#'  ("chr1","start1","end1") then regions from one set, Reg1, is identified
#'  in the other dataset,Reg2, and repoted.
#'  If two coordinate sets ("chr1","start1","end1","chr2","start2","end2")
#'  per input dataframes are reported (eg two interacting regions in test
#'  and benchmark dataset) then both of them are compared,and if both 
#'  coordinates regions overlap they are reported.
#'  Additionaly, an overlap of 2 dataframes can be based on coordinates 
#'  and one column name as well (useful for testing overlaps of eQTL studies). 
#'   
#'       
#'  @param Reg1 - character matrix. Coordinates of the regions which we want
#'  to comape should be stored in the three columns with following names:
#'  "chr1","start1","end1". Additionally "chr2","start2","end2" columns or
#'  any other column can be a part of the matrix (such as SNP ID, proxy SNP 
#'  ID,or reported gene.
#'   
#'  @param Reg2 - character matrix. Coordinates of the regions with which
#'  Reg1 coordinates are compared should be stored in the three columns
#'  with the following names:"chr1","start1","end1". 
#'  Additionally "chr2","start2","end2" columns or any other column
#'  can be a part of the matrix (such as gene names)
#'  
#'   
#'  @param byReg1 - character (default NULL),a name of the column which
#'  will be used for additional filtering of overlapping regions (after
#'  region1 is overlapped with region2).
#'   with {name of another argument} to 
#'  fileter and report only overlapping data 
#'   
#'   @param byReg2 - character (default NULL),a name of the column which
#'  will be used for additional filtering of overlapping regions (after
#'  region1 is overlapped with region2). 
#'  
#'  with {name of another argument} to 
#'  fileter and report only overlapping data 
#'   
#'  
#'  @return A dataframe of coordinates of all overlapping regions.
#'  Each row corresponds to one-SNP-one-regulatory-region overlap 
#'  with addional metadata. Column names indicate whether coordinates and
#'  additional info come from Reg1 or Reg2alog.
#'  
#'  
#'  @import GenomicRanges
#'  @import stringr
#'  
#'  @details Function performs simple overlap between 2 
#'  dataframes of coordinates stored as  "chr1","start1","end1". First,
#'  it creates GRanges object, and then it runs
#'  \code{\link[GenomicRanges]{findOverlaps}}
#'  
#'  @example 
#'    load("~/Reg1.RData"))
#'    load("~/Reg2.RData"))
#'    OverlapGWASC(Reg1,Reg2)
#'  
#' load("~/GWASCat.RData")
#' load("~/RegCat.RData")
#' Benchmark <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Data/Benchmark1/Benchmark1_Dataset_PooledANDProcessed_5sources_Overlapping_EnhPromoterR_SHORT_2kb_17_01_26.rds")
#'
#' OverlapRegions(GWASCat,RegCat)
#' OverlapRegions(GWASCat,RegCat,byReg1="GeneName",byReg2="gene.name")
#  error function
#' OverlapRegions(GWASCat,RegCat,byReg1="gene.name",byReg2="GeneName")
#' OverlapRegions(Benchmark,RegCat)
#' OverlapRegions(RegCat,Benchmark)
#' OverlapRegions(GWASCat,Benchmark)
#'
#' RegCatFull <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Pooled_Analysis_Results_Rdmp_23022017/Pooled_MA_LM_EnhGene_DF_FULL_ADJUSTED_2017-02-28Benchmarking12017-03-02short.rds")
#' BenchmarkFull <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Data/Benchmark1/Benchmark1_Dataset_PooledANDProcessed_5sources_17_01_18.rds")
#'
#' Benchmarked.int <- OverlapRegions(,
#'                                  apply(BenchmarkFull[1:1000,],2,as.character))
#' RegCat=apply(RegCatFull[1:100000,],2,as.character)
#' Bench=apply(BenchmarkFull[1:100000,],2,as.character)
#' Benchmarked.int <- OverlapRegions(Reg1,Reg2)


OverlapRegions <- function(Reg1,Reg2,byReg1=NULL,byReg2=NULL){
  
  
  require(GenomicRanges)
  require(stringr)
  # need to add part where input df will be organized as dataframes
  
  Coordinates1 <- c("chr1","start1","end1")
  Coordinates2 <- c("chr2","start2","end2")
  
  if (!min(Coordinates1%in%colnames(Reg1),
           Coordinates1%in%colnames(Reg2))) {stop("Coordinates are missing, 
                                                  eg colnames chr1,start1,
                                                  end1 are missing !!!")}
  
  # Gr obj
  Reg1.gr <- GRanges(Reg1[,"chr1"],IRanges(as.integer(Reg1[,"start1"]),
                                           as.integer(Reg1[,"end1"])))
  
  Reg2.gr <- GRanges(Reg2[,"chr1"],IRanges(as.integer(Reg2[,"start1"]),
                                           as.integer(Reg2[,"end1"])))
  
  # Reg1-Reg1 overlap            
  Reg1.Reg1.Overlap <- data.frame(findOverlaps(Reg1.gr,Reg2.gr))
  
  colnames(Reg1) <- paste0("Reg1_",colnames(Reg1))
  colnames(Reg2) <- paste0("Reg2_",colnames(Reg2)) 
  
  Reg1Reg2 <- cbind(Reg1[Reg1.Reg1.Overlap$queryHits,],
                    Reg2[Reg1.Reg1.Overlap$subjectHits,])
  
  
  # filter by columns added if requested
  
  if (!(is.null(byReg2)&is.null(byReg1))) {
    
    # select columns to filter
        Reg1FilterC <- which(str_detect(colnames(Reg1Reg2),
                                        byReg1)&str_detect(colnames(Reg1Reg2),"Reg1"))
        Reg2FilterC <- which(str_detect(colnames(Reg1Reg2),
                                        byReg2)&str_detect(colnames(Reg1Reg2),"Reg2"))
        
    if ((length(Reg1FilterC)==0)|(length(Reg2FilterC)==0)){stop(
      'byReg1 or byReg2 arguments were not entered correctly')}
    
    Reg1Reg2 <- Reg1Reg2[which(Reg1Reg2[,Reg1FilterC]==Reg1Reg2[,Reg2FilterC]),]
    
  }
  
  
  # overlap Reg2 with Reg2
  
  
  if ((is.null(byReg2)&is.null(byReg1))&(sum(str_detect(colnames(Reg1Reg2), 
                                                        paste0(Coordinates2,
                                                               collapse="|")))==6)) {
    
    # GR
    Reg1.gr2 <- GRanges(Reg1Reg2[,"Reg1_chr2"],
                        IRanges(as.integer(Reg1Reg2[,"Reg1_start2"]),
                                as.integer(Reg1Reg2[,"Reg1_end2"])))
    Reg2.gr2 <- GRanges(Reg1Reg2[,"Reg2_chr2"],
                        IRanges(as.integer(Reg1Reg2[,"Reg2_start2"]),
                                as.integer(Reg1Reg2[,"Reg2_end2"])))
    
    
    Reg2.Reg2.Overlap2 <- data.frame(findOverlaps(Reg1.gr2,Reg2.gr2))
    
    Reg2Reg2 <- Reg1Reg2[Reg2.Reg2.Overlap2$queryHits,]
    
    # IMORTANT! removing all Region1 entries BOTH coord1 and coord2 
    # overlap the EITHER same Reg2 coord1 OR Reg2 coord1  
    
    
    Reg1.gr2 <- GRanges(Reg2Reg2[,"Reg1_chr2"],
                        IRanges(as.integer(Reg2Reg2[,"Reg1_start2"]),
                                as.integer(Reg2Reg2[,"Reg1_end2"])))
    Reg2.gr2 <- GRanges(Reg2Reg2[,"Reg2_chr1"],
                        IRanges(as.integer(Reg2Reg2[,"Reg2_start1"]),
                                as.integer(Reg2Reg2[,"Reg2_end1"])))
    
    Reg2.Reg1overlap <- data.frame(findOverlaps(Reg1.gr2,Reg2.gr2))
    
    if (nrow(Reg2.Reg1overlap)!=0) {Reg2Reg2 <- Reg2Reg2[-Reg2.Reg1overlap$queryHits,]}
    
    
    Reg1.gr2 <- GRanges(Reg2Reg2[,"Reg1_chr1"],
                        IRanges(as.integer(Reg2Reg2[,"Reg1_start1"]),
                                as.integer(Reg2Reg2[,"Reg1_end1"])))
    Reg2.gr2 <- GRanges(Reg2Reg2[,"Reg2_chr2"],
                        IRanges(as.integer(Reg2Reg2[,"Reg2_start2"]),
                                as.integer(Reg2Reg2[,"Reg2_end2"])))
    
    Reg2.Reg1overlap <- data.frame(findOverlaps(Reg1.gr2,Reg2.gr2))
    
    if (nrow(Reg2.Reg1overlap)!=0) {Reg2Reg2 <- Reg2Reg2[-Reg2.Reg1overlap$queryHits,]}
    
    Reg1Reg2 <- Reg2Reg2
    
  }
  
  return(Reg1Reg2)
  
  }



