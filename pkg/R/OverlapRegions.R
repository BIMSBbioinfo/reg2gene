#' Overlap regions and/or interacting regions and/or their metadata  
#' 
#' The function that compares two input regions or a combinations of two
#' interacting regions in the genome and reportes regions present in both
#' input datasets. If one coordinate set per input dataframe is reported
#' ("chr1","start1","end1") then regions from one set, Reg1, is identified
#' in the other dataset,Reg2, and repoted.
#' If two coordinate sets ("chr1","start1","end1","chr2","start2","end2")
#' per input dataframes are reported (eg two interacting regions in test
#' and benchmark dataset) then both of them are compared,and if both 
#' coordinates regions overlap they are reported.
#' Third, an overlap of 2 dataframes can be identified based on the
#' overlap between genomic coordinates of two regions (regions originated from
#' two different datasets) and meta-data column (useful for testing overlaps
#' with eQTL studies based on the gene name). 
#'   
#'       
#' @param Reg1 - character matrix with at least three columns named:
#' "chr1","start1","end1" which store information about coordinates
#' of the regions which we want to comape. Additionally,
#' "chr2","start2","end2" columns or any other column in the
#' matrix (such as SNP ID, proxy SNP ID,or reported gene).
#'   
#' @param Reg2 - character matrix with at least three columns named:
#' "chr1","start1","end1" which store information about coordinates
#' of the regions which we want to compare Additionally,
#' "chr2","start2","end2" columns or any other column in the
#' matrix (such as SNP ID, proxy SNP ID,or reported gene).
#'  
#'   
#' @param byReg1 - character (default NULL),a name of the column which
#' will be used for additional filtering of the overlapping regions (after
#' region1 is overlapped with region2).
#'   
#' @param byReg2 - character (default NULL),a name of the column which
#' will be used for additional filtering of overlapping regions (after
#' region1 is overlapped with region2). 
#'  
#'   
#'  
#' @return A dataframe of coordinates of all overlapping regions.
#' Each row corresponds to the combination of overlapped rows from 2 input
#' dataframes based on either: a) simple coordinates overlap, b) coordinates overlap + 
#' additionaly filtered based on byReg1 and byReg2 defined columns overlap, c) dual 
#' coordinates overlap between interacting regions. Column names indicate whether 
#' coordinates and additional info originates from Reg1 or Reg2alog.
#'  
#'  
#' @import GenomicRanges
#' @import stringr
#'  
#' @details Function that performs:a) simple overlap between 2 
#' dataframes of coordinates stored as  "chr1","start1","end1". It creates
#' GRanges object, and then it runs \code{\link[GenomicRanges]{findOverlaps}}.
#' If input dataframe contains 2 locations of interacting regions:
#' "chr1","start1","end1","chr2","start2","end2" then it compares:
#' (Reg1coord1 with Reg2coord1)&(Reg1coord2 with Reg2coord2), and returns rows
#' where both combinations of overlaps are TRUE. If byReg1 and byReg1 different 
#' than NULL then overlapped regions are additionally filtered based on the
#' overlap between byReg1 and byReg1 defined columns
#' 
#' @example  load("~/GWASCat.RData")
#' load("~/RegCat.RData")
#'
#' OverlapRegions(GWASCat,RegCat)
#' OverlapRegions(GWASCat,RegCat,byReg1="GeneName",byReg2="gene.name")
#' error function:OverlapRegions(GWASCat,RegCat,byReg1="gene.name",byReg2="GeneName")
#'
#'
#'
#' @export




OverlapRegions <- function(Reg1,Reg2,byReg1=NULL,byReg2=NULL){
  
  
  
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
                                        paste0(byReg1,"$"))&str_detect(colnames(Reg1Reg2),"Reg1"))
        Reg2FilterC <- which(str_detect(colnames(Reg1Reg2),
                                        paste0(byReg2,"$"))&str_detect(colnames(Reg1Reg2),"Reg2"))
        
    if ((length(Reg1FilterC)==0)|(length(Reg2FilterC)==0)){stop(
      'byReg1 or byReg2 arguments were not entered correctly')}
    
    Reg1Reg2 <- Reg1Reg2[which(Reg1Reg2[,Reg1FilterC]==Reg1Reg2[,Reg2FilterC]),]
    
  }
  
  
  # overlap Reg2 with Reg2
  
  
  if ((is.null(byReg2)&is.null(byReg1))&(sum(str_detect(colnames(Reg1Reg2), 
                                                        paste0(Coordinates2,
                                                               collapse="|")))==6)) {
    
    # GR
    # Gr obj
    Reg1.gr <- GRanges(Reg1[,"Reg1_chr2"],IRanges(as.integer(Reg1[,"Reg1_start2"]),
                                             as.integer(Reg1[,"Reg1_end2"])))
    
    Reg2.gr <- GRanges(Reg2[,"Reg2_chr2"],IRanges(as.integer(Reg2[,"Reg2_start2"]),
                                             as.integer(Reg2[,"Reg2_end2"])))
    
    # Reg1-Reg1 overlap            
    Reg2.Reg2.Overlap <- data.frame(findOverlaps(Reg1.gr,Reg2.gr))
    
    # if overlap is duplicated then reg1 overlaps reg1, and reg2 overlaps reg2 -> confirmed
    AllOverlapRegions <- rbind(Reg2.Reg2.Overlap,Reg1.Reg1.Overlap)
    AllOverlapRegions.dpl <- AllOverlapRegions[duplicated(AllOverlapRegions),]
    
    
    Reg1Reg2 <- cbind(Reg1[AllOverlapRegions.dpl$queryHits,],
                      Reg2[AllOverlapRegions.dpl$subjectHits,])
    

    
  }
  
  return(Reg1Reg2)
  
  }


#' Criss-cross overlap of interacting regions in the genome 
#' 
#'  A function that compares four input regions,eg two interacting regions
#'  from one dataset with two interacting regions from the other dataset
#'  in the criss-cross manner. Thus, the input orientation of the 
#'  interacting regions is not strict: the first pair of the region1
#'  from one dataset will be compared to both; the first pair and the
#'  second pair of region2(and its pair is compared in the opposite 
#'  direction). Two coordinate sets ("chr1","start1","end1",
#'  "chr2","start2","end2") per input dataframe are MUST.  
#'       
#'  @param Reg1 - character matrix with at least six columns named
#'  "chr1","start1","end1","chr2","start2","end2" which store 
#'  information about coordinates of two interaction regions in the
#'  genome which we want to comape to two other interacting regions.
#'   
#'  @param Reg2 - character matrix with at least six columns named
#'  "chr1","start1","end1","chr2","start2","end2" which store 
#'  information about coordinates of two interaction regions in the
#'  genome which we want to comape to two other interacting regions.
#'
#'  @return A dataframe of coordinates of all overlapping regions and
#'  metadata that was stored in the input matrices. Each row corresponds
#'  to one combination of Region1-Region2 overlaps. Column names indicate
#'  whether coordinates and additional info come from Reg1 or Reg2alog.
#'  
#'  @import GenomicRanges
#'  @import stringr
#'  
#'  @details Function performs criss-cross overlap between 2 dataframes 
#'  with coordinates from two interacting regions: "chr1","start1","end1",
#'  "chr2","start2","end2". It creates GRanges objects, and performs
#'  \code{\link[reg2gene]{OverlapRegions}} overlaps between 
#'  Reg1Coord1-Reg2Coord1&Reg1Coord2-Reg2Coord2 and vice-versa
#'  Reg1Coord1-Reg2Coord2&Reg1Coord2-Reg2Coord1.
#' 
#'  @example
#'  load(file="~/RegCat.RData")
#'  load(file="~/ComplexReg2.RData")
#'
#'  ComplexOverlaps(RegCat,ComplexReg2)
#'  
#'  
#' @export


ComplexOverlaps <- function(Reg1,Reg2){
  
  # criss-cross overlap when orientation of overlaps is unknown
      Reg11Reg22 <- OverlapRegions(Reg1,Reg2)
            # switch names in Reg2 dataframe
                colnames(Reg1)[1:6] <- c("chr1","start1","end1","chr2","start2","end2")
      Reg12Reg21 <- OverlapRegions(Reg1,Reg2)
  
              colnames(Reg11Reg22) <- colnames(Reg12Reg21)
     
     Reg1Reg2 <- rbind(Reg11Reg22,Reg12Reg21)
     
     
     return(Reg1Reg2)
     
     
}




  
