.switchReg = function(reg, reg.col='reg'){
  vals = values(reg)
  cols = vals[,colnames(vals) != reg.col]
  reg2 = values(reg)[[reg.col]]
  
  
  values(reg2) = cbind(DataFrame(granges(reg)), cols)
  colnames(values(reg2))[1] = reg.col
  reg2
  
}


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
#' GRanges object, and then it runs \code{findOverlaps}.
#' If input dataframe contains 2 locations of interacting regions:
#' "chr1","start1","end1","chr2","start2","end2" then it compares:
#' (Reg1coord1 with Reg2coord1)&(Reg1coord2 with Reg2coord2), and returns rows
#' where both combinations of overlaps are TRUE. If byReg1 and byReg1 different 
#' than NULL then overlapped regions are additionally filtered based on the
#' overlap between byReg1 and byReg1 defined columns
#' 
#' @examples  OverlapRegions(Reg1_toy, Reg2_toy)
#' OverlapRegions(Reg1_toy[2,], Reg2_toy)
#' 
#'
#'
#'
#' @export
OverlapRegions <- function(Reg1,Reg2){


   # Reg1-Reg1 overlap            
  Reg1.Reg1.Overlap <- data.frame(findOverlaps(Reg1,Reg2))
  
  # adjusting names
      colnames(mcols(Reg1)) <- paste0("Reg1_",colnames(mcols(Reg1)))
      colnames(mcols(Reg2)) <- paste0("Reg2_",colnames(mcols(Reg2)))
   
  # adding Overlapped regions and mcols on top of those we want to test       
    Reg1Reg2 <- Reg1[Reg1.Reg1.Overlap$queryHits,]
    Reg1Reg2$Reg2coord1 <- (Reg2[Reg1.Reg1.Overlap$subjectHits,])
         mcols(Reg1Reg2) <- cbind(mcols(Reg1Reg2),mcols(Reg2[Reg1.Reg1.Overlap$subjectHits,]))
    
    ########## 
    # Reg1-Reg1 overlap            
    Reg2.Reg2.Overlap <- data.frame(findOverlaps(Reg1Reg2$Reg1_reg,Reg1Reg2$Reg2_reg))
    
    # if overlap is duplicated then reg1 overlaps reg1, and reg2 overlaps reg2 -> confirmed
    AllOverlapRegions <- which(Reg2.Reg2.Overlap$queryHits == Reg2.Reg2.Overlap$subjectHits)
    InputRows.Overlapped <- Reg2.Reg2.Overlap[AllOverlapRegions,"queryHits"]
    Reg1Reg2 <- Reg1Reg2[InputRows.Overlapped,]
    Reg1Reg2$Coor1Coord2PAIR <- Reg1.Reg1.Overlap$queryHits[InputRows.Overlapped]
    
  
  return(Reg1Reg2)
  
  }



#' Criss-cross overlap of interacting regions 
#' 
#' A function that compares four input regions,eg two interacting regions
#' from one dataset with two interacting regions from the other dataset
#' in the criss-cross manner. Thus, the input orientation of the 
#' interacting regions is not strict: the first pair of the region1
#' from one dataset will be compared to both; the first pair and the
#' second pair of region2(and its pair is compared in the opposite 
#' direction). Two coordinate sets ("chr1","start1","end1",
#' "chr2","start2","end2") per input dataframe are MUST.  
#'       
#' @param Reg1 - character matrix with at least six columns named
#' "chr1","start1","end1","chr2","start2","end2" which store 
#' information about coordinates of two interaction regions in the
#' genome which we want to comape to two other interacting regions.
#'   
#' @param Reg2 - character matrix with at least six columns named
#' "chr1","start1","end1","chr2","start2","end2" which store 
#' information about coordinates of two interaction regions in the
#' genome which we want to comape to two other interacting regions.
#'
#' @return A dataframe of coordinates of all overlapping regions and
#' metadata that was stored in the input matrices. Each row corresponds
#' to one combination of Region1-Region2 overlaps. Column names indicate
#' whether coordinates and additional info come from Reg1 or Reg2alog.
#'  
#' @import GenomicRanges
#' @import stringr
#'  
#' @details Function performs criss-cross overlap between 2 dataframes 
#' with coordinates from two interacting regions: "chr1","start1","end1",
#' "chr2","start2","end2". It creates GRanges objects, and performs
#' \code{\link[reg2gene]{OverlapRegions}} overlaps between 
#' Reg1Coord1-Reg2Coord1&Reg1Coord2-Reg2Coord2 and vice-versa
#' Reg1Coord1-Reg2Coord2&Reg1Coord2-Reg2Coord1.
#' 
#' @examples ComplexOverlaps(GRReg1_toy, GRReg2_toy) 
#' ComplexOverlaps(GRReg1Extended_toy,GRReg2Extended_toy)
#'ComplexOverlaps(Reg1=GRReg1Extended_toy.2,Reg2=GRReg2Extended_toy) 
#' @export
ComplexOverlaps <- function(Reg1,Reg2){
  
  # criss-cross overlap when orientation of overlaps is unknown
      Reg11Reg22 <- OverlapRegions(Reg1,Reg2)
            # switch names in Reg1 dataframe

      Reg12Reg21 <- OverlapRegions(Reg1,.switchReg(Reg2))
     
  
    return(c(Reg11Reg22,Reg12Reg21))


}





  
