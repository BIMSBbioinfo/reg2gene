###################################################
# Examples for Overlap functions

#' Sample file for Region1 imported as GRanges
#'
#'
#' @format GRanges object with 6 ranges and 3 metadata columns
#' @docType data
"GRReg1_toy"


#' Sample file for Region2 imported as data frame
#'
# This dataset contains a set of genomic regions stored in a dataframe format
# Thus, it contains following columns: chr1,start1,end1,chr2,start2,end2.
# Regions are selected in such way that 1st and 2nd row should overlap with the toy example
# Reg1_toy. All others should not. They are selected in a manner that Coord1 of Reg1 and
# Reg2 overlap with each other, but Coord2 of Reg1 and Reg2 DO NOT overlap (and vice-versa).
# Thus, only 1 row of Reg1 finds a matching pairs in Reg2.
#' ComplexOverlaps() should report 2 overlaps since the 2nd coordinates-pair of Reg2 is a
#' flipped version of the 1st coordinates-pair of Reg2. ExpectedOverap column indicates what
#' kind of overlap is expected to observe.
#'
#' @format GRanges object with 14 ranges and 3 metadata columns
#' @docType data
"GRReg2_toy"










