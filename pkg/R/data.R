###################################################
# Examples for Overlap functions

#' Sample file for Region1 imported as GRanges
#'
#'
#' @return GRanges object
#' @format GRanges object with 6 ranges and 3 metadata columns
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
#' @return GRanges
#' @format GRanges object with 14 ranges and 3 metadata columns
"GRReg2_toy"



#' Sample GRanges to test confusionmatrix function
# 
# This dataset contains a set of genomic regions stored in a dataframe format
# Thus, it contains following columns: chr1,start1,end1,chr2,start2,end2.
# Regions are defined in such way that 1st and 2nd row per each chromosome
# should overlap with the toy example in Reg1Extended_toy. All others should not.
# eg They are selected such that Coord1/Reg1 and Coord1/Reg2 overlap with each other,
# but Coord2/Reg1 and Coord2/Reg2 DO NOT overlap (and vice-versa).
# Thus,  4 rows of Reg1Extended_toy finds a matching pairs in Reg2Extended_toy.
# ComplexOverlaps() should report 2 overlaps since the 2nd coordinates-pair of Reg2 is a
# flipped version of the 1st coordinates-pair of Reg2. ExpectedOverap column indicates what
# kind of overlap is expected to observe.
#'
#' @return data.frame object
#' @formatGRanges object with 6 ranges and 5 metadata columns
"BenchMarkedReg2Gene_toy"





