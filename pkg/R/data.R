#' ###################################################
#' # Examples for Overlap functions
#' 
#' #' Sample file for Region1 imported as data frame
#' #' 
#' #' This dataset contains a set of genomic regions stored in a dataframe format
#' #' Thus, it contains following columns: chr1,start1,end1,chr2,start2,end2. 
#' #' Regions are selected in such way that 1st row will overlap with the toy example
#' #' Reg2_toy, whereas 2nd pair of coordinates should not find a matching pair in
#' #' Reg1_toy.
#' #' 
#' #' @return data.frame object
#' #' @format A data frame with 2 obs. and  6 variables
#' "Reg1_toy"
#' 
#' 
#' #' Sample file for Region2 imported as data frame
#' #' 
#' #' This dataset contains a set of genomic regions stored in a dataframe format
#' #' Thus, it contains following columns: chr1,start1,end1,chr2,start2,end2. 
#' #' Regions are selected in such way that 1st and 2nd row should overlap with the toy example
#' #' Reg1_toy. All others should not. They are selected in a manner that Coord1 of Reg1 and
#' #' Reg2 overlap with each other, but Coord2 of Reg1 and Reg2 DO NOT overlap (and vice-versa).
#' #' Thus, only 1 row of Reg1 finds a matching pairs in Reg2. 
#' #' ComplexOverlaps() should report 2 overlaps since the 2nd coordinates-pair of Reg2 is a 
#' #' flipped version of the 1st coordinates-pair of Reg2. ExpectedOverap column indicates what
#' #' kind of overlap is expected to observe.
#' #' 
#' #' @return data.frame object
#' #' @format A data frame with 7 obs. and  7 variables
#' "Reg2_toy"
#' 
#' 
#' #' Sample file for Extended Region1 imported as data frame
#' #' 
#' #' This dataset contains a set of genomic regions stored in a dataframe format
#' #' Thus, it contains following columns: chr1,start1,end1,chr2,start2,end2. 
#' #' Regions are selected in such way that every second row (1,3,5,7) should overlap 
#' #' with the toy examples in Reg2Extended_toy, whereas other pairs of coordinates should 
#' #' not find a matching pair in Reg2Extended_toy. Same regions as Reg1_toy, but new
#' #' chromosomes are added
#' #' 
#' #' @return data.frame object
#' #' @format A data frame with 8 obs. and  6 variables
#' "Reg1Extended_toy"
#' 
#' 
#' #' Sample file for Extended Region2 imported as data frame
#' #' 
#' #' This dataset contains a set of genomic regions stored in a dataframe format
#' #' Thus, it contains following columns: chr1,start1,end1,chr2,start2,end2. 
#' #' Regions are defined in such way that 1st and 2nd row per each chromosome 
#' #' should overlap with the toy example in Reg1Extended_toy. All others should not. 
#' #' eg They are selected such that Coord1/Reg1 and Coord1/Reg2 overlap with each other, 
#' #' but Coord2/Reg1 and Coord2/Reg2 DO NOT overlap (and vice-versa).
#' #' Thus,  4 rows of Reg1Extended_toy finds a matching pairs in Reg2Extended_toy. 
#' #' ComplexOverlaps() should report 2 overlaps since the 2nd coordinates-pair of Reg2 is a 
#' #' flipped version of the 1st coordinates-pair of Reg2. ExpectedOverap column indicates what
#' #' kind of overlap is expected to observe.
#' #' 
#' #' @return data.frame object
#' #' @format A data frame with 28 obs. and  7 variables
#' "Reg2Extended_toy"
#' 
#' 
#' 
