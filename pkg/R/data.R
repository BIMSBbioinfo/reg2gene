###################################################
# Examples for various functions from reg2gene


#' GRanges example to test various functions from reg2gene
#' 
#' @details Genome locations are preselected in such way that all possible outputs of
#' the  benchmarking procedure is captured. *Colums which correspond to the 
#' expected outcome of the benchmarking procedure are:Bench1Exp,Filter2Exp,
#' Bench2Exp,Filter1Exp
#' Where Bench1 andFilter1 corresponds to the benchmarking of GRReg1_toy with 
#' GRReg2_toy, whereas Bench2 and Filter2 corresponds to the benchmarking of 
#' GRReg1_toy with itself.
#' For example,following object GRReg1_toy[7] should be benchmarked 3 times 
#' with  GRReg2_toy dataset
#' 
#' @format GRanges object with 9 ranges and 7 metadata columns
#' @docType data
"GRReg1_toy"


#' Sample file for benchmark dataset
#'
#'
#' @format GRanges object with 14 ranges and 3 metadata columns
#' @docType data
"GRReg2_toy"



###################################################
# Create testing matrix of variables with predefined corr:
# r_vector = c(0.99,0.9,0.8,0.6,0.4,0.3,0.1,0)


#' GRanges example to test modelling of associateReg2gene function from reg2gene
#' 
#' @details cor(y~x1)=0.99; cor(y~x2)=0.9,cor(y~x3)=0.8,...,cor(y~x8)=0
#' Drop in correlation can be used to test results of modelling procedure. 
#' 52 "celltypes" included to approximate for 52 H3K4me1&RNASeq datasets from 
#' Roadmap
#' 
#' @format GRanges object with 9 ranges and 55 metadata columns:
#' @docType data
"ModellingTest"






