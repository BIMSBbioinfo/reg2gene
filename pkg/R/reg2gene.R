#' reg2gene: An R package for predicting the target genes for regulatory elements
#'
#'
#'
#' @docType package
#' @name reg2gene
#' 
#' 
#' @import genomation
#' @import GenomicRanges
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom parallel mclapply
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @import stringr
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom caret confusionMatrix
#' @importFrom caret posPredValue