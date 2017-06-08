#' Selects Reg2gene entries which could be benchmarked with Benchmark dataset
#' 
#' The function that takes as an input results of \code{linkReg2Gene} or
#' any other modelling procedure implemente in \code{reg2gene} package,
#' such as \code{MetaReg2Gene} and predefined Benchmark dataset.
#' For input Reg2Gene object it adds a column with info whether reported interactions is 
#' benchmarked or not. By default it reportes how many times Reg2Gene interactions is
#' observed in the benchmark dataset. If binary is set to TRUE, then only logical vector of
#' TRUE (overlapping benchmark dataset at least once) and FALSE(not overlapping benchmark 
#' dataset at all) entries is reported.
#' 
#' 
#'  @param Reg2GenePreFilter a GRanges object output from \code{linkReg2Gene} or \code{MetaReg2Gene}
#'  function or a character matrix containing at least six columns named:
#'  "chr1","start1","end1","chr2","start2","end2" and where "chr1","start1","end1" corresponds
#'  to regulatory regions whereas "chr2","start2","end2" correpond to TSS of a gene.
#'  
#'  @param Benchmark character matrix with at least three columns named:
#'  "chr1","start1","end1" which store information about coordinates
#'  of the regions which we want to compare and 1) one additional column 
#'  with gene names, or 2)  "chr2","start2","end2" columns
#'   
#' 
#' @param byReg2Gene character (def:NULL), name of the column where gene names are stored 
#' in Reg2Gene object. Optional and useful if filtering should be performed based on the overlap
#' between regions and gene names. 
#' @param byBenchmark character (def:NULL), name of the column where gene names are stored 
#' in Benchmark object. Optional and useful if filtering should be performed based on the overlap
#' between regions and gene names.

library(plyr)
library(stringr)
library(GenomicRanges)
library(reg2gene)


Filter_PreBench <- function(Reg2GenePreFilter,
                            Benchmark,
                            byReg2Gene=NULL,
                            byBenchmark=NULL){
  
      if (class(Reg2GenePreFilter)=="GRanges") { Reg2GenePreFilter <- .ReArrangeGR(GR=Reg2GenePreFilter)}
  

      # setting both Benchmark coordinates on top of one another  
       
         MinNames <- c("chr1","start1","end1")
              Benchmark.AllRegions <- rbind(Benchmark[,MinNames],setNames(Benchmark[,str_replace(MinNames,"1","2")],MinNames))
          # filtering overlap of Coord1
              Reg2GenePreFilter.Coord1 <- OverlapRegions(Reg2GenePreFilter,Benchmark.AllRegions)
              
          # filtering overlap of Coord2, after filtering of overlap Coord1 already performed 
                  Reg2GenePreFilter.Coord1 <- rename(Reg2GenePreFilter.Coord1,"Reg1_chr2"="chr1","Reg1_start2"="start1","Reg1_end2"="end1")
                  Reg2GenePreFilter.Coord2 <- OverlapRegions(Reg2GenePreFilter.Coord1,Benchmark.AllRegions)
          
          #rearranging the format          
                    Reg2GeneFiltered <- Reg2GenePreFilter.Coord2[,1:ncol(Reg2GenePreFilter)]
                          colnames(Reg2GeneFiltered) <- colnames(Reg2GenePreFilter)
                    Reg2GeneFiltered <- unique(Reg2GeneFiltered)
                    
      return(Reg2GeneFiltered)
                    
          
          }
  




Filter_PreModelling <- function()



Filter_PreMetA <- function()
