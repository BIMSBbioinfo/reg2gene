# ---------------------------------------------------------------------------- 
# Unitests for query functions

#########################
# Create test datasets

library(InteractionSet)
library(testthat)
library(GenomicRanges)
library(genomation)
library(reg2gene)

# Create genomicRegions
     genomicRegions <- GRanges(c("chr1:1-2", # 1. overlap prom
                                 "chr2:1-2",  # 2. overlap enh
                                 "chr3:1-2", # 3. overlap tss +/- 1,000,000
                                 "chr4:1-2")) # 4. do not overlap tss +/- 1Mb
    
     annotationsEnh <- GRanges(c("chr1:1-2",
                                 "chr2:1-2",
                                 "chr3:100000-100002",
                                 "chr4:10000001-10000002"))
    
     annotationsGenes <- GRanges(c("chr1:1-2",
                                   "chr2:100000-100002",
                                   "chr3:99999-100002",
                                   "chr4:10000001-10000002"))
    
     seqlengths(annotationsEnh) <- seqlengths(annotationsGenes) <- rep(10000002,
                                                                       4)
     annotations = GInteractions(annotationsEnh,annotationsGenes,
                                 name=c("gen1","gen2","gen3","gen4"))
     
     annotationsGenes$name=annotations$name

     
################################################
####### TEST annotateGenomicRegions     
 
 
 test_that("annotateGenomicRegions INPUT/OUTPUT is correct in form",{
   
   # TEST INPUT:
   expect_true(class(annotations)%in%c("GRanges","GInteractions"))
   expect_is(genomicRegions,"GRanges")
   
   # test expected input meta-data:
   expect_equal(colnames(mcols(annotations)),"name") 
   
    # TEST OUTPUT - QFR
   # is GRanges with lenght>0
   
   QFR <- annotateGenomicRegions(genomicRegions=genomicRegions,
                       annotations=annotations)
  expect_is(QFR,"GInteractions")
  expect_equal(names(mcols(QFR)),c("name","annotatedAs"))
   
   
  # if annotation is GRanges object - works ok
   QFR <- annotateGenomicRegions(genomicRegions=genomicRegions,
                              annotations=annotationsGenes) 
   expect_is(QFR,"GInteractions")
   expect_equal(names(mcols(QFR)),c("name","annotatedAs"))
 })
   
 
 test_that("annotateGenomicRegions arguments are correctly used",{
   
   # if annotation is GRanges object - works ok
   QFR <- annotateGenomicRegions(genomicRegions=genomicRegions,
                              annotations=annotationsGenes) 
   expect_equal(length(QFR),1)
   expect_equal(QFR$name,"gen1")
   expect_equal(second(QFR),
                reduce(trim(suppressWarnings(promoters(annotationsGenes,
                                              upstream=1000,
                                              downstream=1000)))[1]))
   
   # identified=T
   # TEST OUTPUT - QFRIT
           QFRIT <-  annotateGenomicRegions(genomicRegions,
                             annotations =annotations,
                             identified=T)
         
         expect_equal(length(QFRIT),3) # 3 ranges should be identified
         expect_equal(QFRIT$name,c("gen1","gen2","gen3"))
         # check annotatedAs argument
         expect_equal(QFRIT$annotatedAs,c("promoter","enhancer","nearestGene"))
         expect_false("gen4"%in%QFRIT$name)
         
   # location of result is equal to input
           expect_equal(second(QFRIT),
                        reduce(trim(suppressWarnings(promoters(annotationsGenes,
                                upstream=1000, downstream=1000)))[1:3]))
   # identified=F
         QFRIF <- annotateGenomicRegions(genomicRegions,
                             annotations,
                             identified=F)
      
         expect_equal(length(QFRIF),1) # 1 ranges should be identified
         # location of result is equal to input 4th genomicRegions
         expect_equal(QFRIF, genomicRegions[4])
         
   
   # distance
    # increase distance works well
         
         QFRD <-  annotateGenomicRegions(genomicRegions,
                                       annotations =annotations,
                                       distance = 10000000,
                                       identified=T)
         expect_equal(length(QFRD),4) # 3 ranges should be identified
         expect_equal(QFRD$name,c("gen1","gen2","gen3","gen4"))
        
         
    # reduce distance works well   
         QFRD <-  annotateGenomicRegions(genomicRegions,
                                      annotations=annotations,
                                      distance = 100,
                                      identified=T)
         expect_equal(length(QFRD),2) # 3 ranges should be identified
         expect_equal(QFRD$name,c("gen1","gen2"))
         
   # upstream,downstream   [+/-1bp]
         
         QFRDU <-  annotateGenomicRegions(genomicRegions,
                                      annotations=annotations,
                                      downstream = 1,
                                      upstream =  1,
                                      distance = 100)
         
         expect_equal(length(QFRDU),2) # 2 ranges should be identified
         expect_equal(QFRDU$name,c("gen1","gen2"))
         expect_equal(second(QFRDU),
                      reduce(trim(suppressWarnings(promoters(annotationsGenes,
                                                      upstream=1,
                                                      downstream=1)))[1:2]))
         
         # but if distance is 1MB,then the 3rd gene should be identified as well
         expect_equal( annotateGenomicRegions(genomicRegions,
                                           annotations=annotations,
                                           downstream = 1,
                                           upstream =  1)$name,c("gen1",
                                                                  "gen2",
                                                                  "gen3"))
         
 })
         