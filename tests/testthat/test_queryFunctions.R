# ---------------------------------------------------------------------------- 
# Unitests for query functions

#########################
# Create test datasets

library(InteractionSet)
library(testthat)
library(GenomicRanges)
library(genomation)
library(reg2gene)

# Create windows
     windows <- GRanges(c("chr1:1-2", # 1. overlap prom
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
     interactions = GInteractions(annotationsEnh,annotationsGenes,
                                 name=c("gen1","gen2","gen3","gen4"))
     
     
     geneAnnotations=second(interactions)
     mcols(geneAnnotations) <- mcols(interactions) 
     
     
     # Planned combinations
     # 1) intset NULL, geneA 1 #expO: gen1
     # 2) intset full, geneA 1 #expO: gen1,gen2
     # 3) intset full, geneA full #expO: gen1,gen2,gen3
     # 4) intset full, geneA 1&3 #expO: gen1,gen2,gen3
     # 5) intset full -2 , geneA full #expO: gen1,gen3
     # 6) intset 4, geneA 4 no overlap identified
     
     # Combinations when intset does not contain info about gene, whereas
     # it needs to be infered from geneAnnotations
     
     # 7) intset full, geneA 1
     # 8) intset full, geneA full
     # 9) ? intset full, geneA full
     
################################################
####### TEST annotatewindows     
 
 
 test_that("annotatewindows INPUT/OUTPUT is correct in form",{
   
   # TEST INPUT:
   expect_true(class(windows)=="GRanges")
   expect_true(class(geneAnnotations)=="GRanges")
   expect_true(class(interactions)=="GInteractions")
  
   # test expected input meta-data:
   expect_equal(colnames(mcols(interactions)),"name") 
   expect_equal(colnames(mcols(geneAnnotations)),"name") 
   
    # TEST OUTPUT - QFR
         # is GRanges with lenght>0
         
         QFR <- reg2gene(windows = windows,
                         interactions=interactions,
                         geneAnnotations = geneAnnotations)
        
         expect_is(QFR,"GInteractions")
         expect_equal(names(mcols(QFR)),c("name","annotatedAs"))
   
   
  # if interactions is NULL is GRanges object - works ok
         QFR <- reg2gene(windows = windows,
                         interactions=NULL,
                         geneAnnotations = geneAnnotations)
            expect_is(QFR,"GInteractions")
            expect_equal(names(mcols(QFR)),c("name","annotatedAs"))
 
            
            # if interactions have no info about gene - just geneAnnot have info
            
            interactions2 <- interactions
            mcols(interactions2) <- NULL
            
            QFR <- reg2gene(windows = windows,
                            interactions=interactions2,
                            geneAnnotations = geneAnnotations,
                            annotateInteractions = TRUE)
            
            expect_is(QFR,"GInteractions")
            expect_equal(names(mcols(QFR)),c("name","annotatedAs"))     
            
            
            })
   
 
 test_that("reg2gene arguments are correctly used",{
   
   # 1) intset NULL, geneA 1 #expO: gen1
   QFR <- reg2gene(windows = windows,
                   interactions = NULL,
                   geneAnnotations = geneAnnotations[1]) 
   
   expect_equal(length(QFR),1)
   expect_equal(QFR$name,"gen1")
   expect_equal(granges(second(QFR)),granges(geneAnnotations[1]))
       
             # checking identified argument as false
             expect_equal(reg2gene(windows = windows,
                      interactions = NULL,
                      geneAnnotations = geneAnnotations[1],
                      identified = FALSE), windows[-1])
       
   
   # 2) intset full, geneA 1 #expO: gen1,gen2
   QFR <-  reg2gene(windows,
                    interactions = interactions,
                    geneAnnotations = geneAnnotations[1])
   
         expect_equal(length(QFR),2) # 3 ranges should be identified
         expect_equal(QFR$name,c("gen1","gen2"))
         # check annotatedAs argument
         expect_equal(QFR$annotatedAs,c("promoter","enhancer"))
         expect_false("gen4"%in%QFR$name)
         expect_false("gen3"%in%QFR$name)
         
         
         # checking identified argument as false
         expect_equal(reg2gene(windows = windows,
                               interactions = interactions,
                               geneAnnotations = geneAnnotations[1],
                               identified = FALSE), windows[3:4])
   
   # 3) intset full, geneA full #expO: gen1,gen2,gen3
       QFR <-  reg2gene(windows,
                          interactions = interactions,
                          geneAnnotations = geneAnnotations)
       
           expect_equal(length(QFR),3) # 3 ranges should be identified
           expect_equal(QFR$name,c("gen1","gen3","gen2"))
           # check annotatedAs argument
           expect_equal(QFR$annotatedAs,c("promoter","nearestGene","enhancer"))
           expect_false("gen4"%in%QFR$name)
           
           # checking identified argument as false
           expect_equal(reg2gene(windows = windows,
                                 interactions = interactions,
                                 geneAnnotations = geneAnnotations,
                                 identified = FALSE), windows[4])
   
   # 4) intset full, geneA 1&3 #expO: gen1,gen2,gen3
           QFR <-  reg2gene(windows,
                            interactions = interactions,
                            geneAnnotations = geneAnnotations[c(1,3)])
           
           expect_equal(length(QFR),3) # 3 ranges should be identified
           expect_equal(QFR$name,c("gen1","gen3","gen2"))
           # check annotatedAs argument
           expect_equal(QFR$annotatedAs,c("promoter","nearestGene","enhancer"))
           expect_false("gen4"%in%QFR$name)
           
           # checking identified argument as false
           expect_equal(reg2gene(windows = windows,
                                 interactions = interactions,
                                 geneAnnotations = geneAnnotations,
                                 identified = FALSE), windows[4])
           
   # 5) intset full -2 , geneA full #expO: gen1,gen3
           QFR <-  reg2gene(windows,
                            interactions = interactions[-2],
                            geneAnnotations = geneAnnotations)
           
           # Enhancer missing, but promoter failed as well
           
           expect_equal(length(QFR),3) # 3 ranges should be identified
           expect_equal(QFR$name,c("gen1","gen2","gen3"))
           # check annotatedAs argument
           expect_equal(QFR$annotatedAs,c("promoter","nearestGene","nearestGene"))
           expect_false("gen4"%in%QFR$name)
           
           
           # checking identified argument as false
           expect_equal(reg2gene(windows = windows,
                                 interactions = interactions,
                                 geneAnnotations = geneAnnotations,
                                 identified = FALSE), windows[4])
           
           
           
           # 5)a intset full -2 , geneA -2 #expO: gen1,gen3
           QFR <-  reg2gene(windows,
                            interactions = interactions[-2],
                            geneAnnotations = geneAnnotations[-2])
           
    # 5a Enhancer missing, but promoter failed as well
           
           expect_equal(length(QFR),2) # 3 ranges should be identified
           expect_equal(QFR$name,c("gen1","gen3"))
           # check annotatedAs argument
           expect_equal(QFR$annotatedAs,c("promoter","nearestGene"))
           expect_false("gen4"%in%QFR$name)
           expect_false("gen2"%in%QFR$name)
           
           
           # checking identified argument as false
           expect_equal(reg2gene(windows = windows,
                                 interactions = interactions[-2],
                                 geneAnnotations = geneAnnotations[-2],
                                 identified = FALSE), windows[c(2,4)])
           
           
    # 6) intset 4, geneA 4 - no overlap identified
           QFR <-  reg2gene(windows,
                            interactions = interactions[4],
                            geneAnnotations = geneAnnotations[4])
     
           expect_equal(QFR,"No overlap identified!")
   
        
              # check identified=F argument     
                   
                   expect_equal(reg2gene(windows,
                            interactions = interactions[4],
                            geneAnnotations = geneAnnotations[4],
                            identified = FALSE),windows)
                   
            
         
                   
                  # test distance +/-10000000
                   # increase distance works well
                   
                   QFR <-  reg2gene(windows,
                            interactions = interactions,
                            geneAnnotations = geneAnnotations,
                            distance = 10000000)
                   
                   
                   expect_equal(length(QFR),4) # 3 ranges should be identified
                   expect_equal(QFR$name,c("gen1","gen3","gen4","gen2"))
                   # check annotatedAs argument
                   expect_equal(QFR$annotatedAs,c("promoter","nearestGene",
                                                  "nearestGene","enhancer"))
                   expect_true("gen4"%in%QFR$name)
                   
                        expect_equal(length(reg2gene(windows,
                                   interactions = interactions,
                                   geneAnnotations = geneAnnotations,
                                  distance = 10000000,identified = FALSE)),0)
                   
                        
                        
                      # reduce distance -/+ 100  
                        QFR <-  reg2gene(windows,
                                         interactions = interactions,
                                         geneAnnotations = geneAnnotations,
                                         distance = 100)
                        
                        
                        expect_equal(length(QFR),2) # 3 ranges should be identified
                        expect_equal(QFR$name,c("gen1","gen2"))
                        # check annotatedAs argument
                        expect_equal(QFR$annotatedAs,c("promoter","enhancer"))
                        expect_false("gen4"%in%QFR$name)
                        expect_false("gen3"%in%QFR$name)
                        
                        expect_equal(length(reg2gene(windows,
                                                     interactions = interactions,
                                                     geneAnnotations = geneAnnotations,
                                                     distance = 100,identified = FALSE)),2)    
                 
                   
                       
         
                   
                   # upstream,downstream   [+/-1bp]
                   
                        QFR <-  reg2gene(windows,
                                         interactions = interactions,
                                         geneAnnotations = geneAnnotations,
                                         downstream = 1,
                                         upstream =  1,
                                         distance = 100)
                        
                        
                   expect_equal(length(QFR),2) # 2 ranges should be identified
                   expect_equal(QFR$name,c("gen1","gen2"))
               
                   
                   # but if distance is 1MB,then the 3rd gene should be identified as well
                   expect_equal( reg2gene(windows,
                                          interactions = interactions,
                                          geneAnnotations = geneAnnotations,
                                          downstream = 1, upstream =  100000,
                                          distance = 100)$name,c("gen1","gen2",
                                                                 "gen3"))               
                   
    #################################################################
  # ---------------------------------------------------------------                 
   # 7) intset full, geneA 1
                   mcols(interactions) <- NULL
                   
                   QFR <- reg2gene(windows = windows,
                                   interactions = interactions,
                                   geneAnnotations = geneAnnotations[1],
                                    annotateInteractions = TRUE) 
                   
                   expect_equal(length(QFR),1)
                   expect_equal(QFR$name,"gen1")
                   expect_equal(granges(second(QFR)),granges(geneAnnotations[1]))
                   
                   # checking identified argument as false
                   expect_equal(reg2gene(windows = windows,
                                         interactions = NULL,
                                         geneAnnotations = geneAnnotations[1],
                                         identified = FALSE,
                                         annotateInteractions = TRUE), 
                                windows[-1])
                   
                   
                     # if annotateInteractions needed to be TRUE, but it is not
                      # defined
                   
                     expect_true(reg2gene(windows = windows,
                              interactions = interactions,
                              geneAnnotations = geneAnnotations[1])==
                              "ERROR! Try setting annotateInteractions = TRUE")
                     
                   
                   
                   
                   
                    # 8) intset full, geneA full
                     QFR <- reg2gene(windows = windows,
                                     interactions = interactions,
                                     geneAnnotations = geneAnnotations,
                                     annotateInteractions = TRUE) 
                     
                     expect_equal(length(QFR),3) # 3 ranges should be identified
                     expect_equal(QFR$name,c("gen1","gen2","gen3"))
                     # check annotatedAs argument
                     expect_equal(QFR$annotatedAs,
                                  c("promoter","nearestGene","nearestGene"))
                     expect_false("gen4"%in%QFR$name)
                     
                     
             # if annotateInteractions needed to be TRUE, but it is not
             # defined
                   
                   expect_true(reg2gene(windows = windows,
                                          interactions = interactions,
                                        geneAnnotations = geneAnnotations)==
                             "ERROR! Try setting annotateInteractions = TRUE")
                   
                         
        
                   
                   
                   
                   
           
      
         
 })
         