# ---------------------------------------------------------------------------- 
# Unitests for quantification functions
#'
#' # library(testthat)
#' #
#' #
#' ##########################
#' # Create test GenomicRanges
#' 
#' #
#'      #################################
#'      # CREATING EXAMPLES:
#'      
#'      
#'          

library(testthat)
library(GenomicInteractions)
library(stringr)
library(DESeq2)
library(preprocessCore)
library(reg2gene)
# require DESeq2


test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
test3.bw <- system.file("extdata", "test3.bw",package = "reg2gene")


exons <- GRanges(c(rep("chr1",2),"chr2",rep("chr1",3)),
                 IRanges(c(1,7,9,15,1,21),c(4,8,14,20,4,25)),
                 c(rep("+",3),rep("-",3)))
exons$reg <-  exons[c(1,1,3,5,5,5)]
exons$name2 <- exons$name <- paste0("TEST_Reg",c(1,1,3,5,5,5))

exonsGI= GInteractions(exons,exons$reg,
                     name=exons$name,
                     name2=exons$name2)
    


################################################
####### TEST bwToGeneExp

# sampleID - just on example with different names
# normalize  examples below
# check if it works with this or with name as GI anchor1.reg anchor1.name anchor1.name2
# Input exons: INPUT DATA TYPE: all elements present: reg,name,and name2 if GR
# expect geneActSig to be list 
# test that TSS is correctly assigned
# test that N of genes correspond correctly gene name 
# check what happens if no overlap in regions between .bw and exons
# libstand arg works ok - can you come up with example that +/-/* strand produces diff results


 test_that("bwToGeneExp INPUT/OUTPUT is correct in form",{
   
   # TEST INPUT:
   expect_is(exonsGI,"GInteractions")
   expect_is(exons,"GRanges")
   expect_is(test.bw,"character")
   expect_is(test2.bw,"character")
   expect_true(all(str_detect(c(test.bw,test2.bw),"bw")))
   
   expectedMetadata <- c("reg","name","name2")
   
   # test expected input meta-data:
   expect_equal(colnames(mcols(exons)),expectedMetadata) 
   expect_true(all(expectedMetadata[-1]%in%colnames(mcols(exonsGI)))) # test GI
   
   
   # if meta-data column order is different for GRanges it should still work!
   exonsM <- exons
   mcols(exonsM) <- mcols(exonsM)[c(2,3,1)]
   expect_is(bwToGeneExp(exons = exonsM,geneActSignals = test.bw),"GRanges")
   expect_true(length(bwToGeneExp(exons = exons,geneActSignals = test.bw))==
                 length(bwToGeneExp(exons = exonsM,geneActSignals = test.bw)))
   
   # if meta-data order is different for GI it should still work!
   exonsGIM <- exonsGI
   mcols(exonsGIM) <- mcols(exonsGIM)[c(2,3,4,5,1)]
   expect_is(bwToGeneExp(exons = exonsGIM,geneActSignals = test.bw),"GRanges")
   expect_true(length(bwToGeneExp(exons = exonsGIM,geneActSignals = test.bw))==
                 length(bwToGeneExp(exons = exonsGIM,geneActSignals = test.bw)))
   
   
   # if input names incorrect: GI and GRanges:
   # 1st checked if name is missing
   colnames(mcols(exonsM)) <- c("wrong","is","this")
   expect_error(bwToGeneExp(exons = exonsM,geneActSignals = test.bw))
   # 2nd check if reg is missing
   colnames(mcols(exonsM)) <- c("wrong","name","this")
   expect_error(bwToGeneExp(exons = exonsM,geneActSignals = test.bw))
   
   # if reg&name are present then success:
   colnames(mcols(exonsM)) <- c("wrong","name","reg")
   expect_is(bwToGeneExp(exons = exonsM,geneActSignals = test.bw),"GRanges")
   
   # test for GI
   colnames(mcols(exonsGIM)) <- c("wrong","is","this","and","this2")
   expect_error(bwToGeneExp(exons = exonsGIM,geneActSignals = test.bw))
   
   colnames(mcols(exonsGIM)) <- c("wrong","is","this","and","name")
   expect_error(bwToGeneExp(exons = exonsGIM,geneActSignals = test.bw))
   
   # requires 
   colnames(mcols(exonsGIM)) <- c("name2","is","this","and","name")
   expect_is(bwToGeneExp(exons = exonsGIM,geneActSignals = test.bw),"GRanges")
   

   # TEST OUTPUT
   # is GRanges with lenght>0
   expect_is(bwToGeneExp(exons = exons,geneActSignals = test.bw),"GRanges")
   expect_true(length(bwToGeneExp(exons = exons,geneActSignals = test.bw))!=0)
   # and works for >1 .bw file
   expect_true(length(bwToGeneExp(exons = exons,
                                  geneActSignals = c(test.bw,test2.bw)))!=0)
   
   
   
   # TEST OUTPUTs correct cell names:
   
   # adding different sample IDs:
   expect_equal(names(mcols(bwToGeneExp(exons = exons,
                              geneActSignals = c(test.bw,test2.bw),
                              sampleIDs=c("CellType1","CellType2")))[-c(1:2)]),
                  c("CellType1","CellType2"))
   
   
   # adding different sample IDs - swap
   expect_equal(names(mcols(bwToGeneExp(exons = exons,
                              geneActSignals = c(test.bw,test2.bw),
                              sampleIDs=c("CellType2","CellType1")))[-c(1:2)]),
                c("CellType2","CellType1"))
   
   # test default names - basenames expected
    expect_false(all(names(mcols(bwToGeneExp(exons = exons,
                    geneActSignals = c(test.bw,test2.bw)))[-c(1:2)])==
                c("CellType2","CellType1")))
   
   })
 
 
 test_that("bwToGeneExp performs correctly in assoc genes with exons",{
 
   # TEST OUTPUT IF NO OVERLAP BETWEEN
   # check what happens if no overlap in regions between .bw and exons
   # test that TSS is correctly assigned
   
   
   # test that N of genes correspond correctly gene name: for GI and GRanges
   expect_equal(sort(bwToGeneExp(exons = exons,geneActSignals = test.bw)$name),
                unique(exons$name))
   expect_equal(sort(bwToGeneExp(exons = exonsGI,geneActSignals = test.bw)$name),
                unique(exonsGI$name)) 
   # and for >1 .bw file
    expect_equal(sort(bwToGeneExp(exons = exons,
            geneActSignals = c(test.bw,test2.bw))$name),unique(exons$name))
   # by hand  
    expect_equal(sort(bwToGeneExp(exons = exons,geneActSignals = test.bw)$name),
                 c("TEST_Reg1","TEST_Reg3","TEST_Reg5"))
    
    exonsM <- exons
    exonsM$name <- c("this","this","is","wrong","wrong","wrong")
   # names changed: OK!   
    expect_error(expect_equal(bwToGeneExp(exons = exonsM,
                                          geneActSignals = test.bw)$name,
                 c("this","is","wrong")))
    
    expect_equal(bwToGeneExp(exons = exonsM,
                                          geneActSignals = test.bw)$name,
                              c("wrong","this","is"))
    
    
    # names changed 2 > more genes obtained
    exonsM$name <- c("this","is","wrong")
    expect_equal(bwToGeneExp(exons = exonsM,geneActSignals = test.bw)$name,
                 c("this","is","wrong","this","is","wrong"))

   
    # checking that TSS is correctly assigned:
    
    # test that TSS of a gene is correctly reported
   expect_equal(granges(bwToGeneExp(exons = exonsGI,
                                    geneActSignals = test.bw))[c(2,3,1)],
    promoters(unique(second(exonsGI)),1,1))
   
   expect_equal(granges(bwToGeneExp(exons = exonsGI,geneActSignals = test.bw)),
                GRanges(c("chr1","chr1","chr2"),
                        IRanges(c(4,0,8),c(5,1,9)),
                        c("-","+","+")))
    
   # if genes defined differently: as before
   expect_equal(granges(bwToGeneExp(exons = exonsM,geneActSignals = test.bw)),
                GRanges(c("chr1","chr1","chr1","chr1","chr1","chr2"),
                        IRanges(c(4,4,4,0,0,8),c(5,5,5,1,1,9)),
                        c("-","-","-","+","+","+")))
   
   
   # check what happens if no overlap in regions between .bw and exons
   expect_error(bwToGeneExp(exons =  shift(exonsM,100),geneActSignals = test.bw))
   expect_error(bwToGeneExp(exons =  shift(exonsM,100),
                            geneActSignals = c(test.bw,test2.bw)))
   # even if it fails for one .bw it fails for all others as well:
  expect_error(bwToGeneExp(exons =  shift(exonsM,100),
                           geneActSignals = c(test3.bw,test2.bw)))
   # success for shifted .bw file
  expect_is(bwToGeneExp(exons =  shift(exonsM,100),geneActSignals = test3.bw),
            "GRanges")
   
   })
   
  
 test_that("bwToGeneExp performs correctly ",{

  TestExp <- c(1.4667,0.3333,4.8333)
  Test2Exp <- c(0,0.3333,0.8333)
   
  # test2 quantified correctly
    expect_equal(round(
      bwToGeneExp(exons = exonsGI,geneActSignals = test2.bw)$test2,4),Test2Exp)
    
  
  # test quantified correctly
    expect_equal(round(bwToGeneExp(exons = exonsGI,
                                   geneActSignals = test.bw)$test,4),TestExp) 
    
    
  # if strandness is changed: then genes on - strand and + strand get quantified
  # based on different libraries, meaning that result will be 1 region is test
  # and test2.bw used
    # tested example: +/-, -/+, +/-/*,-/+/*, */-/+, */+/-, +/-/*/-/+
    # exons starting with + and -
    
    # libraries: +/-
    expect_equal(round(bwToGeneExp(exons = exonsGI,
                                   geneActSignals =  c(test.bw,test2.bw),
                                   libStrand = c("+","-"))$test,4),
                    c(Test2Exp[1:2],TestExp[3]))
    
    # libraries: -/+
    
    expect_equal(round(bwToGeneExp(exons = exonsGI,
                                   geneActSignals =  c(test.bw,test2.bw),
                                   libStrand = c("-","+"))$test,4),
                 c(TestExp[1:2],Test2Exp[3]))
    
    # + library needs to be followed by - library
   expect_error(bwToGeneExp(exons = exonsGI,
               geneActSignals = c(test.bw,test2.bw),
               libStrand = c("+","+")))
   
   # -/+/*
   # works correctly with combination of stranded and unstranded libraries
   expect_is(bwToGeneExp(exons = exonsGI,
               geneActSignals = c(test.bw,test2.bw,test.bw),
               libStrand = c("-","+","*")),"GRanges")
   
   
   # -/+/* example
   
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test2.bw,test.bw),
                                  libStrand = c("-","+","*"))$test,4),TestExp)
   
   
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test2.bw,test.bw),
                                  libStrand = c("-","+","*"))$test2,4),
                c(TestExp[1:2],Test2Exp[3]))

 
   # +/-/* example  
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test2.bw,test2.bw),
                                  libStrand = c("+","-","*"))$test2,4),Test2Exp)
   
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test2.bw,test2.bw),
                                  libStrand = c("+","-","*"))$test,4),
                c(Test2Exp[1],TestExp[2:3]))
 
   
   # */+/- example  
   
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test2.bw,test.bw),
                                  libStrand = c("*","+","-"))$test,4),TestExp)
   
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test2.bw,test.bw),
                                  libStrand = c("*","+","-"))$test2,4),
                c(TestExp[1],Test2Exp[2:3]))
   
   
   
   # */-/+ example  
   
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test.bw,test2.bw),
                                  libStrand = c("*","-","+"))$test,4),TestExp)
   
   expect_equal(round(bwToGeneExp(exons = exonsGI,
                                  geneActSignals =  c(test.bw,test.bw,test2.bw),
                                  libStrand = c("*","-","+"))$test2,4),
                c(TestExp[1],Test2Exp[2:3]))
   
   
   # +/-/*/-/+ example  
   
   tmp <- bwToGeneExp(exons = exonsGI,
               geneActSignals =  c(test.bw,test2.bw,
                                   test.bw,test.bw,test2.bw),
               sampleIDs = c("name1","name1","name3","name4","name4"),
               libStrand = c("+","-","*","-","+"))

   
   expect_equal(round(c(tmp$name1,tmp$name3,tmp$name4),4),
   c(Test2Exp[1],TestExp[2:3],TestExp,TestExp[1],Test2Exp[2:3]))
   
   
   #test what happens when names differ
   
       tmp2 <- bwToGeneExp(exons = exonsGI,
                          geneActSignals =  c(test.bw,test2.bw,
                                              test.bw,test.bw,test2.bw),
                          sampleIDs = c("name1","name3","name1","name4","name4"),
                          libStrand = c("+","-","*","-","+"))
       
       expect_equal(round(c(tmp$name1,tmp$name3,tmp$name4),4),
                    round(c(tmp2$name1,tmp2$name1.1,tmp2$name4),4))
   
       
   # test exon orientation
         expect_equal (bwToGeneExp(exons = exonsGI,
                     geneActSignals =  c(test.bw,test2.bw),
                     libStrand = c("+","-"))$test,
                    bwToGeneExp(exons = c(exonsGI[4:6],exonsGI[1:3]),
                     geneActSignals =  c(test.bw,test2.bw),
                     libStrand = c("+","-"))$test)
         
         
         expect_equal (bwToGeneExp(exons = exonsGI,
                                   geneActSignals =  c(test2.bw,test.bw),
                                   libStrand = c("+","-"))$test,
                       bwToGeneExp(exons = c(exonsGI[4:6],exonsGI[1:3]),
                                   geneActSignals =  c(test2.bw,test.bw),
                                   libStrand = c("+","-"))$test)
   
   
   # works correctly with combination of stranded and unstranded libraries if 
   # stranded are on behind other
   expect_error(bwToGeneExp(exons = exonsGI,
                         geneActSignals = c(test.bw,test2.bw,test.bw),
                         libStrand = c("-","*","-")))
   
  
 })
 
 
 test_that("bwToGeneExp performs correctly given normalization",{
    
   # expected quntifies results
   TestExp <- c(1.4667,0.3333,4.8333)
   Test2Exp <- c(0,0.3333,0.8333)
   
   TestNorm <- cbind(TestExp,Test2Exp)
   
   # recalculating "quantile" normalization
         QuantileExp <- DataFrame(preprocessCore::normalize.quantiles(TestNorm))
    # recalculating "ratio" normalization
         sizeFactors <- DESeq2::estimateSizeFactorsForMatrix(TestNorm)
         RatioExp <- DataFrame(TestNorm*rep(sizeFactors,each=nrow(TestNorm)))
        
          TestNorm <- DataFrame(TestNorm)
    names(RatioExp) <- names(TestNorm) <-  names(QuantileExp) <- c("test",
                                                                   "test2")
   
  #################################################
    ### TESTING:
    
    # function to test normalization procedure:
          TestNormF <- function(y,norm=NULL){
          
            Res <- bwToGeneExp(exons = exonsGI,
                        geneActSignals = c(test.bw,test2.bw),
                        normalize = norm)
            
            TestR <- apply(mcols(Res)[-c(1:2)],1,round,2)
            
           expect_equal(TestR,apply(y,1,round,2))
          
            
          }
          
    # testing normalization procedure done "by hand" and reported from f()      
          TestNormF(y=TestNorm,norm=NULL)
          TestNormF(y=QuantileExp,norm="quantile")
          TestNormF(y=RatioExp,norm="ratio")
          
    
    # if .bigwig files swapped:
        TestNorm <- cbind(Test2Exp,TestExp)
        QuantileExp <- DataFrame(preprocessCore::normalize.quantiles(TestNorm))
        names(QuantileExp) <- c("test","test2")

        Res <- bwToGeneExp(exons = exonsGI,geneActSignals = c(test2.bw,test.bw),
                           normalize = "quantile")

         expect_true(all(apply(mcols(Res)[-c(1:2)],1,round,2)==
                           apply(QuantileExp,1,round,2)))
 
 })       
         
         

###############################################
####### TEST regActivity          
 
 regRegions <- GRanges(c(rep("chr1",4),rep("chr2",2)),
                       IRanges(c(1,7,9,15,1,15),c(4,8,14,20,4,20)),
                       c(rep("+",3),rep("-",3)))
 regRegions$reg <-  regRegions[c(1,1,3:6)]
 regRegions$name2 <- regRegions$name <- paste0("TEST_Reg",
                                               c(1,1,3:length(regRegions)))
 

 
 # check INPUT/OUTPUT + sampleIDs
 test_that("regActivity INPUT/OUTPUT is correct in form",{
   
   # TEST INPUT:
   expect_is(regRegions,"GRanges")
   expect_is(test.bw,"character")
   expect_is(test2.bw,"character")
   expect_true(all(str_detect(c(test.bw,test2.bw),"bw")))
   
   expectedMetadata <- c("reg","name","name2")
   
   # test expected input meta-data:
   expect_equal(colnames(mcols(regRegions)),expectedMetadata) 
  
   
   
   # TEST OUTPUTs
   
   expect_is(regActivity(regRegions,c(test.bw,test2.bw)),"GRanges")  
   # if meta-data column order is different for GRanges it should still work!
   regRegionsM <- regRegions
   mcols(regRegionsM) <- mcols(regRegionsM)[c(2,3,1)]
   expect_is(regActivity(regRegionsM,c(test.bw,test2.bw)),"GRanges")  
   
   expect_true(length(regActivity(regRegions,c(test.bw,test2.bw)))==
                 length(regActivity(regRegionsM,c(test.bw,test2.bw))))
   
  
   # if input names incorrect: success because gene name not needed, just plain
   # GRanges requested
   colnames(mcols(regRegionsM)) <- c("wrong","is","this")
   expect_is(regActivity(regRegionsM,c(test.bw,test2.bw)),"GRanges")  
   
   
   # TEST OUTPUT
   # is GRanges with lenght>0
   expect_is(regActivity(regRegions, test.bw),"GRanges")
   expect_true(length(regActivity(regRegions, test.bw))!=0)
   # and works for >1 .bw file
   expect_true(length(regActivity(regRegions,c(test.bw,test2.bw)))!=0)
   
   # TEST OUTPUTs correct cell names:
   
   # adding different sample IDs:
   expect_equal(names(mcols(regActivity(regRegions, c(test.bw,test2.bw),
                            sampleIDs=c("CellType1","CellType2")))),
                c("CellType1","CellType2"))
   
   
   # adding different sample IDs - swap
   expect_equal(names(mcols(regActivity(regRegions, c(test.bw,test2.bw),
                                        sampleIDs=c("CellType2","CellType1")))),
                c("CellType2","CellType1"))
   
   # test default names - basenames expected
   expect_false(all(names(mcols(regActivity(regRegions, c(test.bw,test2.bw))))==
                      c("CellType2","CellType1")))
   
 })
 
 # check normalize
 test_that("regActivity performs correctly given normalization",{
   
   # expected quntifies results
   TestExp <- c(0.000000,1.000000,1.833333,2.000000,3.000000,5.000000)
   Test2Exp <- c(0.0000000,1.0000000,0.1666667,0.0000000,0.0000000,1.0000000)
   
   TestNorm <- cbind(TestExp,Test2Exp)
   
   # recalculating "quantile" normalization
   QuantileExp <- DataFrame(preprocessCore::normalize.quantiles(TestNorm))
   # recalculating "ratio" normalization
   sizeFactors <- DESeq2::estimateSizeFactorsForMatrix(TestNorm)
   RatioExp <- DataFrame(TestNorm*rep(sizeFactors,each=nrow(TestNorm)))
   
   TestNorm <- DataFrame(TestNorm)
   names(RatioExp) <- names(TestNorm) <-  names(QuantileExp) <- c("test",
                                                                  "test2")
   
   #################################################
   ### TESTING:
 
   # function to test normalization procedure:
   TestNormF <- function(y,norm=NULL){
     
     Res <- regActivity(regRegions,
                        c(test.bw,test2.bw),
                        normalize = norm)
     
     TestR <- apply(mcols(Res),1,round,2)
     
     expect_equal(TestR,apply(y,1,round,2))
     
     
   }
   
   # testing normalization procedure done "by hand" and reported from f()      
   TestNormF(y=TestNorm,norm=NULL)
   TestNormF(y=QuantileExp,norm="quantile")
   TestNormF(y=RatioExp,norm="ratio")
   
   
   # if .bigwig files swapped:
   TestNorm <- cbind(Test2Exp,TestExp)
   QuantileExp <- DataFrame(preprocessCore::normalize.quantiles(TestNorm))
   names(QuantileExp) <- c("test","test2")
   
   Res <- regActivity(regRegions,c(test2.bw,test.bw),
                      normalize = "quantile")
   
   expect_true(all(apply(mcols(Res),1,round,2)==
                     apply(QuantileExp,1,round,2)))
   
 })     
 
 
 # check performance

 test_that("regActivity performs correctly ",{
   
   TestExp <- c(0.000000,1.000000,1.833333,2.000000,3.000000,5.000000)
   Test2Exp <- c(0.0000000,1.0000000,0.1666667,0.0000000,0.0000000,1.0000000)
   
   # test2 quantified correctly
   expect_equal(round(regActivity(regRegions,test2.bw)$test2,4),
                round(Test2Exp,4))
   
   
   # test quantified correctly
   expect_equal(round(regActivity(regRegions,test.bw)$test,4),
                round(TestExp,4))
   
   # MISSING check isCovNA
   
   #regActivity(regRegions,test.bw,isCovNA = T)

   })

 
 ###############################################
 ####### TEST regActivity   
 
 # Creating Datasets
 EnhActivity <- GRReg1_toy
 mcols(EnhActivity) <- NULL
 EnhActivity$bw1 <- rep(1,length(GRReg1_toy))
 EnhActivity$bw2 <- rep(2,length(GRReg1_toy))
 EnhActivity$bw3 <- rep(3,length(GRReg1_toy))
 
 GeneExpression <- GRReg2_toy
     mcols(GeneExpression) <- mcols(GeneExpression)[c(1,3,2)] 
     GeneExpression$bw1 <- rep(3,length(GeneExpression))
     GeneExpression$bw2 <- rep(4,length(GeneExpression))
     GeneExpression$bw4 <- rep(5,length(GeneExpression))
     

 test_that("regActivityAroundTSS performs correctly ",{
   
   # TEST INPUT:
     expect_is(EnhActivity,"GRanges")
     expect_is(GeneExpression,"GRanges")
     
     expect_is(test.bw,"character")
     expect_is(test2.bw,"character")
     expect_true(all(str_detect(c(test.bw,test2.bw),"bw")))  
             
     expectedMetadataW <- c("reg","name","name2","bw1","bw2","bw4")
     
     # test expected input meta-data:
     expect_equal(colnames(mcols(GeneExpression)),expectedMetadataW) 
     expect_equal(colnames(mcols(EnhActivity)),c("bw1","bw2","bw3")) 
  
     
    
     # TEST OUTPUTs
      # if regActivity and tss arguments are switched - error is reported
     expect_error(regActivityAroundTSS(tss = EnhActivity,
                                       regActivity=GeneExpression,
                                       upstream=1,downstream=1))
     
     # TEST OUTPUT:
     # Successful example:
     regActivityAroundTSSR <- regActivityAroundTSS(regActivity = EnhActivity,
                                      tss=GeneExpression,
                                      upstream=1,downstream=1)
     
     ########################
     expect_is(regActivityAroundTSSR,"GRangesList")
    
     # names of a list correspond to names in GeneExpression object:
     expect_equal(length(regActivityAroundTSSR),14)
     expect_equal(length(regActivityAroundTSSR),
                  length(unique(GeneExpression$name)))
     expect_true(all(names(regActivityAroundTSSR)%in%
                    unique(GeneExpression$name)))
     
     
     
     
     
     # TEST 1 member of the list: 
       expect_is(regActivityAroundTSSR$TEST_Reg1,"GRanges")
       expect_equal(length(regActivityAroundTSSR$TEST_Reg1),3)
       
       # description of object is stored in:"featureType","name","name2"
       
       expect_equal(colnames(mcols(regActivityAroundTSSR$TEST_Reg1))[1:3],
                    c("featureType","name","name2"))
       
       # quantified bw should be: "bw1","bw2" because GeneExpression is missing 
       # "bw3"
       expect_equal(colnames(mcols(regActivityAroundTSSR$TEST_Reg1))[-(1:3)],
                    c("bw1","bw2"))
       
       # non-overlapping bw filtered out if present in both enhActivity object
       # and GeneExpression object
       expect_false("bw3"%in%colnames(mcols(regActivityAroundTSSR$TEST_Reg1)))
       expect_false("bw4"%in%colnames(mcols(regActivityAroundTSSR$TEST_Reg1)))
 
       # gene & regulatory need to be elements of featureType 
       expect_true(all(c("gene","regulatory")%in%
                         regActivityAroundTSSR$TEST_Reg1$featureType))
      # if >2 elements then: 1 gene & other are regulatory
       expect_equal(regActivityAroundTSSR$TEST_Reg1$featureType,c("gene",
                                                    "regulatory","regulatory"))
       
       expect_equal(regActivityAroundTSSR$TEST_Reg1$name,c("TEST_Reg1",
                                                        "chr1:1-4","chr1:1-4"))
       
       # Test results of quantifying genes and enh activity
       expect_equal(regActivityAroundTSSR$TEST_Reg1$bw1,c(3,1,1))
       expect_equal(regActivityAroundTSSR$TEST_Reg1$bw2,c(4,2,2))
       
       
       
       
       
       
    ####################################3  
    # TEST another member of the list: 
       expect_is(regActivityAroundTSSR$TEST_Reg11,"GRanges")
       expect_equal(length(regActivityAroundTSSR$TEST_Reg11),2)
       
       # description of object is stored in:"featureType","name","name2"
       expect_equal(colnames(mcols(regActivityAroundTSSR$TEST_Reg11))[1:3],
                    c("featureType","name","name2"))
     
       # gene & regulatory need to be elements of featureType 
       expect_true(all(c("gene","regulatory")%in%
                         regActivityAroundTSSR$TEST_Reg11$featureType))
       # if >2 elements then: 1 gene & other are regulatory
       expect_equal(regActivityAroundTSSR$TEST_Reg11$featureType,c("gene",
                                                                  "regulatory"))
       
         
       # quantified bw should be: "bw1","bw2" because GeneExpression is missing 
       # "bw3"
       expect_equal(colnames(mcols(regActivityAroundTSSR$TEST_Reg11))[-(1:3)],
                    c("bw1","bw2"))
       
       expect_equal(regActivityAroundTSSR$TEST_Reg11$name,c("TEST_Reg11",
                                                            "chr1:100-101"))
           
       # Test results of quantifying genes and enh activity
       expect_equal(regActivityAroundTSSR$TEST_Reg11$bw1,c(3,1))
       expect_equal(regActivityAroundTSSR$TEST_Reg11$bw2,c(4,2))
       
       
 }) 
 
 
 test_that("regActivity extends TSS correctly ",{
   
   
   regActivityAroundTSSR <- regActivityAroundTSS(EnhActivity,
                                                 GeneExpression,
                                                 upstream=1,
                                                 downstream=1)
   
   regActivityAroundTSSR5 <- regActivityAroundTSS(EnhActivity,
                                                 GeneExpression,
                                                 upstream=5,
                                                 downstream=5)
   
   expect_is(regActivityAroundTSSR5,"GRangesList")
   
   # names of a list correspond to names in GeneExpression object:
   expect_equal(length(regActivityAroundTSSR5),14)
   expect_equal(length(regActivityAroundTSSR5),
                length(unique(GeneExpression$name)))
   expect_true(all(names(regActivityAroundTSSR5)%in%
                     unique(GeneExpression$name)))
   
   
   # TEST 1 member of the list: 
       expect_is(regActivityAroundTSSR5$TEST_Reg1,"GRanges")
       
       # description of object is stored in:"featureType","name","name2"
       
       expect_equal(colnames(mcols(regActivityAroundTSSR5$TEST_Reg1))[1:3],
                    c("featureType","name","name2"))
       
       # quantified bw should be: "bw1","bw2" because GeneExpression is missing 
       # "bw3"
       expect_equal(colnames(mcols(regActivityAroundTSSR5$TEST_Reg1))[-(1:3)],
                    c("bw1","bw2"))
       
       # non-overlapping bw filtered out if present in both enhActivity object
       # and GeneExpression object
       expect_false("bw3"%in%colnames(mcols(regActivityAroundTSSR5$TEST_Reg1)))
       expect_false("bw4"%in%colnames(mcols(regActivityAroundTSSR5$TEST_Reg1)))
       
       # gene & regulatory need to be elements of featureType 
       expect_true(all(c("gene","regulatory")%in%
                         regActivityAroundTSSR5$TEST_Reg1$featureType))
       
       
       # A PART THAT HAS CHANGED:
       # 1. length has changed:
       expect_equal(length(regActivityAroundTSSR5$TEST_Reg1),4)
       # if >2 elements then: 1 gene & other are regulatory
       expect_equal(regActivityAroundTSSR5$TEST_Reg1$featureType,c("gene",
                                                                   "regulatory",
                                                                   "regulatory",
                                                                  "regulatory"))
       
       expect_equal(regActivityAroundTSSR5$TEST_Reg1$name,c("TEST_Reg1",
                                                            "chr1:1-4",
                                                            "chr1:1-4",
                                                            "chr1:5-9"))
       
       # Test results of quantifying genes and enh activity
       expect_equal(regActivityAroundTSSR5$TEST_Reg1$bw1,c(3,1,1,1))
       expect_equal(regActivityAroundTSSR5$TEST_Reg1$bw2,c(4,2,2,2))
       
    # HOWEVER, +/-5 did not bring any changes to  TEST_Reg11
       expect_equal(length(regActivityAroundTSSR5$TEST_Reg11),
                    length(regActivityAroundTSSR5$TEST_Reg11))
})
   
   