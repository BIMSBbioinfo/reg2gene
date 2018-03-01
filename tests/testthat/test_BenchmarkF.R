# #######################################
# # -------- ConfusionMatrix & Benchmark tests
#
#
# Dataset description:
# examples created to cover following possibilities for benchmarking and 
# filtering

# 0) anchor1(re2gene) overlaps anchor1 [1 ranges in benchmark dataset]
#     chr1 [  1,   4] ---      chr1 [ 30,  31]
# 1) anchor1(re2gene) overlaps anchor1 [2 ranges in benchmark dataset]
#     chr1 [200, 201]      * | chr1:202-203
# 2) anchor1 overlaps anchor2 [2 ranges in benchmark dataset]
#     chr1 [  5,   9]      * |   chr1:31-32
# 3) anchor1(re2gene) overlaps anchor2 [1 range in benchmark dataset] but
#    as well anchor1 from another GInteraction pair
#    chr1 [  1,   4] ---      chr1 [  5,   9] 
# 4) no overlap with benchmark dataset but NOT filtered out
#     chr1 [ 24,  29] ---      chr1 [  5,  10]
# 5) no overlap with benchmark dataset but filtered out
#     chr1 [ 15,  20] ---      chr1 [ 21,  25]
# 6) another chr example
#     chr2 [  1,   4] ---      chr2 [  5,   9]
# 7) >2 anchor1 overlaps - a1a2,a1a2,a1a1
#     chr1 [100, 101] ---      chr1 [102, 103]
# 8) input anchors overlap
# chr1 [ 30,  31] ---      chr1 [ 31,  32]
 
library(reg2gene)
require(GenomicRanges)
require(InteractionSet)
require(stringr)
require(testthat)


reg2Gene <- GInteractions(GRReg1_toy,GRReg1_toy$reg)
benchData <- GInteractions(GRReg2_toy,GRReg2_toy$reg)
benchDataList <- list(benchData,reg2Gene)

# for CM data
reg2GeneBench <- reg2Gene
reg2GeneBench$PValue <- seq(0, 1,length.out = length(reg2Gene))
##################################################
# -----------------------------------------------#
# Benchmark testing()


# NOTE: did not test ignore.strands function:

test_that("test that benchmark/filter inputs correct GInteractions object
          form", {

  PredefinedMetaData <- c("anchor1.Bench1Exp","anchor1.Filter2Exp",
                          "anchor1.Bench2Exp","anchor1.Filter1Exp")
  
      
    expect_is(reg2Gene,"GInteractions") 
    expect_is(benchData,"GInteractions")
    expect_true(all(PredefinedMetaData%in%colnames(mcols(reg2Gene))),
                "Predifined filter/bench results missing in reg2gene")
    

    })
    
test_that("test that benchmark/filter inputs correct list of GIs 
          object form", {    
    
  expect_is(benchDataList,"list","List: not inputed correctly to  benchmarkAssociations ")
  expect_failure(expect_is(benchDataList,"GInteractions"))
  expect_equal(length(benchDataList),2)
  expect_equal(sapply(benchDataList, length),c(14,9))
  
})
  
test_that("test that benchmark outputs correct GIs object form", {
    
    benchR <-   benchmarkAssociations(reg2Gene,benchData,binary=TRUE)
    
    expect_is(benchR,"GInteractions")
    
    expect_true("Bench"%in%names(mcols(benchR)),"Bench column not returned")
    
    # test functionality:
    # ------------------
    
    # binary=T returns binary 1 and 0
    expect_false(all(benchR$Bench==benchR$anchor1.Bench1Exp),"should be binary
                 benchmark")
    expect_true(all(as.logical(benchR$Bench)==as.logical(
                        benchR$anchor1.Bench1Exp)),"counts to binary benchmark")
    
    # binary=F returns not >1 as well
    benchR <-   benchmarkAssociations (reg2Gene,benchData,binary=F)
        expect_true(all(benchR$Bench==benchR$anchor1.Bench1Exp),"shouldn't be 
                 binary benchmark")
    
    # Checking what happends when anchor1&anchor2 both overlap only one region
          A1A1A2 <-   benchmarkAssociations (reg2Gene[1],
                    benchData[1]) 
       expect_equal(A1A1A2$Bench,0)
  
    # Check filtering procedure   
       benchR <-   benchmarkAssociations(reg2Gene,benchData,
                                         binary=TRUE,
                                         preFilter = TRUE ) 
       expect_true(all(as.logical(benchR$Filter)==as.logical(
         benchR$anchor1.Filter1Exp)),"filter fails")
       
        
       filter2 <- benchmarkAssociations(reg2Gene,reg2Gene,
                             binary=TRUE,
                             preFilter = TRUE ) 
       
       expect_true(all(as.logical(filter2$Filter)==as.logical(
         benchR$anchor1.Filter2Exp)),"filter fails")

       # FORCE BY NAME
       
       reg2Gene$name <- reg2Gene$anchor1.name
       benchData$name <- benchData$anchor1.name
       
      ForcingName <- benchmarkAssociations(reg2Gene,
                             benchData,
                             binary=FALSE,
                             forceByName = T)
       
      expect_equal(ForcingName$Bench,c(1,2,0,2,1,0,3,1,2))
      
     
      
      })
    
test_that("test that benchmark outputs correct GIs object for list of GIs",{
    
    benchRList <-   benchmarkAssociations (reg2Gene,benchDataList,binary=TRUE)
    
    # OUTPUT OK:
    expect_is(benchRList,"GInteractions")
    
    #if bench input list, then bench column N should be equal to the length of 
    # the list
    expect_equal(length(mcols(reg2Gene))+length(benchDataList),
                 length(mcols(benchRList)))
    
    # names of the list added ok
        names(benchDataList) <- c("B1","B2")
        benchRList <-   benchmarkAssociations (reg2Gene,benchDataList,binary=TRUE)
    
    expect_true(all(c("B1","B2")%in%colnames(mcols(benchRList))),
                "list names adjusted correctly in benchmark f()")
    
    # test functionality:
    # ------------------
    expect_equal(benchRList$B1,c(1,1,0,0,1,1,1,1,1))
    expect_equal(benchRList$B2,rep(1,9))
    
    # testing binary argument
    benchRList <-   benchmarkAssociations (reg2Gene,benchDataList,binary=FALSE)
    expect_equal(benchRList$B1,c(2,1,0,0,1,2,3,2,2))
    
    })
 






#######################################
# -------- ConfusionMatrixReg2Gene 
# 
# context("confusionMatrix function runs correctly")


test_that("Input tests for confusionMatrix", {
  
  
  expect_is(reg2GeneBench,"GInteractions")
  
  # given programmed examples, ppv should be returned as 1
  expect_equal(confusionMatrix(reg2GeneBench,
                               thresholdID="PValue",
                               benchCol = "anchor1.Bench1Exp"),1)
  
  expect_error(confusionMatrix(reg2GeneBench))
  
  
})

# TEST benchCol: DONE
test_that("benchCol for confusionMatrix is correctly used", {                     
  
  # works if benchCol  defined  + PPV
    expect_is(confusionMatrix(reg2GeneBench,
                            thresholdID="PValue",
                            benchCol = "anchor1.Bench1Exp"),"numeric")
  
  # works if benchCol  defined  + PPV + filter column added
  expect_is(confusionMatrix(reg2GeneBench,
                            thresholdID="PValue",
                            benchCol = "anchor1.Bench1Exp",
                            prefilterCol="anchor1.Filter1Exp"),"numeric")
  
  # works if benchCol defined + ConfusionMatrix statistics
    expect_is(confusionMatrix(reg2GeneBench,
                            thresholdID="PValue",
                            benchCol = "anchor1.Bench1Exp",
                            statistics = "ConfusionMatrix"),"list")
 
  # NOT working if benchCol NOT defined & Bench column does NOT exist
    expect_error(confusionMatrix(reg2GeneBench,
                              thresholdID="PValue",
                              statistics = "ConfusionMatrix"))
    
    # works if benchCol NOT defined, but Benchmark column exists
    reg2GeneBench$Bench <- reg2GeneBench$anchor1.Bench1Exp

    expect_is(confusionMatrix(reg2GeneBench,
                              thresholdID="PValue",
                              statistics = "ConfusionMatrix"),"list")
    
    
})

# TEST thresholdID argument: DONE
test_that("thresholdID  for confusionMatrix is correctly used", {                     
  
  # if thresholdID  not defined - expect error
  expect_success(expect_error(confusionMatrix(reg2GeneBench, 
                                              benchCol = "anchor1.Bench1Exp")))
  
  # if thresholdID  wrongly defined - expect error
  expect_success(expect_error(confusionMatrix(reg2GeneBench,
                                              thresholdID="PVal", 
                                              benchCol = "anchor1.Bench1Exp")))
  
  # all fine 
  expect_is(confusionMatrix(reg2GeneBench,
                            thresholdID="PValue",
                            benchCol = "anchor1.Bench1Exp"),"numeric")
  
})

# TEST prefilterCol argument: DONE
test_that("prefilterCol  for confusionMatrix is correctly used", {                     
  
  # function works well with both statistics:ConfusionMatrix and filter on
  expect_failure(expect_error(confusionMatrix(reg2GeneBench, 
                                            benchCol = "anchor1.Bench1Exp",
                                            prefilterCol = "anchor1.Filter1Exp",
                                            thresholdID = "PValue",
                                            statistics="ConfusionMatrix")))
  
  # function works well with both statistics:ConfusionMatrix and filter off
  expect_failure(expect_error(confusionMatrix(reg2GeneBench, 
                                              benchCol = "anchor1.Bench1Exp",
                                              thresholdID = "PValue",
                                              statistics="ConfusionMatrix")))
  
  # function should fail if thresholdID is missing
  expect_error(confusionMatrix(reg2GeneBench,benchCol = "anchor1.Bench1Exp",
                              statistics="ConfusionMatrix"))
  
  
  # function works well with statistics:PPV
  expect_failure(expect_error(confusionMatrix(reg2GeneBench, 
                                              benchCol = "anchor1.Bench1Exp",
                                              thresholdID = "PValue",
                                              statistics="PPV")))
  
  
  # function works when prefilterCol defined
  expect_is(confusionMatrix(reg2GeneBench, 
                            benchCol = "anchor1.Bench1Exp",
                            prefilterCol = "anchor1.Filter1Exp",
                            thresholdID="PValue"),"numeric")
  
  # function works when prefilterCol NOT defined
  expect_is(confusionMatrix(reg2GeneBench, 
                            benchCol = "anchor1.Bench1Exp",
                            thresholdID="PValue"),"numeric")
  
  ######################################
  #---------functionality given prefilterCol
  # functionality when prefilterCol==NULL (not defined)
  ConfusionMatrix.F <-   confusionMatrix(reg2GeneBench, 
                                         benchCol = "anchor1.Bench1Exp",
                                         prefilterCol = "anchor1.Filter1Exp",
                                         thresholdID="PValue",
                                         statistics="ConfusionMatrix")
  
  # functionality when prefilterCol!=NULL
  ConfusionMatrix.NF <- confusionMatrix(reg2GeneBench, 
                                        benchCol = "anchor1.Bench1Exp",
                                        thresholdID="PValue",
                                        statistics="ConfusionMatrix")
  
  # filtering gives different N of TN given toy example
  expect_true(ConfusionMatrix.F$TN!=ConfusionMatrix.NF$TN)
  
  
  
})

# thresholdValue and statistics is calculated ok: DONE

test_that("thresholdValue and statistics for confusionMatrix is correct", {

  # function works when no thresholdValue defined - 0.05 default
  
  expect_is(confusionMatrix(reg2GeneBench, 
                            benchCol = "anchor1.Bench1Exp",
                            prefilterCol = "anchor1.Filter1Exp",
                            thresholdID="PValue"),"numeric")
  
  #########################################
  # Correct functionality:
  
  
  ConfusionMatrix.F <- confusionMatrix(reg2GeneBench,
                  thresholdID = "PValue",
                  thresholdValue = 0.05,
                  benchCol = "anchor1.Bench1Exp",
                  prefilterCol = "anchor1.Filter1Exp",
                  statistics = "ConfusionMatrix")
  
  ConfusionMatrix.NF <- confusionMatrix(reg2GeneBench,
                  thresholdID = "PValue",
                  thresholdValue = 0.05,
                  benchCol = "anchor1.Bench1Exp",
                  statistics = "ConfusionMatrix")
  
  # filtering produces different results in confusion matrix
  expect_true(any(unlist(ConfusionMatrix.F)!=unlist(ConfusionMatrix.NF)))
  
  # different thresholdValue produces different results in confusion matrix
  ConfusionMatrix.FT <- confusionMatrix(reg2GeneBench,
                  thresholdID = "PValue", 
                  thresholdValue = 0.5,
                  benchCol = "anchor1.Bench1Exp",
                  prefilterCol = "anchor1.Filter1Exp",
                  statistics = "ConfusionMatrix")
  expect_true(any(unlist(ConfusionMatrix.F)!=unlist(ConfusionMatrix.FT)))
  
  
  #############################
  # functionality is calculated correctly:
  # TESTING TP,FP,TN,FN:
  
  # benchmarking is correct
  ConfusionMatrix.F$NPV <- round(ConfusionMatrix.F$NPV,3)
    expect_true(all(unlist(ConfusionMatrix.F)==c(1,0,1,6,1,0.143,0.25,0.25)))
  
  ConfusionMatrix.NF$Accuracy <- round(ConfusionMatrix.NF$Accuracy,3)
    expect_true(all(unlist(ConfusionMatrix.NF)==c(1,0,2,6,1,0.25,0.333,0.25)))
  
  ConfusionMatrix.FT$F1 <- round(ConfusionMatrix.FT$F1,3)
    expect_true(all(unlist(ConfusionMatrix.FT)==c(3,1,0,4,.75,0,0.375,0.545)))
    

})





