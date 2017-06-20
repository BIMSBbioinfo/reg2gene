# library(reg2gene)
# 
# #######################################
# # -------- ConfusionMatrixReg2Gene
# 
# 
# context("ConfusionMatrixReg2Gene function runs correctly")
# 
# 
# 
# test_that("ConfusionMatrixReg2Gene arguments are defined correctly", {
#   expect_is(BenchMarkedReg2Gene, 'GRanges')
#   expect_false(is.null(thresholdID),"Column to threshold is not defined")
#   expect_match(colnames(GenomicRanges::mcols(BenchMarkedReg2Gene)),thresholdID,all=F,"Column for thresholding not available")
#   expect_match(colnames(GenomicRanges::mcols(BenchMarkedReg2Gene)),"BenchmarkO",all=F,
#                info = "Benchmark0 column missing, or benchmark column should be named in Benchmark0")
#   expect_is((GenomicRanges::mcols(BenchMarkedReg2Gene)[[thresholdID]]), 'numeric')
#   expect_is(GenomicRanges::mcols(BenchMarkedReg2Gene)$BenchmarkO,"logical")})
# 
# test_that("ConfusionMatrixReg2Gene reports CORRECT results", {
#   expect_equal(ConfusionMatrixReg2Gene(BenchMarkedReg2Gene_toy,thresholdID = "PValue",thresholdValue = 0.5),0.8)
#   expect_equal(GenomicRanges::mcols(BenchMarkedReg2Gene_toy)$BenchmarkO,c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE))
#   expect_equal(GenomicRanges::mcols(BenchMarkedReg2Gene_toy)$PValue,c(0.05,0.06,0.01,0.5,0,1))
# })
#   
# 
# 
# 
# 
# ############################################
# # TEST BenchMarkReg2Gene()
# 
# context(" BenchMarkReg2Gene() functions runs correct")
# 
# test_that("BenchMarkReg2Gene()  correctly performs counting of Benchmared overlaps", {
#   expect_equal(GenomicRanges::mcols(BenchMarkReg2Gene(GRReg2_toy,GRReg1_toy))$BenchmarkO,
#                c(2,2,0,1,1,0,0,1,1,0,0,0,0,0))
#   expect_equal(GenomicRanges::mcols(BenchMarkReg2Gene(GRReg2_toy,GRReg1_toy,binary=T))$BenchmarkO,
#                c(TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))})
# 
# test_that("OUTPUT object is equal to the input object + 1 addtional meta-data: BenchmarkO", {
#           test <- GRReg1_toy; test$BenchmarkO <- c(2,2,2,0,2,0)
#   expect_equal(BenchMarkReg2Gene(GRReg1_toy,GRReg2_toy),test)
#           test$BenchmarkO <- c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE)
#   expect_equal(BenchMarkReg2Gene(GRReg1_toy,GRReg2_toy,binary = T),test)
# 
#   })
# 
# 
# 
# 
# ############################################
# # TEST OverlapRegions()
# 
# context("OverlapRegions() functions runs correct")
# 
# 
# test_that("input: Reg1 is a GRanges object containing correct meta-data", {
#   expect_is(object = Reg1, 'GRanges')
#   #expect_gte(ncol(GenomicRanges::mcols(Benchmark)), 3)
#   expect_match(colnames(GenomicRanges::mcols(Reg1)),"reg",all=F)
#   expect_is(Reg1$reg, 'GRanges')
# })
# 
# test_that("input: Reg2  is a GRanges object containing correct meta-data", {
#   expect_is(object = Reg2, 'GRanges')
#   #expect_gte(ncol(GenomicRanges::mcols(Reg1)), 3)
#   expect_match(colnames(GenomicRanges::mcols(Reg2)),"reg",all=F)
#   expect_is(Reg2$reg, 'GRanges')
# })
# 
# test_that("Reg1 or Reg2 is overlapped correctly",{
#   expect_match(colnames(GenomicRanges::mcols(OverlapRegions(GRReg1_toy,GRReg2_toy))),"Coor1Coord2PAIR",
#                all=F,"Coor1Coord2PAIR created or not")
#   expect_equal(GenomicRanges::mcols(OverlapRegions(GRReg1_toy,GRReg2_toy))$"Coor1Coord2PAIR",
#                c(1,2,3,5)) # enhancer-gene pairs reported correctly
#   expect_equal(GenomicRanges::mcols(OverlapRegions(GRReg2_toy,GRReg1_toy))$"Coor1Coord2PAIR",
#                c(1,1,4,8)) # enhancer-gene pairs reported correctly
#       # genomic ranges selected correclty
#   expect_equal(GRReg1_toy$reg[c(1,2,3,5)],GenomicRanges::mcols(OverlapRegions(GRReg1_toy,GRReg2_toy))$"Reg1_reg")
#   expect_equal(GRReg2_toy$reg[c(1,1,4,8)],GenomicRanges::mcols(OverlapRegions(GRReg2_toy,GRReg1_toy))$"Reg1_reg")
#   })
# 
# 
# 
# 
# ############################################
# # TEST ComplexOverlaps()
# 
# context("ComplexOverlaps() functions runs correct")
# 
# test_that("switch function runs correctly", {expect_equal(GRReg2_toy,.switchReg(.switchReg(GRReg2_toy)))})
# 
# test_that("Reg1 or Reg2 is overlapped correctly",{
#   expect_equal(GenomicRanges::mcols(ComplexOverlaps(GRReg2_toy,GRReg1_toy))$"Coor1Coord2PAIR",
#                c(1,1,4,8,2,2,5,9))
#   # test that all ComplexOverlapped regions come from the Coord1 of GRReg2_toy, no mismatches with Coord2 region
#   expect_equal(unique(data.frame(findOverlaps(ComplexOverlaps(GRReg2_toy,GRReg1_toy),GRReg2_toy))$queryHits),
#                1:length(ComplexOverlaps(GRReg2_toy,GRReg1_toy)))
# })
#   
#   
#   