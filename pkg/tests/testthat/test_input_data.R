
#############################################################
# INPUT DATA UNITESTS
##TESTS to test GRReg1_toy, GRReg1_toy objects and toy.bw input files


##############################################################
# test that GRReg_toy1 and GRReg_toy2 has correct GRanges object form

require(testthat)

 test_that("test that GRReg_toy1 has PREDEFINED GRanges object form", {


    expect_is(GRReg1_toy, 'GRanges',"GRReg1_toy is NOT GRanges")
    expect_equal(length(GRReg1_toy),9)
    expect_equal(colnames(mcols(GRReg1_toy)),
                 c("reg","name2","name","Bench1Exp","Filter2Exp",
                   "Bench2Exp","Filter1Exp"))
    expect_failure(expect_equal(colnames(mcols(GRReg1_toy)),
                 c("name2","name","Bench1Exp","Filter2Exp",
                   "Bench2Exp","Filter1Exp")))

    #Test filtering and benchmarking expected results procedure:
    expect_equal(GRReg1_toy$Filter1Exp,c(1,1,0,1,1,1,1,1,1))
    expect_equal(GRReg1_toy$Filter2Exp,c(1,1,1,1,1,1,1,1,1))
    expect_equal(GRReg1_toy$Bench1Exp,c(2,1,0,0,1,2,3,2,2))
    expect_equal(GRReg1_toy$Bench2Exp,c(1,1,1,1,1,1,1,1,1))

    #Test additional meta-data

    expect_equal(GRReg1_toy$name,paste0("TEST_Reg",1:length(GRReg1_toy)))
    expect_equal(GRReg1_toy$name2,paste0("TEST_Reg",1:length(GRReg1_toy)))

    # test GRanges


    expect_equal(granges(GRReg1_toy),GRanges(c(rep("chr1",8),"chr2"),
            IRanges(c(1,1,15,24,30,5,100,200,1),
                    c(4,4,20,29,31,9,101,201,4))))

    expect_equal(GRReg1_toy$reg , GRanges(c(rep("chr1",8),"chr2"),
                    IRanges(c(5,30,21,5,31,31,102,202,5),
                            c(9,31,25,10,32,32,103,203,9))))

 })


 test_that("test that GRReg_toy2 has PREDEFINED GRanges object form", {


   expect_is(GRReg2_toy, 'GRanges',"GRReg1_toy is NOT GRanges")
   expect_equal(length(GRReg2_toy),14)
   expect_equal(colnames(mcols(GRReg2_toy)),
                c("reg","name2","name"))
   expect_failure(expect_equal(colnames(mcols(GRReg2_toy)),
                               c("name2","name","Bench1Exp","Filter2Exp",
                                 "Bench2Exp","Filter1Exp")))

   #Test additional meta-data

   expect_equal(GRReg2_toy$name,paste0("TEST_Reg",1:length(GRReg2_toy)))
   expect_equal(GRReg2_toy$name2,paste0("TEST_Reg",1:length(GRReg2_toy)))

   # test GRanges
      expect_equal(granges(GRReg2_toy),GRanges(c(rep("chr1",12),rep("chr2",2)),
                            IRanges(c(1,1,5,26,1,31,31,200,200,100,100,102,1,5),
                                c(4,3,9,30,3,33,32,201,201,101,101,103,4,9))))

   expect_equal(GRReg2_toy$reg,GRanges(c(rep("chr1",12),rep("chr2",2)),
                           IRanges(c(1,5,1,31,31,5,5,203,203,102,103,100,5,1),
                             c(3,10,2,32,32,10,10,204,204,103,104,101,9,4))))



 })





################################################################
# test test.bw and test2.bw
# NOT TESTED: due to dependency to rtracklayer, afterwards test will show 
# difference if exists
 
 # test_that("test that test.bw has been correctly inputed", {
 # 
 # test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
 # 
 # 
 #    import.test.bw <- rtracklayer::import(test.bw)
 #    seqlengths (import.test.bw) <- NA
 # 
 #    expect_equal(granges(import.test.bw),
 #                 GRanges(c("chr1","chr1","chr1","chr2","chr2","chr2"),
 #                         IRanges(c(1,5,10,1,5,10),
 #                                 c(4,9,25,4,9,25))))
 #    expect_equal(import.test.bw$score,0:5)
 # 
 # }
 # 
 # 
 # test_that("test that test2.bw has been correctly inputed", {
 # 
 #   test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
 # 
 #   import.test.bw <- rtracklayer::import(test2.bw)
 #   seqlengths (import.test.bw) <- NA
 # 
 #   expect_equal(granges(import.test.bw),
 #                GRanges(c("chr1","chr1","chr1","chr2","chr2","chr2"),
 #                        IRanges(c(1,5,10,1,5,10),
 #                                c(4,9,25,4,9,25))))
 #   expect_equal(import.test.bw$score,seq(0,1,0,0,0,1))
 # 
 # }
 # 




