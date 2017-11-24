#' # # ---------------------------------------------------------------------------- #
#' # # test for associateReg2Gene function
#' #
#' # library(testthat)
#' #
#' #
#'  ##########################
#'  # Create test GenomicRanges
#' #
#' #
#'      #################################
#'      # DESCRIPTION:
#' #
#' #
#'      # 1. gr - example 1 gene + 2 enh
#' #
#'      # all others - 1gene - 1 enhancer
#' #
#'      # 2. gr0 - example 1 gene + 1 enh
#'      # 3. gr0a - example 1 gene + 1 enh (different enh)
#'      # 4. gr0.failure -     example 1 gene + 1 enh (gene name repeated,
#'      #                    failure exp)
#'      # 6. gr1 - example 1 gene + 1 enh, gene contains 1 NA
#'      # 7. gr2 - example 1 gene + 1 enh, gene&enh contain 1 NA in the same cell
#'      # 8. gr3 - example 1 gene + 1 enh, gene&enh contain 1 NA in the diff cell
#'      # 9. gr4 - example 1 gene + 1 enh, enh contain 1 NA
#' #
#'  gr <- GRanges(seqnames=rep("chr1",3),IRanges(1:3,3:5))
#' #
#'      x <- 1:5
#'      y <- 2:6
#'      z <- 10:14
#' #
#'      GeneInfo <- as.data.frame(matrix(c(rep("gene",3),rep("regulatory",6)),
#'                                       ncol = 3,byrow = T),stringsAsFactors=FALSE)
#'                  colnames(GeneInfo) <- c("featureType","name","name2")
#' #
#'    mcols(gr) <- DataFrame(cbind(GeneInfo,rbind(x,y,z)))
#' #
#' #
#' #
#'  ###########################
#'  # Create different problems in genes and enhancers
#' #
#'        # one gene - one enh example
#'            gr0 <- gr[1:2]
#'        # one gene - one enh example -other - good example for failure due to
#'        # gene - enh names
#'            gr0.failure <- gr[2:3]
#'        # another one gene - one enh example
#'            gr0a <- gr[c(1,3)]
#' #
#'    # Adding missing data
#'    gr4 <- gr3 <-  gr2 <- gr1 <- gr0
#'          # 1  gene contains NA
#'            mcols(gr1)[1,4] <- NA
#'        # 2 gene and enh contain NA in the same cell
#'            mcols(gr2)[,4] <- rep(NA,2) #
#'        # 3 gene and enh contain NA in the different cell
#'            mcols(gr3)[1,4] <- NA
#'            mcols(gr3)[2,5] <- NA
#'        # 4 enh contain NA
#'            mcols(gr4)[2,5] <- NA
#' #
#' #
#' #
#'  # -----------------------------------------------#
#' #
#'  # test that output is Granges
#' test_that("test that associateReg2Gene outputs correct GRanges object form", {
#'              
#'              
#'     correct_names <- c("reg","name","name2","n","coefs",
#'                        "pval","pval2","qval","qval2")
#'              
#'           # the same for GRanges
#'              r0 <- associateReg2Gene(gr0,asGInteractions = F)  
#'              r1 <- associateReg2Gene(gr1,asGInteractions = F)  
#'              r2 <- associateReg2Gene(gr2,asGInteractions = F)  
#'              r3 <- associateReg2Gene(gr3,asGInteractions = F)  
#'              r4 <- associateReg2Gene(gr4,asGInteractions = F)  
#'              r <- associateReg2Gene(gr,asGInteractions = F)  # test >1 enhancer
#'              r0a <- associateReg2Gene(gr0a)  # test diff gene~enh pair
#'              
#'           # the same for GRangesLists
#'           # test output for GRList
#'              grl0 <- associateReg2Gene(GRangesList(gr0),asGInteractions = F)
#'              grl1 <- associateReg2Gene(GRangesList(gr0,gr1),
#'                                        asGInteractions = F) # GRList >1 GRanges
#'              grl.all <- associateReg2Gene(GRangesList(gr0,gr1,gr2,gr3,gr4),
#'                                       asGInteractions = F) # GRList >1 GRanges
#'              
#'              
#'     #################################
#'     # --------TEST------------------#      
#'              
#'          expect_is(r0, 'GRanges')
#'          expect_is(r1, 'GRanges')
#'          expect_is(r2, 'GRanges')
#'          expect_is(r3,'GRanges')
#'          expect_is(r4, 'GRanges')
#'          expect_is(grl0, 'GRanges')
#'          expect_is(grl1,'GRanges')
#'          expect_is(grl.all,'GRanges')
#'          expect_is(r, 'GRanges')
#'          expect_is(r0a, 'GRanges')
#'              
#'              # test names
#'          expect_equal(names(mcols(r0)),correct_names)
#'          expect_equal(names(mcols(r1)),correct_names)
#'          expect_equal(names(mcols(r2)),correct_names)
#'          expect_equal(names(mcols(r3)),correct_names)
#'          expect_equal(names(mcols(r4)),correct_names)
#'          expect_equal(names(mcols(grl0)),correct_names)
#'          expect_equal(names(mcols(grl1)),correct_names)
#'          expect_equal(names(mcols(grl.all)),correct_names)
#'          expect_equal(names(mcols(r)),correct_names)
#'          expect_equal(names(mcols(r0a)),correct_names)
#'              
#'            })
#' 
#' 
#' 
#' test_that("test that associateReg2Gene outputs correct GInteractions 
#'           class form", {
#'             
#'             # expected metadata
#'             correct_names_GI <- c("name","name2","n","coefs",
#'                                   "pval","pval2","qval","qval2")
#'             
#'             
#'             r0 <- associateReg2Gene(gr0)  # test output for Granges
#'             # the same for GRangesLists
#'             grl0 <- associateReg2Gene(GRangesList(gr0)) # test output for GRList
#'             # GRList >1 GR
#'             grl.all <- associateReg2Gene(GRangesList(gr0,gr1,gr2,gr3,gr4))
#'             #  for >1 enhancer
#'             r <- associateReg2Gene(gr)  # test output for Granges
#'             #  different gene~enh pair
#'             r0a <- associateReg2Gene(gr0a)  # test output for Granges
#'             
#'             
#'   #################################
#'   # --------TEST------------------#
#'             
#'             # test objec type
#'         expect_is(r0, 'GInteractions',"GRanges failure")
#'         expect_is(grl0, 'GInteractions',"GRangesList failure")
#'         expect_is(grl.all, 'GInteractions')
#'         expect_is(r, 'GInteractions',"more >1 enh failure")
#'         expect_is(r0a, 'GInteractions')
#'             
#'             # test names
#'         expect_equal(names(mcols(r0)),correct_names_GI,"GRanges names 
#'                                               failure for GInteractions")
#'         expect_equal(names(mcols(grl.all)),correct_names_GI,"GRangesList names 
#'                                                     failure for GInteractions")
#'         expect_equal(names(mcols(r)),correct_names_GI,"GRanges names 
#'                                 failure for GInteractions when >1 enh present")
#'         expect_equal(names(mcols(r0a)),correct_names_GI)
#'             #
#'           }) 
#' 
#' 
#' test_that("test that associateReg2Gene correctly prefilters data", {
#' 
#'    #############################
#'  # check for GRanges
#' #
#'      r0 <- associateReg2Gene(gr0)  # test output for Granges, should return V
#'      r1 <- associateReg2Gene(gr1)  # gene NA, should return NA
#'      r2 <- associateReg2Gene(gr2)  # same cell type NA, returns value
#'      r3 <- associateReg2Gene(gr3)  # diff cell NA, NA
#'      r4 <- associateReg2Gene(gr4)  # enh contain 1 NA, NA
#'      r0a <- associateReg2Gene(gr0a)# test output for Granges, should return V
#'      r <- associateReg2Gene(gr)    # 2 enh contain, should return V
#' #
#' #
#'      # a. these examples should be filtered prior modelling thus NA is reported
#'            expect_equal(r1$coefs, NA)
#'            expect_equal(r3$coefs, NA)
#'            expect_equal(r4$coefs, NA)
#' #
#' #
#'      # b. these examples should NOT be filtered prior modelling
#'      # thus NA should NOT be reported
#' #
#'            expect_failure(expect_equal(r0$coefs, NA))
#'            expect_failure(expect_equal(r0a$coefs, NA)) # another gene~enh ex
#'            expect_failure(expect_equal(r2$coefs, NA))
#'            expect_false(all(is.na(r$coefs))) # 2 enh example
#' #
#' #
#'  #############################
#'  # check for GRangesList
#' #
#'      grl0 <- associateReg2Gene(GRangesList(gr0)) # test output for GRList
#'      grl1 <- associateReg2Gene(GRangesList(gr1)) # test output for GRList
#'      grl2 <- associateReg2Gene(GRangesList(gr0,gr1)) # GRList >1 GRanges
#'      grl2a <- associateReg2Gene(GRangesList(gr1,gr0)) # GRList >1 GRanges
#'      grl.all <- associateReg2Gene(GRangesList(gr0,gr1,gr2,gr3,gr4))
#' 
#'      
#'      
#'   #################################
#'   # --------TEST------------------#
#'      
#'      
#'    # a. these examples should be filtered prior modelling thus NA is reported
#'            expect_failure(expect_equal(grl0$coefs, NA))
#' 
#'    # b. these examples should NOT be filtered prior modelling
#'    # thus NA should NOT be reported
#'            expect_equal(grl1$coefs, NA)
#' 
#'    # expect both FALSE&TRUE entries
#'            expect_equal(is.na(grl2$coefs),c(FALSE,TRUE))
#'            expect_equal(is.na(grl2a$coefs),c(TRUE,FALSE)) #other way around
#'            expect_equal(is.na(grl.all$coefs),c(FALSE,TRUE,FALSE,TRUE,TRUE))
#' 
#' 
#'  })
#' 
#' 
#'  
#'  
#'  
#'  
#'  #
#' #
#' #
#' #
#' #
#' #
#' #
#' #
#' #
#' #
#' #
#' #
#' # testing_that_code runs correctly for different modelling procedures
#' # testing tag=NULL,scaleData=TRUE,cores=1,B=1000
#' # test that grangesList is reordered correctly - both - correct positions of enhancer and genes and all other elements
#' # testing that there is a failure when no featureType        name       name2
#' # testing that zeroVar works correctly and other function as well
#' # testing stepwise()
#' # scaling
#' # gamma calculations
#' # etc...
#' #
#' #
#' #
#' #
#' #
#' #
#' #
#'  set.seed(153)
#'  associateReg2Gene(gr0,cores = 1)
#' #
#'  associateReg2Gene(GRangesList(gr1),cores = 1,tag="TAGADDED")
#'  associateReg2Gene(GRangesList(gr1),cores = 1,tag="TAGADDED",B=500)
#'  associateReg2Gene(GRangesList(gr0),cores = 1)
#' #
#'  # -----------------------------------------------#
#'  # usage
#' #
#' #
#' #
#' #
#'  # Corr = 0.8 #corr of interest
#'  # N_samples <- 10
#'  # set.seed(46468645)
#'  # mu_mvnorm <- runif(2, min=1.5, max=3.5)
#'  # RandomMatrix <- matrix(Corr,ncol = 2, nrow=2);diag(RandomMatrix) <- 1
#'  #
#'  # # generating a random matrix
#'  # set.seed(46468645)
#'  # corrM = MASS::mvrnorm(n=N_samples,mu=mu_mvnorm,Sigma=RandomMatrix,
#'                       empirical=TRUE)
#'  # corrM <- data.frame(t(corrM),stringsAsFactors = F)
#' 
#' 
#' 
#' 
#' # TO TEST: grlist2gr
#'  regTSS_toy <- GRReg1_toy
#'    regTSS_toy$bw1 <- rep(1,length(GRReg1_toy))
#'    regTSS_toy$bw2 <- rep(2,length(GRReg1_toy))
#'    regTSS_toy$bw3 <- rep(3,length(GRReg1_toy))
#'  regReg_toy <- GRReg2_toy
#'     regReg_toy$bw1 <- rep(3,length(regReg_toy))
#'     regReg_toy$bw2 <- rep(4,length(regReg_toy))
#' #'
#'  grlist2gr(regActivityAroundTSS(regReg_toy,regTSS_toy,upstream=5,downstream=5))
#' #'
#' #'
#' #'
#' #'
#' #'
#' #' grReorg <- function(x){
#' 
#'  grReorg(test[[1]])
#' 
#'  test <- regActivityAroundTSS(regReg_toy,regTSS_toy,upstream=5,downstream=5)
#'  x=test[[1]]
#' x
#' 
#' 
#' 
#' # you can use the same test as written for assocreg2gene before
#' reg=granges(x[x$featureType != "gene",])
#' gene=x[x$featureType == "gene",c("name" ,"name2")]
#' 
#' grpairs=rep(gene,length(reg))
#' reg$reg=grpairs
#' reg$name=grpairs$name
#' reg$name2=grpairs$name2
#' 
#' 
#' reg
#' }