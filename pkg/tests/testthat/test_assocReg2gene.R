# # ---------------------------------------------------------------------------- #
# # test for associateReg2Gene function
#
require(testthat)
require(energy)
require(glmnet)
require(ranger)
require(qvalue)
require(InteractionSet)


#################################
#------------------------------
# Defining test datasets:
#################################

# DESCRIPTION:
     # 1. gr - example 1 gene + 2 enh
     # all others - 1gene - 1 enhancer
           # 2. gr0 - example 1 gene + 1 enh
           # 3. gr0a - example 1 gene + 1 enh (different enh)
           # 4. gr0.failure -     example 1 gene + 1 enh (gene name repeated,
           #                   failure exp)
           # 6. gr1 - example 1 gene + 1 enh, gene contains 1 NA
           # 7. gr2 - example 1 gene + 1 enh, gene&enh contain 1 NA in the same cell
           # 8. gr3 - example 1 gene + 1 enh, gene&enh contain 1 NA in the diff cell
           # 9. gr4 - example 1 gene + 1 enh, enh contain 1 NA


 x <- 1:5
 y <- 2:6
 z <- 10:14
 k <- c(rep(1,4),4)
 
 GeneInfo <- as.data.frame(matrix(c(rep("gene",3),rep("regulatory",6)),
                                  ncol = 3,byrow = T),stringsAsFactors=FALSE)
             colnames(GeneInfo) <- c("featureType","name","name2")


           gr <- GRanges(seqnames=rep("chr1",3),IRanges(1:3,3:5)) #2Enhancers
           grE2 <- GRanges(seqnames=rep("chr1",3),IRanges(11:13,13:15))#2DiffE
           mcols(grE2) <- DataFrame(cbind(GeneInfo,rbind(x,y,k)))
          mcols(gr) <- DataFrame(cbind(GeneInfo,rbind(x,y,z)))
       # Create different problems in genes and enhancers
           gr0 <- gr[1:2] # one gene - one enh example
           gr0.failure <- gr[2:3]#1G-1E[good example for failure due G-E names]
           gr0a <- gr[c(1,3)] # another one gene - one enh example
           gr4 <- gr3 <-  gr2 <- gr1 <- gr0 # Adding missing data
             mcols(gr1)[1,4] <- NA # 1  gene contains NA
             mcols(gr2)[,4] <- rep(NA,2) # G and E contain NA in the same cell
             mcols(gr3)[1,4] <- NA # G and E contain NA in the different cell
             mcols(gr3)[2,5] <- NA # enh contain NA
             mcols(gr4)[2,5] <- NA
#
# -----------------------------------------------#
#
# test that output is Granges
#
             
test_that("test that associateReg2Gene OUTPUT is in correct form for 1 gene",{

   GR.EXAMPLES <- list(gr,gr0,gr0a,gr1,gr2,gr3,gr4)

   correct_names <- c("reg","name","name2","n","coefs",
                      "pval","qval")

    lapply(GR.EXAMPLES,function(x){
     
     assocR <- associateReg2Gene(x,asGInteractions = F)
     # test object type
     expect_is(assocR, 'GRanges')
     # test names
     expect_equal(names(mcols(assocR)),correct_names)

    })

    # FAILURE EXAMPLE:
    # expect failure due to problems in featureType
        expect_null(associateReg2Gene(gr0.failure))
        expect_message(associateReg2Gene(gr0.failure))
              #"Error in INPUT,featureType should contain 1 gene + min 1 enh")

    # if object is not GRanges nor GRangesList expect error
        expect_error(associateReg2Gene(GInteractions(gr,shift(gr,1))))
           })

test_that("test that associateReg2Gene OUTPUT is correct form for >1 gene", {

 # create GRangesList
  grl0 <- GRangesList(gr0)
  grl1 <- GRangesList(gr0,gr1) # GRList >1 GRanges
  grl.all <- GRangesList(gr0,gr1,gr2,gr3,gr4)# GRList >1 GRanges

    GR.EXAMPLES <- list(grl0,grl1,grl.all)

  correct_names <- c("reg","name","name2","n","coefs",
                     "pval","qval")

  lapply(GR.EXAMPLES,function(x){

    assocR <- associateReg2Gene(gr0,asGInteractions = F)

    # test object type
    expect_is(assocR, 'GRanges')
    # test names
    expect_equal(names(mcols(assocR)),correct_names)

  })

  # test failure example as list()
  expect_null(associateReg2Gene(GRangesList(gr0.failure)))
  expect_message(associateReg2Gene(GRangesList(gr0.failure)))
  

  })

test_that("test that associateReg2Gene outputs correct GInteractions
          class form", {

           # creating GRangesLists:
            grl0 <- GRangesList(gr0) # test output for GRList
            # GRList >1 GR
            grl.all <- GRangesList(gr0,gr1,gr2,gr3,gr4)
            
            # expected metadata
       correct_names_GI <- c("name","name2","n","coefs","pval","qval")

            GR.EXAMPLES <- list(gr,gr0,gr0a,grl0,grl.all)

     # TEST!!!
     lapply(GR.EXAMPLES,function(x){

              assocR <- associateReg2Gene(x,asGInteractions = T)

              # test object type
              expect_is(assocR, 'GInteractions')
              # test names
              expect_equal(names(mcols(assocR)),correct_names_GI)

            })

     # FAILURE EXAMPLE:
     # expect failure due to problems in featureType
     expect_null(associateReg2Gene(gr0.failure,asGInteractions = T))
     expect_message(associateReg2Gene(gr0.failure,asGInteractions = T))
     
    
     # test failure example as list()
     expect_null(associateReg2Gene(GRangesList(gr0.failure),
                                    asGInteractions = T))
     expect_message(associateReg2Gene(GRangesList(gr0.failure),
                                                  asGInteractions = T))
     
     
               })

test_that("test that associateReg2Gene correctly prefilters data - 1gene", {

  # examples which should be successful
  GR.EXAMPLES.success <- list(gr,gr0,gr0a,gr2)

  # TEST!!!
  lapply(GR.EXAMPLES.success,function(x){

    assocR <- associateReg2Gene(x,asGInteractions = T)
    
    expect_failure(expect_equal(assocR$coefs, NA))

  })

  # examples which should return NA
  GR.EXAMPLES.NA <- list(gr1,gr3,gr4)
  #gr1: gene NA, should return NA; gr3: diff cell NA, NA;
  #gr4: enh contain 1 NA, NA

   lapply(GR.EXAMPLES.NA,function(x){
    
    assocR <- associateReg2Gene(x,asGInteractions = T)
    expect_equal(assocR$coefs, NA)
    expect_failure(expect_is(assocR$coefs, "numeric"))

  })

})

test_that("test that associateReg2Gene correctly prefilters data >1gene", {

  grl0 <- GRangesList(gr0) # test output for GRList
  grl1 <- associateReg2Gene(GRangesList(gr1)) # test output for GRList
  grl2 <- associateReg2Gene(GRangesList(gr0,gr1)) # GRList >1 GRanges
  grl2a <- associateReg2Gene(GRangesList(gr1,gr0)) # GRList >1 GRanges
  grl.all <- associateReg2Gene(GRangesList(gr0,gr1,gr2,gr3,gr4))

    #################################
  # --------TEST------------------#


   # a. these examples should be filtered prior modelling thus NA is reported
           expect_failure(expect_equal(grl0$coefs, NA))

   # b. these examples should NOT be filtered prior modelling
   # thus NA should NOT be reported
           expect_equal(grl1$coefs, NA)

   # expect both FALSE&TRUE entries
           expect_equal(is.na(grl2$coefs),c(FALSE,TRUE))
           expect_equal(is.na(grl2a$coefs),c(TRUE,FALSE)) #other way around
           expect_equal(is.na(grl.all$coefs),c(FALSE,TRUE,FALSE,TRUE,TRUE))
#
#
 })

test_that("test that associateReg2Gene correctly assigns E to G", {

  grl.all <- GRangesList(gr0,gr0a,gr0.failure,gr1,gr2,gr3,gr4,grE2,gr)

  # GENERATING NAMES
  
        set.seed(987654)
        GenesRandom <-  c("gr0","gr0a","gr0.failure","gr1","gr2","gr3","gr4",
                          "grE2","gr")

  # mutating GE example to get gene names
        grl.all <-  lapply(seq_along(1:length(grl.all)), function(x){
      
            tmp <- grl.all[[x]]
      
            tmp$name[1] <- GenesRandom[x]
      
          return(tmp)
  })
        names(grl.all) <- GenesRandom

          TestGE.Names <- associateReg2Gene(GRangesList(grl.all))

    ###########
    # TESTING    
    
    # expected filtering results: if wrong DATA TYPE - filter out COMPLETELY!!
          
    expect_false(all(GenesRandom%in%(TestGE.Names$name)),"gr0.failure should 
                                                    be missing - filtered OUT")
    expect_true(all(GenesRandom%in%c(TestGE.Names$name,"gr0.failure")),
                                "gr0.failure should be missing - filtered OUT")
    
    # cases that should succed & fail
    expect_false(any(TestGE.Names$name[complete.cases(TestGE.Names$n)]%in%
                       c("gr1","gr3",'gr4')))
    expect_true(any(TestGE.Names$name[!complete.cases(TestGE.Names$n)]%in%
                       c("gr1","gr3",'gr4')))

    # expected gene~enhancer order ACROSS GENES
    expect_equal(first(TestGE.Names),GRanges(rep("chr1",10),
                                          IRanges(c(2,3,2,2,2,2,12,13,2,3)
                                          ,c(4,5,4,4,4,4,14,15,4,5))))  
    expect_equal(TestGE.Names$name,c("gr0","gr0a","gr1","gr2","gr3",
                                     "gr4","grE2","grE2","gr","gr"))
          

    # WITHIN GENE:
    # test that order assignments to enhancers is correct:
     # 1st example is expected to be different numbers,2 enh have != acitivity
     # 2nd example is expected to be equal numbers, 2 enh have == acitivity
    
    TestGE.Names$coefs <- round(TestGE.Names$coefs,2)
    expect_equal(TestGE.Names$coefs[TestGE.Names$name=="grE2"],c(1,0.71))
    expect_equal(TestGE.Names$coefs[TestGE.Names$name=="gr"],c(1,1))
    
    })       
          
#################################      
# TEST MODELLING:
# defining test data:

      x <- c(2.000346,2.166255,0.7372374,0.9380581,2.423209, 
             2.599857,4.216959,2.589133,1.848172,3.039659)
      y <- c(2.866875,2.817145,2.1434456,2.9039771,3.819091,5.009990,
             5.048476,2.884551,2.780067,4.053136)
      z <- 1:length(y)
      
      corrM <- rbind(x,y)
      GR <- GRanges(seqnames=rep("chr1",2),IRanges(1:2,3:4))
      
      GeneInfo <- as.data.frame(matrix(rep(c("gene","regulatory"),each=3),
                                       ncol = 3,byrow = TRUE),stringsAsFactors=FALSE)
      colnames(GeneInfo) <- c("featureType","name","name2")
      mcols(GR) <- DataFrame(cbind(GeneInfo,corrM))


      # create a list of genes
      GR1 <- GR2 <- GR
      mcols(GR) <- DataFrame(cbind(GeneInfo,rbind(x,y)))
      mcols(GR1) <- DataFrame(cbind(GeneInfo,rbind(x,z)))
      mcols(GR2) <- DataFrame(cbind(GeneInfo,rbind(y,z)))
      GR.list <- GRangesList(GR1,GR2,GR)
      

test_that("test that associateReg2Gene gamma P-value is correctly calculated", {
    
  # Function to test Pvalue from associateReg2Gene and recalculated p-value
  Test_B_pvalue <-  function(GR=GR,
                             diffB=5,
                             seeding=818440,
                             method="pearson"){
    
    corrM <-  as.data.frame(mcols(GR))[,-c(1:3)]
    corrM <- apply(corrM,2,as.numeric)
    
    # Function rewritten
    set.seed(seeding)
    #get resampling Ys
    Ys=replicate(diffB,sample(corrM[1,],ncol(corrM)))
    
    if (method!="dcor"){
          RepeatCorr<-as.vector(abs(cor(Ys,corrM[2,],method=method)))
          orig <- cor(corrM[1,],corrM[2,],method=method)
          }
    
    if (method=="dcor"){ 
          RepeatCorr<- apply(Ys,2, function(y){energy::dcor(corrM[2,],y)})
          orig <- energy::dcor(corrM[1,],corrM[2,])}
    
    set.seed(seeding)
    param <-  fitdistrplus::fitdist(RepeatCorr,"gamma",method="mme") # fit gamma
    # p-val based on gamma
    set.seed(seeding)
    pval=1-pgamma(orig,param$estimate[1],param$estimate[2])
    
    
    #######################
    # TESTING
    set.seed(seeding)
    expect_equal(pval,associateReg2Gene(GR,B=diffB,method=method)$pval)
    
    
  }      
  
  ####################################
  # ------------------TEST ----------
  ###################################
  
  # test P-value from gamma distribution is correctly calculated
       Test_B_pvalue(GR=GR,diffB=10,seeding=818440)
  
  # Test it is true for different B's
        B=c(5,10,15,100,1000)
          lapply(B,function(x) {Test_B_pvalue(GR=GR,diffB=x,seeding=818440)})
  
  # Test that it is true regardelss of seeds used
        seeds=c(6548,22882,1006666,13787)
          lapply(seeds,function(x) {Test_B_pvalue(GR=GR,diffB=100,seeding=x)})
  
  # Test that it is true regardelss of method used
        meth <- c("pearson","spearman","dcor")
        lapply(meth,function(x) {Test_B_pvalue(GR=GR,
                                                diffB=100,
                                                seeding=818440,
                                                method=x)})
        
  # testing that it is true for diff correlation example
        Test_B_pvalue(GR=GR1,diffB=10,seeding=818440)
        Test_B_pvalue(GR=GR2,diffB=10,seeding=818440)
        
        
  # testing that it is true for list a of genes
       set.seed(818440)
       expect_equal(round(associateReg2Gene(GR.list)$pval,3),
                    c(0.125,0.220,0.018))
                                
       
    # testing that it is true for >2 enhancer
       GRE3 <- c(GR,shift(GR1[2],3),shift(GR[2],4))
       
       set.seed(818440)
       expect_equal(round(associateReg2Gene(GRE3)$pval,3),
                    c(0.019,0.125,0.019))
       
})
  
test_that("test that associateReg2Gene coef are correctly calculated", {
  
  # test coef are correctly calculated for GRanges
  expect_equal(associateReg2Gene(GR)$coefs,0.8) # test for 1 GENE
  expect_equal(round(associateReg2Gene(GR1)$coefs,3),0.484) # test 1 gene/diff corr
  # test coef are correctly calculated for GRangesList
  expect_equal(round(associateReg2Gene(GR.list)$coefs,3),c(0.484,0.399,0.8))
  
  # test coef are correctly calculated across methods:
  meth=c("pearson",'spearman',"elasticnet","dcor","randomForest")
 
   CoefsTested <- sapply(meth,function(x){associateReg2Gene(GR,method=x)$coefs})
      expect_equal(as.numeric(round(CoefsTested,2)),c(0.80,0.83,NA,0.77,NA))
      
  # test coef are correctly calculated across methods + diff X~Y relationship
  CoefsTested <- sapply(meth,function(x){associateReg2Gene(GR2,method=x)$coefs})
      expect_equal(as.numeric(round(CoefsTested,2)),c(0.40,0.38,NA,0.54,NA))
  
  
  # testing that it is true for >2 enhancer
  GRE3 <- c(GR,shift(GR1[2],3),shift(GR[2],4))
  expect_equal(round(associateReg2Gene(GRE3)$coefs,3),
               c(0.80,0.484,0.80))
  
})

test_that("test that associateReg2Gene qval are correctly calculated", {
  
  data("ModellingTest")
  set.seed(665454)
  P <- associateReg2Gene(ModellingTest)
  expect_equal(qvalue::qvalue(P$pval)$qvalue,P$qval)
  expect_equal(round(P$qval,5),c(0.00001,0.00002,0.00006,0.00106,
                                 0.01550,0.05145,0.49656,0.99948))
  
  
})
  
test_that("test that associateReg2Gene elasticnet/rf works as expected", {

  require(glmnet)
  require(ranger)
  data("ModellingTest")

  set.seed(7777)
  expect_equal(round(associateReg2Gene(ModellingTest,"elasticnet",B=15)$coefs,3),
  c(0.599,0.170,0.156,0.102,0.073,0.045,0.009,0.000))
  
  set.seed(7777)
  expect_equal(round(associateReg2Gene(ModellingTest,"elasticnet",B=25)$coefs,3),
               c(0.599,0.170,0.156,0.102,0.073,0.045,0.009,0.000))
  
  # Note, depending on the B agrument used, slightly different results
   set.seed(7777)
  expect_failure(expect_equal(
     round(associateReg2Gene(ModellingTest,"elasticnet",B=50)$coefs,3),
                c(0.599,0.170,0.156,0.102,0.073,0.045,0.009,0.000)))
   

  set.seed(7777)
  expect_equal(round(associateReg2Gene(ModellingTest,"randomForest",B=15)$coefs,3),
               c(13.490,12.507,9.395,5.893,2.607,2.187,1.137,1.568))
  
  
  })
  
  
  
 