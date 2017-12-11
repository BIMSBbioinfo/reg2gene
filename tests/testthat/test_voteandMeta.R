
library(testthat)
require(GenomicRanges)
require(InteractionSet)

######################
# CREATING EXAMPLE:

gr2 <- gr <- GRanges(seqnames=rep("chr1",3),IRanges(1:3,3:5))
      
      x <- 1:5
      y <- 2:6
      z <- 10:14
      a <- rep(0,length(x))
      GeneInfo <- as.data.frame(matrix(c(rep("gene",3),rep("regulatory",6)),
                                ncol = 3,byrow = TRUE),stringsAsFactors=FALSE)
      colnames(GeneInfo) <- c("featureType","name","name2")

  mcols(gr) <- DataFrame(cbind(GeneInfo,rbind(x,y,z)))
  mcols(gr2) <- DataFrame(cbind(GeneInfo,rbind(x,y,a)))

# create associateReg2Gene output objects, GInteractions will all 
# output results

    AssocObject <- reg2gene::associateReg2Gene(gr)
    AssocObject2 <- reg2gene::associateReg2Gene(gr2)

# input for meta-analysis is list of such objects

associations <- list(AssocObject,AssocObject2)
names(associations) <- c("AssocObject","AssocObject2")


test_that("test that voteAssociations IN/OUTPUT is in correct form",{
  
    expect_is(associations,"list") # expected input is list of
    expect_equal(unique(sapply(associations,class)),"GInteractions") # GIs
  
  # not working for individual GInteractions or GRange
    expect_error(voteAssociations(AssocObject))
    expect_error(voteAssociations(AssocObject2))
    expect_error(voteAssociations(gr))
  
  # CHECKING output
  expect_is(voteAssociations(associations),"GInteractions")
  
  # CHECKING output - for objects without overlap empty GI
    expect_is(voteAssociations(list(associations$AssocObject[1],
                          associations$AssocObject[2])),"GInteractions")
  
  #with different threshold still GI is returned
    expect_is(voteAssociations(list(associations$AssocObject[1],
                                    associations$AssocObject[2]),
                               vote.threshold=0.51),"GInteractions")
    
    
    ##################################
    # Test if true for:
    #  example #2
          set.seed(6878); x=rnorm(15)
          set.seed(444);  y=rnorm(15)
          set.seed(6848);  z=rnorm(15)
          
  
          example <-  example2 <-   GRanges(GRReg1_toy[1:2],
                                          featureType=c("gene","regulatory"),
                                          name=c("gene","regulatory"),
                                          name2=c("gene","regulatory"))
                    
                mcols(example2) <-  cbind(mcols(example2)[,1:3],
                                          DataFrame(rbind(x,y)))
                mcols(example) <-  cbind(mcols(example)[,1:3],
                                         DataFrame(rbind(x,z)))
                
          AssocExample2 <- associateReg2Gene(example2)
          AssocExample <- associateReg2Gene(example)
     
    ######################################
    # TESTING
          
     # both ranges fail at filtering    
        expect_error(voteAssociations(list(AssocExample2,AssocExample)))
     # should succed with changed threshold to 1 both examples are succesfull
        expect_equal(length(voteAssociations(list(AssocExample2,AssocExample),
                                           cutoff.val = 1)),1)
    
      # Equal results with different order:  
        # both ranges fail at filtering    
        expect_error(voteAssociations(list(AssocExample,AssocExample2)))
        # should succed with changed threshold to 1 both examples are succesfull
        expect_equal(length(voteAssociations(list(AssocExample,AssocExample2),
                                             cutoff.val = 1)),1)  
        
      # MUTATING: one range fails at filtering    
          AssocExample2$pval <- 0.01
      expect_equal(length(voteAssociations(list(AssocExample2,AssocExample))),0)
    
      # MUTATING: TWO ranges are filtered    
        AssocExample$pval <- 0.01
        expect_equal(length(voteAssociations(list(AssocExample,AssocExample2),
                                             cutoff.val = 1)),1)  
        
        
        # MUTATING:   cutoff.stat -> qval (all NA)
        expect_error(voteAssociations(list(AssocExample2,AssocExample),
                                      cutoff.stat = "qval"))
        
        # MUTATING:   cutoff.stat -> qval (all 0.1)
        AssocExample2$qval <-  AssocExample$qval <- 0.01
        expect_equal(length(voteAssociations(list(AssocExample2,AssocExample),
                                      cutoff.stat = "qval")),1)
        
        
        # testing cutoff.stat: if NOT existing
        expect_error(voteAssociations(list(AssocExample2,AssocExample),
                                      cutoff.stat="blablacar"))
        
        
        
        })


test_that("test that voteAssociations output is in correct values",{
  
  
  # CHECKING output
  expect_equal(voteAssociations(associations)$votes,2)
  
  # CHECKING output >1 overlap
  expect_equal(voteAssociations(list(AssocObject,
                                     AssocObject,
                                     AssocObject2))$votes,c(3,2))
  # MUTATING: check: vote.threshold
          # CHECKING output >1 overlap: different threshold, above 2 expected
          expect_equal(voteAssociations(list(AssocObject,
                                             AssocObject,
                                             AssocObject2),
                                        vote.threshold = 0.66)$votes,c(3,2))
          expect_equal(voteAssociations(list(AssocObject,
                                             AssocObject,
                                             AssocObject2),
                                        vote.threshold = 0.67)$votes,3)
         # even if threshold 0, all entries with vote 1 are filtered out 
          expect_equal(voteAssociations(list(AssocObject,
                                             AssocObject,
                                             AssocObject2),
                                        vote.threshold = 0)$votes,c(3,2))
  
  
  # CHECKING output >1 overlap: different order
  expect_equal(voteAssociations(list(AssocObject2,
                                     AssocObject2,
                                     AssocObject))$votes,3)
  
  ##################################
  # Test if true for:
  #  example #2
  set.seed(6878); x=rnorm(15)
  set.seed(444);  y=rnorm(15)
  set.seed(6848); z=rnorm(15)
  
  
          example <-  example2 <-   GRanges(GRReg1_toy[1:2],
                                            featureType=c("gene","regulatory"),
                                            name=c("gene","regulatory"),
                                            name2=c("gene","regulatory"))
          
          mcols(example2) <-  cbind(mcols(example2)[,1:3],
                                    DataFrame(rbind(x,y)))
          mcols(example) <-  cbind(mcols(example)[,1:3],
                                   DataFrame(rbind(x,z)))
          
          AssocExample2 <- associateReg2Gene(example2)
          AssocExample <- associateReg2Gene(example)
          
  ######################################
  # TESTING
          
    expect_error(expect_equal(voteAssociations(list(AssocExample2,
                                                    AssocExample))$votes,2))
  
 # MUTATING: cutoff.val, if set 1, success
          
      expect_equal(voteAssociations(list(AssocExample2,AssocExample),
                                   cutoff.val=1)$votes,2)
      
          
   # MUTATING: one range fails at filtering    
        AssocExample2$pval <- 0.01
       
    # MUTATING: TWO ranges are filtered    
        AssocExample$pval <- 0.01
        expect_equal(voteAssociations(list(AssocExample2,AssocExample))$votes,2)
  
  # MUTATING:   cutoff.stat
      AssocExample2$qval <-  AssocExample$qval <- 0.01
        expect_equal(voteAssociations(list(AssocExample2,AssocExample),
                                cutoff.stat = "qval")$votes,2)

  # MUTATING: cutoff.val, if set >0.01, failure again
        
        expect_error(expect_equal(voteAssociations(list(AssocExample2,
                                  AssocExample),cutoff.val=0.001)$votes,2))    
        
  

})





############################################################3
# TEST META-ANALYSIS


test_that("test that metaAssociations IN/OUTPUT is in correct form",{
  
  expect_is(associations,"list") # expected input is list of
  expect_equal(unique(sapply(associations,class)),"GInteractions") # GIs
  
  # not working for individual GInteractions or GRange
  expect_error(metaAssociations(AssocObject))
  expect_error(metaAssociations(AssocObject2))
  expect_error(metaAssociations(gr))
  
  # CHECKING output
  expect_is(metaAssociations(associations),"GInteractions")
  
  # CHECKING output - for objects without overlap - NULL
  expect_null(metaAssociations(list(associations$AssocObject[1],
                                  associations$AssocObject[2])))
  
  
  ##################################
  # Test if true for:
  #  example #2
  set.seed(6878); x=rnorm(15)
  set.seed(444);  y=rnorm(15)
  set.seed(6848);  z=rnorm(15)
  
  
  example <-  example2 <-   GRanges(GRReg1_toy[1:2],
                                    featureType=c("gene","regulatory"),
                                    name=c("gene","regulatory"),
                                    name2=c("gene","regulatory"))
  
  mcols(example2) <-  cbind(mcols(example2)[,1:3],
                            DataFrame(rbind(x,y)))
  mcols(example) <-  cbind(mcols(example)[,1:3],
                           DataFrame(rbind(x,z)))
  
  AssocExample2 <- associateReg2Gene(example2)
  AssocExample <- associateReg2Gene(example)
  
  ######################################
  # TESTING
  
  # both ranges fail at filtering    
  expect_is(metaAssociations(list(AssocExample2,AssocExample)),"GInteractions")
  # Equal results with different order:  
  expect_is(metaAssociations(list(AssocExample,AssocExample2)),"GInteractions")
  
 
  
  
})

test_that("test that metaAssociations output is in correct values",{
  
  
  # CHECKING output
  set.seed(4555)
 
  expect_equal(metaAssociations(associations)$coefs,c(1,0.5))
  expect_equal(metaAssociations(associations)$n,c(10,10))
  
  # CHECKING output >1 overlap
  expect_equal(metaAssociations(list(AssocObject,
                                     AssocObject,
                                     AssocObject2))$n,c(15,15))
  
  expect_equal(round(metaAssociations(list(AssocObject,
                                     AssocObject,
                                     AssocObject2))$coefs,2),c(1.00,0.67))
  

  
  #################################
  # Test if true for:
  #  example #2
  set.seed(6878); x=rnorm(15)
  set.seed(444);  y=rnorm(15)
  set.seed(6848); z=rnorm(15)
  
  
  example <-  example2 <-   GRanges(GRReg1_toy[1:2],
                                    featureType=c("gene","regulatory"),
                                    name=c("gene","regulatory"),
                                    name2=c("gene","regulatory"))
  
  mcols(example2) <-  cbind(mcols(example2)[,1:3],
                            DataFrame(rbind(x,y)))
  mcols(example) <-  cbind(mcols(example)[,1:3],
                           DataFrame(rbind(x,z)))
  
  AssocExample2 <- associateReg2Gene(example2)
  AssocExample <- associateReg2Gene(example)
  
  ######################################
  # TESTING
  
  expect_equal(metaAssociations(list(AssocExample2,AssocExample))$coefs,
               -0.07232527)
  expect_equal(metaAssociations(list(AssocExample2,AssocExample))$n,30)
  
  # check procedure of calculating coeficients
  expect_equal(metaAssociations(list(AssocExample2,AssocExample))$coefs,
            sum(AssocExample2$n*(AssocExample2$coefs)+(AssocExample$n*
                        AssocExample$coefs))/(AssocExample$n+AssocExample2$n))
  
  
  expect_equal(metaAssociations(list(AssocExample2,AssocExample))$coefs,
               sum(15*-0.1006173+15*-0.04403323)/30)
  
    comb <- -2 * sum(log(c(0.766882,0.9262192)))
    pval=1-pchisq(comb,df=4)
  
    expect_equal(metaAssociations(list(AssocExample2,AssocExample))$pval,pval)
    
  
})









