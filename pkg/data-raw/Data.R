setwd("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/pkg/")
################################################################
# Create examples to test OverlapRegions() and ComplexOverlaps() 

  MinNames <- c("chr1","start1","end1","chr2","start2","end2")
              
# creating example for testing ComplexOverlaps
    #   -------Reg1---------------------------------     
              Reg1.ranges <- matrix(c(1,4,5,9,
                                      15,20,21,25),
                                    byrow = T, ncol=4)
              
              
              Reg1TEST <- cbind(rep("chr1",nrow(Reg1.ranges)),Reg1.ranges[,1:2],
                                    rep("chr1",nrow(Reg1.ranges)),Reg1.ranges[,3:4])
              
              
     #   -------Reg2---------------------------------        
              Reg2.ranges <- matrix(c(1,3,5,10,
                                      5,9,1,2,
                                      26,30,31,32,
                                      1,3,31,32,
                                      31,32,1,3,
                                      5,10,31,32,
                                      31,32,5,10),
                                    byrow = T, ncol=4)
                
              Reg2TEST <- cbind(rep("chr1",nrow(Reg2.ranges)),Reg2.ranges[,1:2],
                                rep("chr1",nrow(Reg2.ranges)),Reg2.ranges[,3:4])
                   ExpectedOverap <- c("1122","1221","noO","11","12","21","12")  
              Reg2TEST <- cbind(Reg2TEST,ExpectedOverap)  
              
          colnames(Reg1TEST)[1:6] <- colnames(Reg2TEST)[1:6] <-  MinNames

          Reg1_toy <- data.frame(Reg1TEST,stringsAsFactors = F)
          Reg2_toy <- data.frame(Reg2TEST,stringsAsFactors = F)
          
         
          
              save(Reg1_toy,file="/data/Reg1_toy.RData")
              save(Reg2_toy,file="/data/Reg2_toy.RData")
      #########################################        
      # Add Granges toy version 
              library(GenomicRanges)
              
              #GRReg1_toy
                      GRReg1_toy <- GRanges(Reg1_toy[,"chr1"],IRanges(as.integer(Reg1_toy[,"start1"]),
                                                                     as.integer(Reg1_toy[,"end1"])))
                      
                    
                      GRReg1_toy$reg <- GRanges(Reg1_toy[,"chr2"],IRanges(as.integer(Reg1_toy[,"start2"]),
                                                                         as.integer(Reg1_toy[,"end2"])))
                      
                      GRReg1_toy$name <- GRReg1_toy$name2 <- "TEST_Reg1"
                      
              #GRReg2_toy
                      GRReg2_toy <- GRanges(Reg2_toy[,"chr1"],IRanges(as.integer(Reg2_toy[,"start1"]),
                                                                      as.integer(Reg2_toy[,"end1"])))
                      
                      
                      GRReg2_toy$reg <- GRanges(Reg2_toy[,"chr2"],IRanges(as.integer(Reg2_toy[,"start2"]),
                                                                          as.integer(Reg2_toy[,"end2"])))
                      
                      GRReg2_toy$name <- GRReg2_toy$name2 <- "TEST_Reg2"
                      
               save(GRReg1_toy,file="/data/GRReg1_toy.RData")
               save(GRReg2_toy,file="/data/GRReg2_toy.RData")
              
################################################################
# Create examples to test OverlapRegions() and ComplexOverlaps() 
              
              Reg1Extended <- rbind(Reg1_toy,
                                    do.call("rbind",lapply(c("chr4","chr2","chr3"),
                                                           function(x){apply(Reg1_toy,2,
                                                                             str_replace,"chr1",x)})))
              Reg1Extended_toy <- data.frame(Reg1Extended,stringsAsFactors = F)
              
             
              
              Reg2Extended <- rbind(Reg2_toy,
                                    do.call("rbind",lapply(c("chr4","chr2","chr3"),
                                                           function(x){apply(Reg2_toy,2,
                                                                             str_replace,"chr1",x)})))
              Reg2Extended_toy <- data.frame(Reg2Extended,stringsAsFactors = F)
              
              
              
          save(Reg2Extended_toy,file="/data/Reg2Extended_toy.RData")
          save(Reg1Extended_toy,file="/data/Reg1Extended_toy.RData")
          
          
        #########################################        
        # Add Granges toy version 
              library(GenomicRanges)
              
              #GRReg1Extended_toy
              GRReg1Extended_toy <- GRanges(Reg1Extended_toy[,"chr1"],IRanges(as.integer(Reg1Extended_toy[,"start1"]),
                                                                              as.integer(Reg1Extended_toy[,"end1"])))
              
              
              GRReg1Extended_toy$reg <- GRanges(Reg1Extended_toy[,"chr2"],IRanges(as.integer(Reg1Extended_toy[,"start2"]),
                                                                                  as.integer(Reg1Extended_toy[,"end2"])))
              
              GRReg1Extended_toy$name <- GRReg1Extended_toy$name2 <- "TEST_Reg1Ext"
              
              #GRReg2Extended_toy
              GRReg2Extended_toy <- GRanges(Reg2Extended_toy[,"chr1"],IRanges(as.integer(Reg2Extended_toy[,"start1"]),
                                                                              as.integer(Reg2Extended_toy[,"end1"])))
              
              
              GRReg2Extended_toy$reg <- GRanges(Reg2Extended_toy[,"chr2"],IRanges(as.integer(Reg2Extended_toy[,"start2"]),
                                                                                  as.integer(Reg2Extended_toy[,"end2"])))
              
              GRReg2Extended_toy$name <- GRReg2Extended_toy$name2 <- "TEST_Reg2Ext"
              
        save(GRReg1Extended_toy,file="/data/GRReg1Extended_toy.RData")
        save(GRReg2Extended_toy,file="/data/GRReg2Extended_toy.RData")
        
        
        GRReg1Extended_toy.2 <- c(GRReg1Extended_toy,GRReg1Extended_toy[1])
        GRReg1Extended_toy.2$reg[1] <- resize(GRReg1Extended_toy$reg[1],2)
        
        
        #save(GRReg1Extended_toy.2,file="/data/akalin/Projects/AAkalin_reg2gene/reg2gene/pkg/data/GRReg2Extended_toy2.RData")
        save(GRReg1Extended_toy.2,file="/data/GRReg2Extended_toy2.RData")

        GRReg1Extended_toy.3 <- c(GRReg1Extended_toy.2,GRReg1Extended_toy.2[1])
        GRReg1Extended_toy.3$reg[10] <- GRanges("chr1",IRanges(30,31))
        
        #save(GRReg1Extended_toy.3,file="/data/akalin/Projects/AAkalin_reg2gene/reg2gene/pkg/data/GRReg2Extended_toy3.RData")
        save(GRReg1Extended_toy.3,file="/data/GRReg2Extended_toy3.RData")
             
################################################################
# Create examples to test ConfusionMatrixReg2Gene functions                            

              
        BenchmarkedReg2Gene_toy <- BenchMarkReg2Gene(GRReg1Extended_toy.3,GRReg2Extended_toy,binary=T)
              
                  Pvalues <- seq(0,1,1/(length(BenchmarkedReg2Gene_toy)-1))
                  BenchmarkedReg2Gene_toy$PValues <- Pvalues
             
             #save(BenchmarkedReg2Gene_toy,file="/data/akalin/Projects/AAkalin_reg2gene/reg2gene/pkg/data/BenchmarkedReg2Gene_toy.RData")
             save(BenchmarkedReg2Gene_toy,file="data/BenchmarkedReg2Gene_toy.RData")
                 

             
           
              
      
             
             