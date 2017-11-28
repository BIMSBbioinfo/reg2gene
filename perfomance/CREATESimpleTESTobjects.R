

# ComplexOverlap and OverlapRegions Examples

              MinNames <- c("chr1","start1","end1","chr2","start2","end2")
              
              # creating example for testing ComplexOverlaps
              
              Reg1.ranges <- matrix(c(1,4,5,9,
                                      15,20,21,25),
                                    byrow = T, ncol=4)
              
              
              Reg1TEST <- cbind(rep("chr1",nrow(Reg1.ranges)),Reg1.ranges[,1:2],
                                    rep("chr1",nrow(Reg1.ranges)),Reg1.ranges[,3:4])
              
              
              
              Reg2.ranges <- matrix(c(1,3,5,10,
                                      5,10,1,3,
                                      26,30,31,32,
                                      1,3,31,32,
                                      31,32,1,3,
                                      5,10,31,32,
                                      31,32,5,10),
                                    byrow = T, ncol=4)
                
              Reg2TEST <- cbind(rep("chr1",nrow(Reg2.ranges)),Reg2.ranges[,1:2],
                                    rep("chr1",nrow(Reg2.ranges)),Reg2.ranges[,3:4])
              ExpectedOverap <- c("1122",
                              "1221",
                              "noO",
                              "11",
                              "12",
                              "21",
                              "12")  
              Reg2TEST <- cbind(Reg2TEST,ExpectedOverap)  
              
          colnames(Reg1TEST)[1:6] <- colnames(Reg2TEST)[1:6] <-  MinNames

              
              
             
