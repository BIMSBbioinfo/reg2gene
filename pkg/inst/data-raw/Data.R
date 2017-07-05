################################################################
# Create examples to test OverlapRegions() and ComplexOverlaps() 

MinNames <- c("chr1","start1","end1","chr2","start2","end2")

# creating example for testing ComplexOverlaps
#   -------Reg1---------------------------------     

GRReg1_toy <- GRanges(c(rep("chr1",4),rep("chr2",2)),
                      IRanges(c(1,1,1,15,1,15),c(4,4,4,20,4,20)))
GRReg1_toy$reg <- GRanges(c(rep("chr1",4),rep("chr2",2)),
                          IRanges(c(5,5,30,21,5,21),c(6,9,31,25,9,25)))
GRReg1_toy$name <- GRReg1_toy$name2 <- paste0("TEST_Reg",1:length(GRReg1_toy))



GRReg2_toy <- GRanges(c(rep("chr1",7),rep("chr2",7)),
                      IRanges(c(1,5,26,1,31,5,31,1,5,26,1,31,5,31),
                            c(3,9,30,3,32,10,32,3,9,30,3,32,10,32)))
GRReg2_toy$reg <- GRanges(c(rep("chr1",7),rep("chr2",7)),
                          IRanges(c(5,1,31,31,1,31,5,5,1,31,31,1,31,5),
                                 c(10,2,32,32,3,32,10,10,2,32,32,3,32,10)))
GRReg2_toy$name <- GRReg2_toy$name2 <- paste0("TEST_Reg",1:length(GRReg2_toy))


save(GRReg2_toy,file="data/GRReg2_toy.RData")
save(GRReg1_toy,file="data/GRReg1_toy.RData")




################################################################
# Create examples to test ConfusionMatrixReg2Gene()

BenchMarkedReg2Gene_toy <- GRReg1_toy
BenchMarkedReg2Gene_toy$BenchmarkO <- c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE)
BenchMarkedReg2Gene_toy$PValue <- c(0.05,0.06,0.01,0.5,0,1)


save(BenchMarkedReg2Gene_toy,file="data/BenchMarkedReg2Gene_toy.RData")

################################################################
# Create examples to test Filter_PreModelling()

GR_exp_toy <- GRanges(rep("chr1",4),IRanges(1:4,2:5))
 
mcols(GR_exp_toy) <-  matrix(rep("test",12),nrow=4,
                             dimnames = list(1:4,c("featureType",
                                                   "name","name2")))

mcols(GR_exp_toy) <- cbind(mcols(GR_exp_toy),data.frame(matrix(c(rep(0,10),
                                                              c(rep(0,9),1),
                                                          c(rep(0,3),rep(1,7)),
                                                              c(rep(1,9),NA)),
                                        byrow = T,nrow=4),stringsAsFactors = F))



################################################################
# Create examples to test Filter_PreModelling()

GR_exp_toy <- GRanges(rep("chr1",4),IRanges(1:4,2:5))

mcols(GR_exp_toy) <-  matrix(rep("test",12),nrow=4,
                             dimnames = list(1:4,c("featureType",
                                                   "name","name2")))

mcols(GR_exp_toy) <- cbind(mcols(GR_exp_toy),data.frame(matrix(
                                              c(rep(0,10),c(rep(0,9),1),
                                              c(rep(0,3),rep(1,7)),
                                              c(rep(1,9),NA)),
                                        byrow = T,nrow=4),stringsAsFactors = F))




################################################################
# Create examples to test bwToGeneExp() and regActivity()

test2.bw <- test.bw <- GRanges(c("chr1","chr1","chr1","chr2","chr2","chr2"),
                               IRanges(c(1,5,10,1,5,10),
                                       c(4,9,25,4,9,25)))

test.bw$score <- c(0,1:(length(test.bw)-1))
test2.bw$score <- 2*(1:(length(test.bw)))

seqlengths(test2.bw) <- seqlengths(test.bw) <- c(249250621,243199373)

#rtracklayer::export.bw(test.bw,"/inst/extdata/test.bw")
#rtracklayer::export.bw(test2.bw,"inst/extdata/test2.bw")