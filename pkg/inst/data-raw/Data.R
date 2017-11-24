################################################################
# Creating example for testing benchmark f()
#   -------Reg1---------------------------------     


GRReg1_toy <- GRanges(c(rep("chr1",8),"chr2"),
                      IRanges(c(1,1,15,24,30,5,100,200,1),
                              c(4,4,20,29,31,9,101,201,4)))


GRReg1_toy$reg <- GRanges(c(rep("chr1",8),"chr2"),
                          IRanges(c(5,30,21,5,31,31,102,202,5),
                                  c(9,31,25,10,32,32,103,203,9)))

GRReg1_toy$name <- GRReg1_toy$name2 <- paste0("TEST_Reg",1:length(GRReg1_toy))
GRReg1_toy$Bench1Exp <- c(2,1,0,0,1,2,3,2,2)
GRReg1_toy$Bench2Exp <- GRReg1_toy$Filter2Exp <- rep(1,9)
GRReg1_toy$Filter1Exp <- c(1,1,0,1,1,1,1,1,1)




GRReg2_toy <- GRanges(c(rep("chr1",12),rep("chr2",2)),
                      IRanges(c(1,1,5,26,1,31,31,200,200,100,100,102,1,5),
                              c(4,3,9,30,3,33,32,201,201,101,101,103,4,9)))

GRReg2_toy$reg <- GRanges(c(rep("chr1",12),rep("chr2",2)),
                          IRanges(c(1,5,1,31,31,5,5,203,203,102,103,100,5,1),
                                  c(3,10,2,32,32,10,10,204,204,103,104,101,9,4)))
GRReg2_toy$name <- GRReg2_toy$name2 <- paste0("TEST_Reg",1:length(GRReg2_toy))


#save(GRReg2_toy,file="pkg/data/GRReg2_toy.RData")
#save(GRReg1_toy,file="pkg/data/GRReg1_toy.RData")
save(GRReg2_toy,file="data/GRReg2_toy.RData")
save(GRReg1_toy,file="data/GRReg1_toy.RData")



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