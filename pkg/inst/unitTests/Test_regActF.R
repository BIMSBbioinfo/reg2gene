# 
# 
# 
# # tests for data integration f()
# 
# 
#  test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
#  test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
#  regTSS_toy <- GRanges(c(rep("chr1",4),rep("chr2",2)),
#                        IRanges(c(1,7,9,15,1,15),c(4,8,14,20,4,20)),
#                                              c(rep("+",3),rep("-",3)))
#  regTSS_toy$reg <-  regTSS_toy[c(1,1,3:6)]
#  regTSS_toy$name2 <- regTSS_toy$name <- paste0("TEST_Reg",
#                                          c(1,1,3:length(regTSS_toy)))
#  regActivity(regTSS_toy,c(test.bw,test2.bw)) 
#  
#  
#  
#  
#  
#  
#  
# 
#  
#  
#  
#  test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
#   test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
#   regTSS_toy <- GRanges(c(rep("chr1",4),rep("chr2",2)),
#                         IRanges(c(1,7,9,15,1,15),c(4,8,14,20,4,20)),
#                                               c(rep("+",3),rep("-",3)))
#   regTSS_toy$reg <-  regTSS_toy[c(1,1,3:6)]
#   regTSS_toy$name2 <- regTSS_toy$name <- paste0("TEST_Reg",
#                                           c(1,1,3:length(regTSS_toy)))
#                                           
#                                           
#   bwToGeneExp(exons = regTSS_toy,geneActSignals = c(test.bw,test2.bw))
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  # regActivityAroundTSS
#  
#  
#   regTSS_toy <- GRReg1_toy
#     regTSS_toy$bw1 <- rep(1,length(GRReg1_toy))
#     regTSS_toy$bw2 <- rep(2,length(GRReg1_toy))
#     regTSS_toy$bw3 <- rep(3,length(GRReg1_toy))
#   regReg_toy <- GRReg2_toy
#      regReg_toy$bw1 <- rep(3,length(regReg_toy))
#      regReg_toy$bw2 <- rep(4,length(regReg_toy))
#   
#   regActivityAroundTSS(regReg_toy,regTSS_toy,upstream=1,downstream=1)
#   regActivityAroundTSS(regReg_toy,regTSS_toy,upstream=5,downstream=5)
#   
#   