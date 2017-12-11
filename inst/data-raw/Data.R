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
#test2.bw$score <- 2*(1:(length(test.bw)))
test2.bw$score <- c(0,1,0,0,0,1)

  
seqlengths(test2.bw) <- seqlengths(test.bw) <- c(249250621,243199373)


test3.bw <- shift(test2.bw,100)
#rtracklayer::export.bw(test.bw,"/inst/extdata/test.bw")
#rtracklayer::export.bw(test2.bw,"inst/extdata/test2.bw")
#rtracklayer::export.bw(test3.bw,"inst/extdata/test3.bw")



########################################
# Create testing matrix of variables with predefined corr 
# library('MASS')
# library(corpcor)
# library(GenomicRanges)
# 
# r_vector = c(0.99,0.9,0.8,0.6,0.4,0.3,0.1,0)
# samples <- 52
# 
# # generate random seeds for all 14 matrices
# set.seed(390711161)
# 
# mu_mvnorm <- runif(length(r_vector)+1, min=1.5, max=3.5) 
# 
#   # generating a random matrix
#       RandomMatrix <- matrix(0,ncol = length(r_vector)+1, nrow=length(r_vector)+1)
#       diag(RandomMatrix) <- 1
#       set.seed(390711161)
#       RandomMatrix = mvrnorm(n=max_samples, mu=mu_mvnorm,
#                              Sigma=RandomMatrix, empirical=TRUE)
#       
#  
#   # generating matrix with predefined corr
#   # but corr matrix needs to be positive semi-definite matrix, thus
#   # make.positive.definite() is used to compute the nearest positive definite of a real symmetric matrix
#   # but such that wanted correlation remains the same - diag=1, and edges equal to predefined value
#   
#   Var_matrix <- matrix(0,ncol = length(r_vector)+1, nrow=length(r_vector)+1)
#   diag(Var_matrix) <- 1
#   Var_matrix[1,2:ncol(Var_matrix)] <-  r_vector
#   Var_matrix[2:nrow(Var_matrix),1] <-  r_vector
#   
#   
#   # finding random positive-semidefinite correlation matrices with wanted corr
#   set.seed(390711161)      
#   while(is.positive.definite(Var_matrix)!=T){
#     
#     Var_matrix <- make.positive.definite(Var_matrix,tol = )
#     diag(Var_matrix) <- 1
#     Var_matrix[1,2:ncol(Var_matrix)] <-  r_vector
#     Var_matrix[2:nrow(Var_matrix),1] <-  r_vector
#   
#   }
#   set.seed(390711161)      
#   data = mvrnorm(n=samples, mu=mu_mvnorm,Sigma=Var_matrix, empirical=TRUE)
#   
#     data[data<=0] <- 0.00000001
# 
#     cor(data[,1],data[,-1])
#     
#     # Creating GRanges object
#     
#     ModellingTest <- GRanges(rep("chr1",9),IRanges(1:9,2:10),
#                             featureType=c("gene",rep("regulatory",8)),
#                             name=rep("test",9),
#                             name2=rep("test",9))
#                     
#   
#     mcols(ModellingTest) <- cbind(mcols(ModellingTest),DataFrame(t(data)))
#   
#     
# save(ModellingTest,file="data/ModellingTest.RData")
  
