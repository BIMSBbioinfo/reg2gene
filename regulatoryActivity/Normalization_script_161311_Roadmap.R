
# Script that takes COHORT : Blueprint,CEMT,McGill or Roadmap
# and Experiment Type as an Input
# Calculates size factors using DESeq2 estimateSizeFactorsForMatrix()
# And adjust calculated scores for library size
# Outputs both GR object and DF




rm(list = ls())


args <- commandArgs(T)
COHORT <- as.character(args[1])
EXPERIMENT <- as.character(args[2]) # DNAMethylation,DNase-Hypersensitivity,H3K27ac,H3K4me1,RNA-Seq




#####################################
# Setting Log files output
#####################################
sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/Normalization_Scripts/",Sys.Date(),EXPERIMENT,COHORT,"_Norm_OUTPUT.Rout"),
          open="wt"))
print(paste("Normalized COHORT is ",COHORT," with experiment type ",EXPERIMENT, "using version of enhancers/genes"," performed on",Sys.Date()))
sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/Normalization_Scripts/",Sys.Date(),EXPERIMENT,COHORT,"_Norm_ERROR.Rout"),
          open="wt"),type="message",append=TRUE)


# COHORT <- "Blueprint"
# EXPERIMENT<-"H3K4me1"


library(DESeq2)
library(stringr)




INPATH <- "/data/akalin/Projects/AAkalin_Catalog_RI/Data/ScoreMatrixBinResults_Roadmap/"



#####################
# Reading and analysis

# detect all GRanges objects created by Calculating_Scores_Enhancer_Regions.R script
All.experiments.input.path <- list.files(paste0(INPATH,COHORT,"/"),recursive = T,full.names = T,pattern = ".GR.rds")
# for RNASeq skip results reported as exons and orientate on results reportea on the level of genes
All.experiments.input.path <-All.experiments.input.path[!str_detect(All.experiments.input.path,"exon")]

# selecting experiment and reading corresponding GR object
Exp.COHORT.Inpath <- All.experiments.input.path[str_detect(All.experiments.input.path,EXPERIMENT)]
GRanges.object <- readRDS(Exp.COHORT.Inpath)

# adjusting for difference in experiment names due to the difference in loops in CalculatingScoreMatrix ()
if (EXPERIMENT=="RNA-Seq") {GRanges.object.names <- readRDS(str_replace(Exp.COHORT.Inpath,".GR.rds",".df.rds"))}


# removing info about genes and estimating size factors
scores <- as.matrix(mcols(GRanges.object)) 



# removing cell types that have  0 counts for all genes
  # identifying cells with 0 scores
    cells.with.zero.score <- which(apply(scores,2,function(x){((sum(x>0,na.rm = T)==0))}))
    
  # calculating size factors for all cells that have at least one count different from 0
    DF.object.esf <- estimateSizeFactorsForMatrix(scores[,-cells.with.zero.score])

    
  # creating vector of estimated size factors, in the case all 0, sizefactor is set to be 1  
    Size.Factors.Vector <- matrix(rep(1,ncol(scores)),ncol = ncol(scores))
    colnames(Size.Factors.Vector) <- colnames(scores)
    Size.Factors.Vector[,-cells.with.zero.score] <- DF.object.esf


print(Size.Factors.Vector)
# multiplying results of Calculating_Scores_Enhancer_Regions.R script with estimated size factors  
mcols(GRanges.object) <- cbind(as.character(GRanges.object$gene.id),as.character(GRanges.object$gene.name),as.vector(Size.Factors.Vector)*scores)


#saving data as GRanges object
saveRDS(GRanges.object,str_replace(Exp.COHORT.Inpath,".rds",".normalized.rds"))


# saving Data as data frame object to be easier to work with later on  
DataFrame.object <- DataFrame(GRanges.object)
DF <- data.frame(DataFrame.object)
DF <- DF[,colnames(DF)%in%colnames(mcols(GRanges.object))]
df.ss <-apply(DF,2,as.numeric)
if (EXPERIMENT=="RNA-Seq") {rownames(df.ss) <- row.names(GRanges.object.names)}

# saving data as DF  
saveRDS(df.ss,str_replace(Exp.COHORT.Inpath,".rds",".normalizedDF.rds"))




sink()
sink(type="message")

# sessionInfo()
