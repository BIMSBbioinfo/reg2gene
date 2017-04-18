
# IDEA: Calculate ScoreMatrixBin() for any genomic region file entered (GRanges object stores as .rds file - prefered) or .bed file
# Sofar 2 possible files were used as either gene or enhancer regions:
# GRanges_object.path <- "/data/akalin/Projects/AAkalin_Catalog_RI/Data/Enhancer_regions_def/Roadmap/Pooled_unique_stack_tiled_mnemonics_GRanges_124cells_6_EnhG_7_Enh_121116.rds"
# GRanges_object.path <- "/data/akalin/Base/Annotation/hg19/GENCODE/v24/gencode.v24lift37.basicAnnAndNoncodingGRanges.FilteringExLngth.ExonsReduced16_06_07.rds"
# INPUT EXAMPLES
# Cohort="Roadmap"
# Method="RNA-Seq"
# GRanges_object.path <- "/data/akalin/Projects/AAkalin_Catalog_RI/Data/Enhancer_regions_def/EnhRegions_tiled_and_resized.rds"
# GRanges_object.path <- "/data/akalin/Base/Annotation/hg19/GENCODE/v24/gencode.v24lift37.basicAnnAndNoncodingGRanges.FilteringExLngth.ExonsReduced16_06_07.rds"
# 






args <- commandArgs(TRUE)
Cohort <- as.character(args[1]) # Blueprint,CEMT,McGill or Roadmap
Method <- as.character(args[2]) # H3K27ac, "DNA Methylation","RNA-Seq", "DNase-Hypersensitivity"
GRanges_object.path <- as.character(args[3])



#####################################
# Setting Log files output
#####################################
sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/",Sys.Date(),Method,Cohort,"_ScoreMatrixBinScript_OUTPUT.Rout"),
          open="wt"))
print(paste("Tested cohort is ",Cohort," with experiment type ",Method, "using version of enhancers/genes",GRanges_object.path," performed on",Sys.Date()))
sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/",Sys.Date(),Method,Cohort,"_ScoreMatrixBinScript_ERROR.Rout"),
          open="wt"),type="message",append=TRUE)






#####################################
# Libraries and functions
#####################################

print(Sys.Date())
library(rtracklayer)
library(genomation)
library(stringr)
library(parallel)
library(GenomicRanges)
source("/data/akalin/Projects/AAkalin_Catalog_RI/Scripts/16_05_27_Consortium_data_Integration/Getting_Full_Paths_for_IndexFile_Function_2.R")
source("/data/akalin/Projects/AAkalin_Catalog_RI/Scripts/16_06_30_ScoreMatrixBin_Scripts/ScoreMatrixListAsGRangesMcols_ver2.R")



#Cohort="Roadmap"
#Method="H3K27ac"

#####################################
# Paths
#####################################

directory.path <- "/data/akalin/Base/"
out.path <- "/data/akalin/Projects/AAkalin_Catalog_RI/Data/ScoreMatrixBinResults_Roadmap/"
Index.files=readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Index_file_ALL_4_COHORT/16_06_27_All_4cohortsIndex_file.rds")    # detecting index files for all cohorts
All.possible.experiments <- names(table(as.character(Index.files[,"Experiment"])))
Gene_ranges_path <- "/data/akalin/Base/Annotation/hg19/GENCODE/v24/gencode.v24lift37.basicannotationAndNoncodingGRanges.Genes.160602.rds"


#####################################
# Tests
#####################################

if (!exists("Method")) {print("No Method/Experiment entered")}
if (!exists("Cohort")) {print("No Cohort entered")}
if (!exists("GRanges_object.path")) {print("No path to GRanges.rds object entered")}


# testing if Experiment and Cohort name was entered properly

if ((Cohort%in%c("Blueprint","CEMT","McGill","Roadmap"))==0) {stop("Wrong cohort name! Enter either: Blueprint,CEMT,McGill or Roadmap")}
if ((Method%in%All.possible.experiments)==0) {print("Wrong experiment name! Enter either:",
                                                    paste(All.possible.experiments,collapse=" "))}
# testing if cohort and Method directories exists, and I not they are created (Method nested within Cohort)
if (!dir.exists(paste0(out.path,Cohort))) {(dir.create(paste0(out.path,Cohort)))}
if (!dir.exists(paste0(out.path,Cohort,"/",str_replace(Method," ","")))) {(dir.create(paste0(out.path,Cohort,"/",str_replace(Method," ",""))))}






#####################################
# Import and set out directory
##################################### 

setwd(paste0(out.path,Cohort,"/",str_replace(Method," ","")))

# Import GRanges_object.path as either bed file or rds file
if (str_detect(GRanges_object.path,"rds")==T) { GRanges_object <- readRDS(GRanges_object.path)}
if (str_detect(GRanges_object.path,"bed")==T) { GRanges_object <- import.bed(GRanges_object.path,as="GRanges")}

#removing info from chrY and chrM 
GRanges_object <- GRanges_object[(!str_detect(chrom(GRanges_object),"chrY|M"))]

# getting full paths to Index file    
Index.files.full.paths <- Complete.paths.function(directory.path=directory.path,Index.file=Index.files,COHORT=Cohort,EXPERIMENT=Method)




#####################################
# Analysis
#####################################        





#####################################    
# Arranging log2Fold change if available for histone modifications    
###################################    

# if method/experiment tested is Histone modfication or DNaze then log2Fold change is by default set to be used
if (Method%in%c("DNA Methylation","RNA-Seq")!=T){
  # replace full paths with full paths for log2Fold .bw files
  Index.files.full.paths[,"bw.full.names"] <- str_replace(Index.files.full.paths[,"bw.full.names"],
                                                          Index.files.full.paths[,"bw"],Index.files.full.paths[,"Log2Fold"])
  # test if all log2Fold .bw files are actually present in directory with all the data and if not set to be NA
  All.possible.paths <- list.files(directory.path,recursive = T,full.names = T)
  Index.files.full.paths[!Index.files.full.paths[,"bw.full.names"]%in%All.possible.paths,"bw.full.names"] <- NA
  Index.files.full.paths <- Index.files.full.paths[complete.cases(Index.files.full.paths[,"bw.full.names"]),]
}



#####################################    
# Calculating coverage for every enhancer region (with arranging strand awareness for RNASeq data and c)
###################################    


# 1.st testing which regions fall out of the chromosome regions of .bw files
Not_covered_regions=which(end(GRanges_object) > seqlengths(BigWigFile(Index.files.full.paths[1,"bw.full.names"]))[as.character(seqnames(GRanges_object))])
# subsetting granges object for those that do not overlap
if (length(Not_covered_regions)>0) {GRanges_object <- GRanges_object[-Not_covered_regions]}


# 2. calculating scores
# splitting function for strand info
# if method is not RNASeq that there in no need to specially
if (Method!="RNA-Seq"){
  
  
  scores.exp=mclapply(Index.files.full.paths[,"bw.full.names"],function(x) { try(ScoreMatrixBin(x,windows = GRanges_object,bin.num = 1,type = "bigWig", is.noCovNA=T),silent = T)},mc.cores=round(length(Index.files.full.paths)/5))
  
  saveRDS(scores.exp,file=paste0("./",Sys.Date(),"_",Cohort,"_",str_replace(Method," ","_"),".rds"))
  
  
  # arranging everything as GRanges object
  GRanges_object <- ScoreMatrixList_as_GRanges_mcols(scores.exp,GRanges_object,Index.files.full.paths)
  
  saveRDS(GRanges_object,file=paste0("./",Sys.Date(),"_",Cohort,"_",str_replace(Method," ","_"),".GR.rds"))
  
  
}

# if method/experiment tested is RNASeq than strand needs to be taken into account
if (Method=="RNA-Seq"){ 
  # separate analysis for all strands, assumption is that it needs to be strictly defined from which strand gene 
  # is transcribed, if it is not known such genes will be eliminated
  
  
  # reading GRanges object with info about genes:
  GRanges_genes <- readRDS(Gene_ranges_path)
  GRanges_genes <- GRanges_genes[(!str_detect(chrom(GRanges_genes),"chrY|M"))]
  # elimination of GRanges object with unknown strand
  
  
  GRanges_object <- GRanges_object[!strand(GRanges_object)=="*"]
  # Separate strand analysis  
  Strand <- c("+","-")
  # Creating container list for results  
  Per.Strand.GrangesList <- list()
  
  for (i in 1:length(Strand)){
    
    GRanges_object.ss <- GRanges_object[strand(GRanges_object)==Strand[i]]
    #Subsetted index file based on strand
    Index.files.full.paths.ss <- Index.files.full.paths[Index.files.full.paths[,"Strand"]%in%c(Strand[i],"*"),]
    
    scores.exp <- mclapply(Index.files.full.paths.ss[,"bw.full.names"], function(x) { try(ScoreMatrixBin(x,windows = GRanges_object.ss,bin.num = 1,type = "bigWig", is.noCovNA=T),silent = T)},mc.cores=round(nrow(Index.files.full.paths)/5))
    
    
    GRanges_object.ss <- ScoreMatrixList_as_GRanges_mcols(scores.exp= scores.exp,GRanges_object = GRanges_object.ss,Index.files.full.paths = Index.files.full.paths.ss)
    Per.Strand.GrangesList[[Strand[i]]] <- GRanges_object.ss
    
  }
  
  
  GRanges_object <- c(Per.Strand.GrangesList[[1]],Per.Strand.GrangesList[[2]])
  # mcols(GRanges_object) <- mcols(GRanges_object)[!str_detect(names(mcols(GRanges_object)),".1")]
  saveRDS(GRanges_object,file=paste0("./",Sys.Date(),"_",Cohort,"_",str_replace(Method," ","_"),"exons.GR.rds"))
  
  
  
  # Calculating per gene RPKM
  # IDEA: split results according to gene and calculate RPKM (sum across exons(mean score across exon multiply with exon with))/gene width
  Container.RPKM=mclapply(split(GRanges_object,as.character(GRanges_object$sample)),function(x) {
    Cell.type=as.data.frame(mcols(x)[-1])
    scores.per.gene=(apply(Cell.type*width(x),2,sum,na.rm=T))/sum(width(x),na.rm=T) # sum across exons (mean value across exon*exon width)/ gene_length(sum_of_all_exon_lengths)
  }, mc.cores=10)
  
  
  #  creating data frame 
  Container.RPKM.df=do.call("rbind.data.frame",Container.RPKM)
  # saving the data
  colnames(Container.RPKM.df)=colnames(mcols(GRanges_object)[-1])
  
  
  GRanges_genes <- GRanges_genes[(GRanges_genes$gene.id%in%rownames(Container.RPKM.df))]
  GRanges_genes.ordered <- GRanges_genes[order(GRanges_genes$gene.id)]
  Container.RPKM.df <- Container.RPKM.df[order(rownames(Container.RPKM.df)),]
  # subsetting names that are extra in RNASeq file, and I have no idea why!!!!!
  Container.RPKM.df <- Container.RPKM.df[!(!rownames(Container.RPKM.df)%in%GRanges_genes.ordered$gene.id),]
  # testing if everything is okay with order of data    
  if (sum(rownames(Container.RPKM.df)==GRanges_genes.ordered$gene.id)!=length(GRanges_genes.ordered)) {stop}
  
  values(GRanges_genes.ordered) <- Container.RPKM.df
  saveRDS(Container.RPKM.df,file=paste0("./",Sys.Date(),"_",Cohort,"_",str_replace(Method," ","_"),"genes.df.rds"))
  saveRDS(GRanges_genes.ordered,file=paste0("./",Sys.Date(),"_",Cohort,"_",str_replace(Method," ","_"),"genes.GR.rds"))
  
}




sink()
sink(type="message")

