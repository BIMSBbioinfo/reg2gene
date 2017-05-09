

setwd("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/")

# GET enhancer regions; GRanges object with 10 enhancer ranges, for FAS example such that results can overlap
       library(stringr)
       library(GenomicRanges)

        GenExpEnhancerReg <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Data/RPKM_EnhActivity_Pairs_Rdmp/Roadmap/H3K27ac/ENSG00000026103_FAS.rds")
        #GenExpEnhancerReg <- readRDS("~/ENSG00000026103_FAS.rds")
        
      # get GRanges for enhancer regions around FAS
        EnhancerReg <- colnames(GenExpEnhancerReg)[-1]
        EnhancerReg <- unlist(str_split(EnhancerReg,"_"))
        EnhancerReg.gr <-  GRanges(EnhancerReg[ seq(1,length(EnhancerReg),3)],IRanges(start=as.integer(EnhancerReg[seq(2,length(EnhancerReg),3)]),end=as.integer(EnhancerReg[seq(3,length(EnhancerReg),3)])))
       
      # sample 20 GRanges   
        regRegions <- sample(EnhancerReg.gr,20)
        
      # save an example  
        save(regRegions,file="pkg/inst/extdata/regRegions.RData")
        

# GET 2 bw files with as an example.

         library(rtracklayer) 
        
        BW.files.LIVER <- "/data/akalin/Base/RoadmapEpigenomics/Experiment/H3K27ac/bigwig/E066-H3K27ac.fc.signal.bigwig"
      
           Chr10.bw = import(BW.files.LIVER, which=seqinfo(BigWigFile(BW.files.LIVER))["chr10"])  
               export.bw(Chr10.bw,"pkg/inst/extdata/E066-H3K27ac.chr10.fc.signal.bigwig")
       
               # E085 Fetal_Intestine_Small
        BW.files.FIS <- "/data/akalin/Base/RoadmapEpigenomics/Experiment/H3K27ac/bigwig/E085-H3K27ac.fc.signal.bigwig"
               
          Chr10.bw = import(BW.files.FIS, which=seqinfo(BigWigFile(BW.files.FIS))["chr10"])  
               export.bw(Chr10.bw,"pkg/inst/extdata/E085-H3K27ac.chr10.fc.signal.bigwig")
       
        
        

 #     GET TSS example   - FAS gene on chr 10  
  #     Chromosome 10: 90,750,414-90,775,542 forward strand.
               # TSS a GRanges object that have the TSS location and associated
               #' gene expression values per cell type or condition as meta data. Each
               #' row should have a "name" and "name2" columns for unique id or name/symbol
               #' for the gene which the TSS is associated with. One could be Ensembl id and the
               #' other could be used for gene symbol.
               #' Other metadata column names should represent sample names/ids and should
               #' match the GRanges object provided via regActivity argument.
        
               library(stringr)
               library(GenomicRanges)
                                    #' 
                  TSS <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Data/RPKM_EnhActivity_Pairs_Rdmp/Roadmap/H3K27ac/ENSG00000026103_FAS.rds")
                       
                  # create GRanges
                        TSS.location <- colnames(TSS)[1]
                        TSS.location <- unlist(str_split(TSS.location,"_"))
                                 TSS.location.gr <-  GRanges(TSS.location[1],IRanges(start=as.integer(TSS.location[2]),end=as.integer(TSS.location[3])))
                         
                   # add mcols
                         name <-  "ENSG00000026103"
                         name2 <- "FAS"
                         gene.expression <- t(TSS[,1])
                         
                                 mcols(TSS.location.gr) <- cbind(name,name2,gene.expression)
                        TSS <- promoters(TSS.location.gr,1,1)  
                 save(TSS,file="pkg/inst/extdata/TSS.RData")


 # get regActivity example
                 
                 library(genomation)
                 library(GenomicRanges)
                 load("pkg/inst/extdata/regRegions.RData")
                 activitySignals <- c("pkg/inst/extdata/E085-H3K27ac.chr10.fc.signal.bigwig",
                                      "pkg/inst/extdata/E066-H3K27ac.chr10.fc.signal.bigwig")
                 
                 regActivity <- regActivity(regRegions,activitySignals)
                 
                 colnames(mcols(regActivity)) <- c("E085","EO66")
                 save(regActivity,file="pkg/inst/extdata/regActivity.RData")
                 

                 
                                  
# Get an example of  GeneExpSignals and LibStrand
        
            #  paths to RNA-Seq .bw files and corresponding strands 
                
                 library(stringr)
                 source("/data/akalin/Projects/AAkalin_Catalog_RI/Scripts/Data_Integration/16_05_27_Consortium_data_Integration/Getting_Full_Paths_for_IndexFile_Function_2.R")
                 
                 # detecting index files for all cohorts - master table with all .bw files integrated
                 Index.files=readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Index_file_ALL_4_COHORT/6_12_08_All_4cohortsIndex_file_ver2.rds")    
                 
                 
                 # extracting full paths to RNA-Seq .bw files from Roadmap
                 
                       RoadMap_RNASeqExample <- Complete.paths.function("/data/akalin/Base/RoadmapEpigenomics/Experiment/",Index.files,COHORT = "Roadmap","RNA-Seq")
                       
                       
                       RoadMap_RNASeqExample.paths <- RoadMap_RNASeqExample[,"bw.full.names"]
                       RoadMap_RNASeqExample.strands <- RoadMap_RNASeqExample[,"Strand"]
                       
                        # subsetting
                           GeneExpSignals <- RoadMap_RNASeqExample.paths[19:24]
                           LibStrand <- RoadMap_RNASeqExample.strands[19:24]
                           
                           save(GeneExpSignals,file="~/GeneExpSignals.RData")
                           save(LibStrand,file="~/LibStrand.RData")
                           
                 
                 #save(GeneExpSignals,file="pkg/inst/extdata/GeneExpSignals.RData")
                 #save(LibStrand,file="pkg/inst/extdata/LibStrand.RData")
                 

# Get an example of Exons GRanges object
                           
      library(GenomicRanges)             
                          
           Exons=readRDS("/data/akalin/Base/Annotation/hg19/GENCODE/v24/gencode.v24lift37.basicAnnAndNoncodingGRanges.FilteringExLngth.ExonsReduced16_06_07.rds")
           GRanges_genes <- readRDS("/data/akalin/Base/Annotation/hg19/GENCODE/v24/gencode.v24lift37.basicannotationAndNoncodingGRanges.Genes.160602.rds")
            
           GENES <- data.frame(GRanges_genes) 
                  GENES <- GENES[,c("seqnames","start","end","strand","gene.id","gene.name")]
           
              Genes.order <- match(as.character(Exons$sample),GENES$gene.id)
           
        mcols(Exons) <- cbind(mcols(Exons),GENES[Genes.order,])
           
           
            # add 2 genes from both strands as an example 
           Exons <- Exons[which(Exons$sample%in%c("ENSG00000026103","ENSG00000113119","ENSG00000025039","ENSG00000261469"))]
           
           
           
                         save(Exons,file="~/Exons.RData")
#                         save(Exons,file="pkg/inst/extdata/Exons.RData")
  
                         
                         
                         
# Get an example of  GeneExpSignals 
                         load("~/GeneExpSignals.RData")        

                         sampleNames <- str_extract(basename(GeneExpSignals),"E[0-9]{3}")
                             
                         
                         
                         save(sampleNames,file="~/sampleNames.RData")
#                                             save(sampleNames,file="pkg/inst/extdata/sampleNames.RData")
                         