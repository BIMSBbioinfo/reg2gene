

setwd("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/")


# GET enhancer regions; GRanges object with 10 enhancer ranges

        Enhancer_regions <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Data/Enhancer_regions_def/Roadmap/Pooled_unique_stack_tiled_mnemonics_GRanges_124cells_6_EnhG_7_Enh_121116.rds")
        enhRegions <- GRanges(Enhancer_regions[1:10])
        
        save(enhRegions,file="pkg/inst/extdata/enhRegions.RData")
        

# GET 2 bw files with as an example.

         library(rtracklayer) 
        
        BW.files.LIVER <- "/data/akalin/Base/RoadmapEpigenomics/Experiment/H3K27ac/bigwig/E066-H3K27ac.fc.signal.bigwig"
      
           Chr10.bw = import(BW.files.LIVER, which=seqinfo(BigWigFile(BW.files.LIVER))["chr10"])  
               export.bw(Chr10.bw,"pkg/inst/extdata/E066-H3K27ac.chr10.fc.signal.bigwig")
       
               # E085 Fetal_Intestine_Small
        BW.files.FIS <- "/data/akalin/Base/RoadmapEpigenomics/Experiment/H3K27ac/bigwig/E085-H3K27ac.fc.signal.bigwig"
               
          Chr10.bw = import(BW.files.LIVER, which=seqinfo(BigWigFile(BW.files.LIVER))["chr10"])  
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
              
                            
                 save(TSS.location.gr,file="pkg/inst/extdata/TSS.RData")
