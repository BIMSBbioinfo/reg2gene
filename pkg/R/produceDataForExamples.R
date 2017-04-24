

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
       
        
        

TSS <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Data/RPKM_EnhActivity_Pairs_Rdmp/Roadmap/H3K27ac/E")
