# Function that organizes results of linear models into df that are used as inputs for METAL program
# It runs METAL
# And stores all output of METAL program: results of MA and what metal actually reported
# Runs for every cohort, method, gene (hiearchically for any gene, get all available Methods, for given Method get available cohorts)


args <- commandArgs(TRUE)
MA.Type <- as.character(args[1]) # "P.value" or SE
in.path <- as.character(args[2]) # "/data/akalin/Projects/AAkalin_Catalog_RI/Results/GWAS_approach_results_Rdmp/"
out.path <- as.character(args[3]) # "/data/akalin/Projects/AAkalin_Catalog_RI/Results/MetaAnalysis/Rdmp/"







 #out.path <- "/data/akalin/Projects/AAkalin_Catalog_RI/Results/MetaAnalysis/Roadmap/"
 #in.path <- "/data/akalin/Projects/AAkalin_Catalog_RI/Results/GWAS_approach_results_Rdmp/"
 #MA.Type <- "P.value" # "P.value" or SE

 
 
 
# LIBRARIES
library(stringr)
library(parallel)
library(GenomicRanges)


 # ADJUSTING WHERE OUTPUT WILL BE SAVED 
sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/MA/",Sys.Date(),MA.Type,"MA.Rout"),
          open="wt"))
print(paste0("MA of all genes performed on ",Sys.Date()))
sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/MA/",Sys.Date(),MA.Type,"_MA2.Rout"),
          open="wt"),type="message",append=TRUE)

# adjusting outpath for MA type
out.path <- paste0(str_replace(out.path,"/$",""),"_",MA.Type,"/")
if (!dir.exists(paste0(out.path))) {(dir.create(paste0(out.path)))}








# defining variables
    Methods <- c("H3K27ac", "DNAMethylation","DNase-Hypersensitivity","H3K4me1")
    Cohorts <- c("Blueprint","CEMT","McGill","Roadmap")


  # adding info about sample sizes for cohorts    
  
    #N=matrix(c(48,32,29,56,51,39,NA,56,44,24,NA,44,29,39,14,32),nrow=4,ncol=4)
    # defining table
    N <- matrix(NA,nrow=4,ncol=4)
      colnames(N) <- Cohorts
      rownames(N) <- Methods
    
    # for one gene I have reported cell types that were used for the analysis
    One.gene.example.paths <- list.files(paste0(in.path,"statistics"),full.names = T,recursive=T)
    
    # extracting info about cell types used for the LM,and filling a table used in MA step
    for (i in 1:length(One.gene.example.paths)){
      
      TMP <- try(readRDS(One.gene.example.paths[i]))
    # chencking necessary because DNase not available for all cohorts
        if (class(TMP)!="try-error"){
         n.cell.types <- sum(!str_detect(TMP,"FILTERING|Error")) # extracting N of cell types
            # filling appropriate rows and columns 
              row.name <- Methods[str_detect(str_replace_all(One.gene.example.paths[i],"[[:punct:]]",""),str_replace_all(Methods,"[[:punct:]]",""))]
              col.name <- Cohorts[str_detect(One.gene.example.paths[i],Cohorts)]
      N[row.name,col.name] <- n.cell.types
      
      }
    
    }
    


# define it by this method to be reliable
  All.genes <- readRDS("/data/akalin/Base/Annotation/hg19/GENCODE/v24/gencode.v24lift37.basicAnnAndNoncodingGRanges.FilteringExLngth.ExonsReduced16_06_07.rds")
  All.genes <- as.character(unique(All.genes$sample))



# getting a list of all existing results
    All.GWAS.performed.list <- list.files(in.path,pattern="lm.output.categories.df",recursive = T, full.names = T)

    
    
##############################
# Function that performs MA
##############################
  MA_function=function(gene,All.GWAS.performed=All.GWAS.performed.list, matrix.N=N){
     #gene=All.genes[3];All.GWAS.performed=All.GWAS.performed.list; matrix.N=N;i=1;j=1
  
          Available.GWAS.per.gene <- All.GWAS.performed[str_detect(All.GWAS.performed,gene)]
        # identification of cohorts and methods that exist for given gene
          Cohorts.per.gene <- unique(unlist(str_extract_all(Available.GWAS.per.gene,"Roadmap|McGill|Blueprint|CEMT")))
          Methodes.per.gene <- unique(str_extract(Available.GWAS.per.gene,"H3K27ac|DNAMethylation|DNase-Hypersensitivity|H3K4me1"))
        
        
        GWAS.Results <- list()
        
        # running analysis per gene per method
                for (i in 1:length(Methodes.per.gene)){
          #i=1
          Method <- Methodes.per.gene[i]
          
               for (j in 1:length(Cohorts.per.gene)){
          #j=2  
           Cohort <- Cohorts.per.gene[j]
                    # add part where some gene might not be analyzed in different cohorts - due to filtering procedure
                    # test whether exists or not
                    gene.results <- try(readRDS(Available.GWAS.per.gene[str_detect(Available.GWAS.per.gene,paste0(Cohort,"/",Method))]),silent=T)
                               print(basename(Available.GWAS.per.gene[str_detect(Available.GWAS.per.gene,paste0(Cohort,"/",Method))]))
                    
                  # check if this gene actually exists
                if (class(gene.results)!="try-error"){
                      
                      tmp <- t((gene.results))
              # extracting info for MA
                    tmp.adjusting <- cbind(rownames(tmp),paste0(rownames(tmp),"_1"),paste0(rownames(tmp),"_2"),tmp[,c("Beta","SE.Beta","P.value")],matrix.N[Method,Cohort])
                    # adjusting names    
                    colnames(tmp.adjusting) <- c("EnhancerR","EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","P.value","N")
                    # adjusting paths
                    if (!dir.exists(paste0(out.path,Method))) {(dir.create(paste0(out.path,Method)))}
                        if (!dir.exists(paste0(out.path,Method,"/",str_replace(gene,"lm.output.categories.df|rds","")))) {(dir.create(paste0(out.path,Method,"/",str_replace(gene,"lm.output.categories.df|rds",""))))}
                    
                      
                    # saving per cohort per gene per method input for METAL program
                         setwd(paste0(out.path,Method,"/",str_replace(gene,"lm.output.categories.df|rds","")))
                         write.csv(tmp.adjusting,file=paste0("MA_",Cohort),quote=F, row.names = F)
          ##########
          # Running MA 
          # adjusting MA scrips
          RunningScript <- readLines("/data/akalin/Projects/AAkalin_Catalog_RI/Scripts/16_09_21_MetaAnalysis/MA_batch_Script.txt")

          # adjsuting Script that runs MA for existance or not of MA files
                Existing.MA <- list.files()
                    PROCESSING_LINES <- which(str_detect(RunningScript,"PROCESS")) # identfy Processed Cohorts
                    To_comment_lines <- !str_detect(RunningScript[PROCESSING_LINES],paste0(Existing.MA,collapse = "|"))
                # comment lines that were not analyzed - filtered due to some reason
                      RunningScript[PROCESSING_LINES[To_comment_lines]] <- paste0("#",RunningScript[PROCESSING_LINES[To_comment_lines]])
                # deciding whether you want P-values combined + direction of effect or effect size + SE
                      if (MA.Type=="SE"){
                        To_comment_lines <- str_detect(RunningScript,"SCHEME SAMPLESIZE")
                        RunningScript[To_comment_lines] <- paste0("#",RunningScript[To_comment_lines])}
                      if (MA.Type=="P.value"){
                        To_comment_lines <- str_detect(RunningScript,"SCHEME STDERR")
                        RunningScript[To_comment_lines] <- paste0("#",RunningScript[To_comment_lines])
                      }
                      
                # saving script      
                  writeLines(RunningScript,paste0(getwd(),"/MA_batch_Script.txt"))
                # run MA   
                  output.metal <- system("/gnu/store/8y1fic332nfqigfc77w4a3752fcy9smx-metal-2011-03-25/bin/metal MA_batch_Script.txt",intern=T)
              # store output of MA
                        writeLines(output.metal,"output.metal.txt")
                   }
        }
        
  }

}
##############################
# Running MA
##############################   
    
    
    

mclapply(All.genes,MA_function,mc.cores=20)



sink()
sink(type="message")