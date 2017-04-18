# Script that take .rds file of Enhancer activity and Gene Expression df, checks if there is duplicated entries on the level of Donor and TIssue,
# performs linear modelling and reports statistics
# 2 loops: for Blueprint data and Not blueprint data
# Not Blueprint data: paste combination of Donor and Cell type id and identify duplicates
# Select !duplicates and filter for unique entries
# Blueprint data: select row names, extract Donor ID (grep 1st "_", everything before is Donor ID)
# and compare it with predifined duplicated possible entries ("Tissue_description","Donor")
# performs filtering, association studies, Bonferroni correction
# as well for selected gene it saves cell types after they were filtered and lm was runned in out.path/statistics (easy to run script and calculate n of cell types per method per cohort)

# Cohort <- "Roadmap"
# Method <- "H3K4me1"
# inpath <- "/data/akalin/Projects/AAkalin_Catalog_RI/Data/RPKM_EnhActivity_Pairs_Rdmp/"
# out.path <- "/data/akalin/Projects/AAkalin_Catalog_RI/Results/GWAS_approach_results_Rdmp/log2exp/"
#gene.example <- "ENSG00000001497_LAS1L.rds"



args <- commandArgs(TRUE)
Cohort <- as.character(args[1]) # "Blueprint",CEMT,McGill or Roadmap
Method <- as.character(args[2]) # H3K27ac, "DNAMethylation","RNA-Seq", "DNase-Hypersensitivity","H3K4me1"
inpath <- as.character(args[3])# list.files("/data/akalin/Projects/AAkalin_Catalog_RI/Data/RPKM_EnhActivity_Pairs/",full.names=T)
out.path <- as.character(args[4])#"/data/akalin/Projects/AAkalin_Catalog_RI/Results/GWAS_approach_results_testing_071216/"
gene.example <- as.character(args[5])   #ENSG00000001497_LAS1L.rds"









# ------------------------------
# saving info about analysis
# ------------------------------

    sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/GWASA/",Sys.Date(),Method,Cohort,"_GWAS_APP_Roadmap.Rout"),
              open="wt"))
    print(paste("Tested cohort is ",Cohort," with experiment type ",Method," performed on",Sys.Date()))
    sink(file(paste0("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Log_files/GWASA",Sys.Date(),Method,Cohort,"_GWAS_APPe_Roadmap.Rout"),
              open="wt"),type="message",append=TRUE)




# ------------------------------
# packages and functions
# ------------------------------

    library(reshape2)
    library(GenomicRanges)
    library(stringr)
    library(parallel)
    library(stringr)


# function that reports if no overlap with enhancers identified ((ncol(INPUT) <= 2 ))
# if all cells have gene expression 0
# if matrix is very sparse >90% observations are 0
# additionaly writen within GWAS function is that if var of Enh region is ==0 elimninate that region - creates problems 
        filtering_funtion = function(INPUT){
  
  if(ncol(INPUT) <= 2 ){next}
  
  if(all(INPUT[,1]==0)){next}
  
  if(sum(INPUT[,1]==0)/nrow(INPUT) > 0.9){next}
  
  INPUT=INPUT[complete.cases(INPUT),,drop=FALSE]
  
  return(INPUT)
}

# function that takes gene expression and enhancer activity matrix and runs linear model function, and separately stores different statistics
       lm.statistics.function <- function(gene.expression.log.scale,x,result.list=list(),Promoter.name=Promoter.name){ 
  # function that runs linear model per gene for all nerby enhancers
  tmp=data.frame(cbind(gene.expression.log.scale,x))
  colnames(tmp)[1]=Promoter.name
  # run linear model gene-enhancer pair    
  test.statistics=summary.lm(lm(formula= tmp[,1] ~ tmp[,2]))
  # report separately statistics and run bonferoni adjustment
  Intercept=test.statistics$coefficients[1,"Estimate"]
  Beta=test.statistics$coefficients[2,"Estimate"]
  SE.Beta=test.statistics$coefficients[2,"Std. Error"]
  P.value=test.statistics$coefficients[2,"Pr(>|t|)"]
  t.value=test.statistics$coefficients[2,"t value"]
  
  lm.statistics=cbind(Intercept,Beta,SE.Beta,P.value,t.value)
  # reporting as a list and in 2nd step divide it to 2 df and list
  result.list[[1]] <- test.statistics; result.list[[2]] <- lm.statistics
  names(result.list) <- Promoter.name
  return(result.list)
  
}


# function that take .rds file of Enhancer activity and Gene Expression df, performs linear modelling and reports statistics
# 2 loops: for Blueprint data and Not blueprint data
# Not Blueprint data: paste combination of Donor and Cell type id and identify duplicates
# Select !duplicates and filter for unique entries
# Blueprint data: select row names, extract Donor ID (grep 1st "_", everything before is Donor ID)
# and compare it with predifined duplicated possible entries ("Tissue_description","Donor")
# performs filtering, association studies, Bonferroni correction
        Function_GWAS_APPROACH_PER_GENE_ENH = function(RPKM.files.single){
  
  # input data
  # RDS.files <- readRDS(RPKM.files[10])
  #  name.INPUT <- basename(RPKM.files[10])
  
  
  ###########
  # 1. DATA PREPARATION STEP 
  RDS.files <- readRDS(RPKM.files.single)
  INPUT <- apply(RDS.files,2,as.numeric)
  INPUT <- try(filtering_funtion(INPUT))
  
  if (class(INPUT)=="try-error"){ return(print("Did NOT PASS FILTERING"))}
  
  if (class(INPUT)!="try-error"){
    
    INPUT=INPUT[,!(apply(INPUT,2,var,na.rm=T)==0)] # eliminate enhancers with zero variability
    #naming
    name.INPUT=basename(RPKM.files.single)
    #naming
    Promoter.name=colnames(INPUT)[1] 
    rownames(INPUT) <- rownames(RDS.files)
    
    # function that eliminates Duplicated cell lines
    # identifying expected duplicated entries - by Tissue description and Donor columns
    # not necessary duplicated in INPUT data, since some duplicated entries were reoveved in prestep of creatind EnhProm .rds files (RNASeq-Exp .rds files)
    Duplicated.donors <- (Index.files.Cohort.Exp[duplicated(Index.files.Cohort.Exp[,c("Tissue_description","Donor")]),"Donor"])
    Duplicated.entries <- Index.files.Cohort.Exp[Index.files.Cohort.Exp[,"Donor"]%in%Duplicated.donors,]
    
    if (Cohort!="Blueprint"){
      
      which.unique <- duplicated(paste0(Duplicated.entries[,"Tissue_description"],Duplicated.entries[,"Donor"]))
      Duplicated_ID <- Duplicated.entries[which.unique,"Unique_id"]
      
      
      
      # removing duplicated entries
      INPUT <- INPUT[!rownames(INPUT)%in%Duplicated_ID,]
      
    }
    
    if (Cohort=="Blueprint"){
      # ordering
      INPUT <- INPUT[order(rownames(INPUT)),]
      Duplicated.entries <- Duplicated.entries[order(Duplicated.entries[,"Donor"]),]
      
      # Selecting Donor_ID from rownames of INPUT .rds files (if short name add 2nd part)
      Donors.per.INPUT <- unlist(lapply(str_split(rownames(INPUT),"_"),function(x) {
        # x <- str_split(rownames(INPUT),"_")[17]    
        tmp <- unlist(x)[1]
        tmp2  <- unlist(x)[2]
        if (nchar(tmp)<5) { tmp <-  paste(tmp,tmp2)}
        return(tmp)
      }))
      
      
      # all rows of INPUT data that overlap with expected duplicated entries
      Donors.per.INPUT.all <- Donors.per.INPUT[Donors.per.INPUT%in%Duplicated.entries[,"Donor"]]
      # however nor necessary all are duplicated rownames, thus identifying only true duplicated rows
      Donors.per.INPUT.duplicated <- unique(Donors.per.INPUT.all[duplicated(Donors.per.INPUT.all)])
      
      
      
      # Row that do not have any duplication
      Correct_row_names <- which(!Donors.per.INPUT%in%Donors.per.INPUT.duplicated)
      # Rows that are duplicated
      Duplicate_row_names <- which(Donors.per.INPUT%in%Donors.per.INPUT.duplicated)
      # Selecting of duplicated rows unique
      Corrected_Duplicated_Rows <- Duplicate_row_names[which(!duplicated(Donors.per.INPUT[Duplicate_row_names]))]
      # Poolin Correct and Corrected duplicated rows together
      Corrected_rows <- sort(c(Correct_row_names,Corrected_Duplicated_Rows))
      
      
      
      INPUT <- INPUT[Corrected_rows,]
      
      
      
    }
    
    # starting analysis
    
    #INPUT=readRDS(inpath[1])
    if (var(INPUT[,1])!=0){ #eliminate genes with zero var
      INPUT=INPUT[,!(apply(INPUT,2,var)==0)] # eliminate enhancers with zero variability
      Promoter.name=colnames(INPUT)[1]  
      gene.expression.log=(log2(INPUT[,1]+1)) # adjusting gene expression
      gene.expression.log.scale=scale(gene.expression.log) # scaling gene expression
      
      #gene.expression.log.scale=scale(INPUT[,1])
      Enhancer_scaling=scale(INPUT[,-1]) # scaling input data
      
      ###########
      # 2. LINEAR MODELING STEP            
      
      # running linear model function            
      lm.statistics.results=apply(Enhancer_scaling,2,lm.statistics.function,gene.expression.log.scale=gene.expression.log.scale,Promoter.name=Promoter.name,result.list=list())
      # separating results into lm.statistics and lm.statistics.subset
      lm.statistics.subset <- matrix(NA,ncol=5,nrow=length(names(lm.statistics.results)),dimnames = list( names(lm.statistics.results),
                                                                                                          c("Intercept","Beta","SE.Beta","P.value","t.value")))
      lm.statistics <- list()
      for (i in 1:length(lm.statistics.results)){
        lm.statistics[i] <- lm.statistics.results[[i]][1]
        #names(lm.statistics)
        lm.statistics.subset[i,] <- unlist(lm.statistics.results[[i]][2])
      } 
      
      
      
      # getting GRanges from colnames TO get info about distance
      tmp=str_split(colnames(INPUT[,-1]),pattern = "_")
      Enhancer.gr=stack(GRangesList((lapply(tmp,function(x){GRanges(seqnames = x[1], IRanges(start=as.numeric(x[2]),end=as.numeric(x[3])))}))))
      x=unlist(str_split(colnames(INPUT)[1],pattern = "_"))
      Promoter.gr=GRanges(seqnames = x[1], IRanges(start=as.numeric(x[2]),end=as.numeric(x[3])))
      # getting distance between promoter and enhancers - enhancer regions remained 2Mb large, so I needed to trim them
      options(warn=-1)
      Dist.Prom.Enh.gr=tryCatch(as.integer(abs(distance(promoters(Promoter.gr,start=1,end=1),Enhancer.gr))))
      options(warn=0)
      
      # adding info about distance to our results
      lm.statistics.subset=round(cbind(lm.statistics.subset,Dist.Prom.Enh.gr),4)
      lm.statistics.subset <- t(lm.statistics.subset)
      
      # saving info about cell types that were used in this analysis for one gene example
      if (gene.example==name.INPUT) { return(list(rownames(INPUT)))}
      ###########
      # 3. SAVING RESULTS and statistics from analyzed cell types     
      saveRDS(lm.statistics,paste0("lm.output.",name.INPUT))
      saveRDS(lm.statistics.subset,paste0("lm.output.categories.df",name.INPUT))
      
    }
  }
  
}



        
        
        
        
        
        

# ------------------------------
# reading and adjusting input parameters
# ------------------------------

        
        
        
# listing files from folder that contains RPKM per gene calculations
      inpath <- list.files(inpath,full.names=T) 

# creating directories where per gene per cohort per method datawill be stored
      if (!dir.exists(out.path)) {(dir.create(out.path))}
      if (!dir.exists(paste0(out.path,Cohort))) {(dir.create(paste0(out.path,Cohort)))}
      if (!dir.exists(paste0(out.path,Cohort,"/",str_replace(Method," ","")))) {(dir.create(paste0(out.path,Cohort,"/",str_replace(Method," ",""))))}
      # create folder where statistics about tested cell types used in analysis will be saved
      if (!dir.exists(paste0(out.path,"statistics/"))){dir.create(paste0(out.path,"statistics/"))}
  setwd(paste0(out.path,Cohort,"/",str_replace(Method," ",""),"/"))

# getting data paths
    Method.path <- paste0(inpath,"/",Method) 
    Inpath.RPKM.Directory <- Method.path[str_detect(Method.path,Cohort)]
    RPKM.files <- list.files(Inpath.RPKM.Directory,full.names = T, recursive = T) 


# Index ile required for additional filtering step - combination of Donor and Cell Type ID    
    Index.files <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/Index_file_ALL_4_COHORT/6_12_08_All_4cohortsIndex_file_ver2.rds")

# filterin index file for cohort of interest      
      if (Method=="DNAMethylation"){
        Index.files.Cohort.Exp <- Index.files[(Index.files[,"Cohort"]==Cohort)&(Index.files[,"Experiment"]=="DNA Methylation"),]   
}
      if (Method!="DNAMethylation"){
        Index.files.Cohort.Exp <- Index.files[(Index.files[,"Cohort"]==Cohort)&(Index.files[,"Experiment"]==Method),]   
}







# ------------------------------
# Analysis
# ------------------------------


  statistics.container <- unlist(mclapply(RPKM.files,function (x) try(Function_GWAS_APPROACH_PER_GENE_ENH(x)),mc.cores = 15))

# saving info about cell types that passed filtering and were used for linear modelling step
  saveRDS(statistics.container,paste0(out.path,"statistics/",gene.example,Cohort,Method,"cell_type.statistic.rds"))





    sink()
    sink(type="message")