
Complete.paths.function <- function(directory.path, Index.file, COHORT, EXPERIMENT) {

  # IDEA: function that outputs Index file with full paths (otherwise only basename of .bw files is stored)
  # and select cohort and experiment of intrests, and subset table according to that
  # INPUT:
  # directory.path - path to the folder where all data is stored
  # Index.file - index file
  # COHORT - cohort of interest (possible entries: Blueprint,CEMT,McGill or Roadmap)
  # EXPERIMENT - experiment of interest
  # directory.path <- "/data/akalin/Base/"
  # Index.file <- Index.files
  # COHORT <- "CEMT"
  # EXPERIMENT <- "H3K27ac"
  
  
  
    #########################
    # cohort orientated analysis
     
      # selecting, testing cohort
        if (length(COHORT)!=0) { Cohort.ids <- which(Index.file[,"Cohort"]==COHORT) }
      # if no cohort indicated, analyze all. cohorts  
        if (length(COHORT)==0) { Cohort.ids <- 1:nrow(Index.file)
                                  print("All cohort analyzed")}
      # if wrong name of cohort is indicated, STOP
        if (length(Cohort.ids)==0) {stop("Wrong cohort name! Enter either: Blueprint,CEMT,McGill or Roadmap")}
      # subsetting by cohort name
        Index.ss <- (Index.file[Cohort.ids,])
      # ordering Index file by Unique_id
        Index.ss <- Index.ss[order(Index.ss[,"Unique_id"]),]
        
        
        
     #########################
     # experiment orientated analysis
        
      # subsetting by Experiment
        Experiment.overlapping.rows <- Index.ss[,"Experiment"]%in%c("Index.file",EXPERIMENT)
        
        Index.ss <- Index.ss[Experiment.overlapping.rows,]       
       # testing if the name of experiment was correctly entered         
        if (length(Experiment.overlapping.rows)==0) {print("Wrong experiment name! Enter either:",
                                                           paste(names(table(as.character(Index.file[,"Experiment"]))),collapse=" "))}
        
        # get full names og .bw files         
        
        # adjusting for problematic names - Roadmap E011,E012,E013: + to plus
           Index.ss[,"bw"] <-  str_replace(Index.ss[,"bw"],"\\+_","plus_")
           Index.ss[,"bw"] <-  str_replace(Index.ss[,"bw"],"\\-","minus_sign")
        
        bw.file.names.collapsed <- paste0(Index.ss[,"bw"],collapse="|") # setting names of bw files as a search string with which I will search full paths
        bw.full.names.paths <- (list.files(directory.path,recursive = T, full.names = T)) # path to folder of analyzed cohort
        
        # adjusting for problematic names - Roadmap E011,E012,E013: + to plus
        bw.full.names.paths <- str_replace(bw.full.names.paths,"\\+_","plus_")
        bw.full.names.paths <- str_replace(bw.full.names.paths,"\\-","minus_sign")
        
        bw.full.names <- bw.full.names.paths[str_detect(bw.full.names.paths,bw.file.names.collapsed)] # subsetting full paths according to the basenames
        
        # ordering .bw files to correspond to each other
        Ordering.file <- cbind(basename(bw.full.names),bw.full.names)
        colnames(Ordering.file) <- c("bw","bw.full.names")
        Index.ss <- merge(Index.ss,Ordering.file,by="bw")
        Index.ss <- apply(Index.ss,2,as.character)
        # reversing  problematic names - Roadmap E011,E012,E013: plus to +
        
              Index.ss[,"bw"] <- str_replace(Index.ss[,"bw"],"plus_","\\+_")
              Index.ss[,"bw.full.names"] <- str_replace(Index.ss[,"bw.full.names"],"plus_","\\+_")
              Index.ss[,"bw"] <- str_replace(Index.ss[,"bw"],"minus_sign","\\-")
              Index.ss[,"bw.full.names"] <- str_replace(Index.ss[,"bw.full.names"],"minus_sign","\\-")

        
 return(Index.ss)       
}

