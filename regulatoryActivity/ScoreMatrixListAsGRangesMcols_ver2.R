# function that takes output of ScoreMatrixBin function and corresponding GRanges_object, 
# and outputs the same GRanges object but with values of ScoreMatrixBin as mcols()

ScoreMatrixList_as_GRanges_mcols = function(scores.exp,GRanges_object,Index.files.full.paths) {
  
    ################################################  
    # function that takes output of ScoreMatrixBin function and corresponding GRanges_object, 
    # and outputs the same GRanges object but with values of ScoreMatrixBin as mcols()
   ################################################  
  
      # setting scores.ext to be mcols of input GRanges object        
      scores.exp.df=rbind.data.frame(scores.exp)
      # identying unique ids     
      Cell.types=Index.files.full.paths[,"Unique_id"]
      Method.available.bigwig.tracks=which(!is.na(Index.files.full.paths[,"bw.full.names"]))
      
      # adjusting for the fact that some cell type have more than one Experiment Result
                  # IDENTIFYING SUCH CELL TYPES, eg pz 284_MYELOID CELL
                    Duplicated.cell.types <- unique(Cell.types[Method.available.bigwig.tracks][duplicated(Cell.types[Method.available.bigwig.tracks])])
                
                        positions.where.adjusted.naming.required <- Cell.types[Method.available.bigwig.tracks]%in%Duplicated.cell.types
                # Adding the 1st 5 characters of name of bw file    
                        Cell.types[positions.where.adjusted.naming.required] <- paste0(Cell.types[positions.where.adjusted.naming.required],
                                                                                   str_extract(Index.files.full.paths[Method.available.bigwig.tracks][positions.where.adjusted.naming.required],".{5}"))
                    
      # Removing all \n to _, to adjust for possible problem caused by that
                        Cell.types[Method.available.bigwig.tracks] <- str_replace_all(Cell.types[Method.available.bigwig.tracks]," ","_")
                        
   # adding colnames to df  
       colnames(scores.exp.df)=Cell.types[Method.available.bigwig.tracks]# adding cell type names
              
      # taking abs value because genes from neg strand have neg values
       scores.exp.df=abs(scores.exp.df) 
      # adding data to GRanges_object and saving     
      mcols(GRanges_object)=cbind(mcols(GRanges_object),scores.exp.df) # adding info to GRanges object
      
      return(GRanges_object)
      
}




