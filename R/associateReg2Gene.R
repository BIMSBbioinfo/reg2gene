#' Associate regulatory regions to genes in CHUNKS
#'
#' The function runs \code{\link{associateReg2Gene}} function for chunks of 
#' genes to avoid reaching a memory quota in R
#'
#' @param input a GRangesList that contains regulatory activity and
#' gene expression results. It has to have specific columns, see the
#' output of \code{\link{associateReg2Gene}} function.
#' @param out.dir (character) Path to the directory where intermediate results
#' will be stored.
#' @param method one of "pearson","spearman","dcor","elasticnet",
#' "randomForest". Default: "pearson". 
#' @param scaleData if TRUE (default) the the values for gene expression and
#' regulatory activity will be scaled to have 0 mean and unit variance using
#' \code{\link[base]{scale}} function.
#' @param mc.cores number of cores to be used
#' @param B number of randomizations, default 1000. This procedure
#' is used to estimate
#' P-values for coeficients returned by different methods.
#' @param chunks Default: 1 - everything is analyzed together. If >1 then
#' regulatory regions are divided into N chunks (useful for >>N of genes/
#' regulatory regions to avoid reaching R memory quota). For each chunk of 
#' genes separate .rds object is saved into memory, and when analysis is 
#' finished for separate chunks, results are pooled together and saved. 
#' @param saveTag (character, defult: "ChunkAnalysis") Unique naming for created
#' .rds files which correspond to different gene chunks. This character is used 
#' to pull together all .rds files from one round of analysis.
#' @param remove.chunks (default T). To remove all chunks of .rds data created
#' while running a function
#' @return a GInteraction object containing every potential association and
#' between a regulatory region and TSS, and the estimated association statistic
#' , its P-value and Q-value.
#'
#' @details
#' The function implements \code{\link{associateReg2Gene}} function, but 
#' regulatory regions are dividied into predefined number of chunks; 
#' anaysis is performed separately for each chunk; and results are saved in 
#' predefined output directory. After analysis, all results are pulled
#' together, Q-value is calculated, and all intermediate files are removed.
#'
#' @author Inga Patarcic
#'
#' @import foreach 
#' @import stringr
#' @importFrom doMC registerDoMC
#' @importFrom ranger ranger
#' @import glmnet 
#' @importFrom qvalue qvalue
#' @importFrom fitdistrplus fitdist
#' @import InteractionSet
#' 
.associateReg2GeneChunks <- function(input,
                                    out.dir,
                                    method="pearson",
                                    saveTag="ChunkAnalysis",
                                    scaleData=TRUE,
                                    mc.cores=1,
                                    B=1000,
                                    chunks=1,
                                    remove.chunks=T){
  
  # save warnings, output,etc. 
  sink(paste0(out.dir,"log.txt"), append=TRUE, split=TRUE)
  
  # splitting in gene chunks
  Split.factor <- split(1:length(input), sort(1:length(input)%%chunks))
  
  
  # run associateReg2Gene for gene chunks and save them separately
  lapply(Split.factor,function(x){
    print(x[1])
    
    # adjusting for the fact that code breakes from time to time
    
    counter=1
      while(counter<=10){
        Res_associateReg2Gene <- try(associateReg2Gene(input[x],
                                                       method=method,
                                                       B=B,
                                                       scaleData=scaleData,
                                                       mc.cores=mc.cores),silent=TRUE)
        if(!inherits(Res_associateReg2Gene, 'try-error')) break
        if(inherits(Res_associateReg2Gene, 'try-error')) {counter=counter+1}
      }
    
    # report error even after 10 iterations of trying
    if (inherits(Res_associateReg2Gene, 'try-error')&counter>10){print(paste0(out.dir,
                                                          saveTag,x[1],
                                                    ".rds was not successful"))}
    
    saveRDS(Res_associateReg2Gene,paste0(out.dir,saveTag,x[1],".rds"))})
  
  
  # input all associateReg2Gene gene chunks in correct order
  GeneChunks <- sort(list.files(out.dir,pattern = saveTag,full.names = T,
                                recursive=T))
  GeneOrder <- order(as.numeric(str_replace(str_extract(GeneChunks,
                                                  "[0-9]*.rds"),".rds",""))) 
  
  GeneEnhancerResults <- do.call(c,lapply(GeneChunks[GeneOrder], 
                                          function(x) readRDS(x)))
  
  if (typeof(GeneEnhancerResults)=="list"){
    GeneEnhancerResults <- do.call(c,GeneEnhancerResults)}
  
  # calculating appropriate Qvalue for all reported Pvalues
  # assumption Pvalue is always stored in the 6th column
  
  QValue <- qvalue(mcols(GeneEnhancerResults)[,6])$qvalues
  mcols(GeneEnhancerResults)[,7] <- QValue
  
  # remove all input files
  if (remove.chunks==T){lapply(GeneChunks,function(x){system(paste0("rm ",x))})}
  
  return(GeneEnhancerResults)
  
  sink()
  
  
} 

#' Associate regulatory regions to genes
#'
#' The function associates/links regulatory regions to genes based
#' on regulatory activity and gene expression relationship.
#'
#' @param input a GRangesList that contains regulatory activity and
#' gene expression results. It has to have specific columns, see the example 
#' below or the output of \code{\link{regActivityAroundTSS}} function.
#' @param method one of "pearson","spearman","dcor","elasticnet",
#' "randomForest". Default: "pearson". See Details for more.
#' @param scaleData if TRUE (default) the the values for gene expression and
#' regulatory activity will be scaled to have 0 mean and unit variance using
#' \code{\link[base]{scale}} function.
#' @param mc.cores number of cores to be used
#' @param B number of randomizations, default 1000. This procedure
#' is used to estimate P-values for coeficients returned by different methods.
#' @param asGInteractions if TRUE (default) results are reported as 
#' \code{\link[InteractionSet]{GInteractions}}. Otherwise, a GRanges object with 
#' regulatory region coordinates and an additional column "reg" containing gene
#' GRanges is reported
#' @param ... further arguments to methods, not implemented yet
#' @return a GInteraction object containing every potential association and
#' between a regulatory region and TSS, and the estimated association statistic
#' , its P-value and Q-value.
#'
#'
#' @details
#' The function implements four methods to associate regulatory
#' activity to genes: correlation (pearson or spearman), distance
#' correlation ("dcor"), elasticnet and random forests. Simple linear regression on
#' scaled data is equivalent to Pearson correlation. Distance correlation is
#' a metric for statistical dependence between two variables. It can capture
#' non-linear relationships different from correlation.
#' Elasticnet and randomForests are methods that
#' can capture additivity between regulatory activities and their relationship
#' to gene expression. For all the methods, P-values are estimated by shuffling
#' the gene expression vector \code{B} times and getting a null distribution
#' for the estimated statistic, and comparing the original statistic to the null
#' distribution. Estimated statistics are
#' correlation coefficients for correlation and distance correlation methods,
#' Regression coefficients for elasticnet and gini importance scores for random
#' forests.
#'
#' @author Altuna Akalin
#'
#' @import foreach 
#' @import stringr
#' @importFrom doMC registerDoMC
#' @importFrom ranger ranger
#' @import glmnet
#' @importFrom qvalue qvalue
#' @importFrom fitdistrplus fitdist
#' 
#' @examples 
#' ###############################
#' #STEP 1.  Getting random and predefined .8 correlation
#'  
#'  require(GenomicRanges)
#'  require(doMC)
#'  require(glmnet)
#'  require(foreach)
#'  require(stringr)
#'  require(qvalue)
#'  
#'  ####################################
#'  # INPUT1: GRanges
#'  # create an example: 2 vectors with correlation 0.8
#'  
#'    x <- c(2.000346,2.166255,0.7372374,0.9380581,2.423209, 
#'         2.599857,4.216959,2.589133,1.848172,3.039659)
#'    y <- c(2.866875,2.817145,2.1434456,2.9039771,3.819091,5.009990,
#'         5.048476,2.884551,2.780067,4.053136)
#'      corrM <- rbind(x,y)
#'  
#'  # define Granges object:
#'    gr0 <- GRanges(seqnames=rep("chr1",2),IRanges(1:2,3:4))
#'     
#'       GeneInfo <- as.data.frame(matrix(rep(c("gene","regulatory"),each=3),
#'                  ncol = 3,byrow = TRUE),stringsAsFactors=FALSE)
#'
#'         colnames(GeneInfo) <- c("featureType","name","name2")
#'
#'    mcols(gr0) <- DataFrame(cbind(GeneInfo,corrM))
#'  
#'  
#'  ####################################
#'  # OUTPUT: associateReg2Gene
#'  
#'  
#'     associateReg2Gene(gr0,mc.cores = 1)
#'     associateReg2Gene(gr0,mc.cores = 1,asGInteractions=FALSE) # report as GRanges
#'     associateReg2Gene(GRangesList(gr0),mc.cores = 1)
#'      
#'  # change N of resampling rounds   
#'     associateReg2Gene(GRangesList(gr0),mc.cores = 1,B=100) 
#'  # elastic net should return all NA values because 
#'  # only one predictor variable(x) is used to predict y  
#'   
#'     associateReg2Gene(GRangesList(gr0),mc.cores = 1,method="elasticnet") 
#' 
#' @export
associateReg2Gene<-function(input,
                            method="pearson",
                            scaleData=TRUE,
                            mc.cores=1,
                            B=1000,
                            asGInteractions=T,
                            ...){
 
   # drop NULL genes in the list
  nulls=which(sapply(input,is.null))
  if(length(nulls)>0){
      input=input[-c(nulls)]
  }

  # set up cores for multicore shit
  if(mc.cores>2){
    registerDoMC(mc.cores)
  }

  # decide which prediction/association method you want to use
  
  if(method=="pearson"){
    model.func=quote(corResample(mat,scaleData=scaleData,method="pearson",col=1,B=B))
  }else if(method == "spearman"){
    model.func=quote( corResample(mat,scaleData=scaleData,method="spearman",col=1,B=B) )
  }else if(method == "dcor"){
    model.func=quote( dcorResample(mat,scaleData=scaleData,col=1,B=B) )
  }else if(method == "elasticnet"){
    model.func=quote( glmnetResample(mat,scaleData=scaleData,col=1,B=B) )
  }else if(method == "randomForest"){
    model.func=quote( rfResample(mat,scaleData=scaleData,col=1,B=B) )
  }
 

  
  # function runs for GRangesList and GRanges
  if ((class(input)!="GRangesList")&(class(input)!="GRanges")){
    
    stop("Error in INPUT: neither GRanges nor GRangesList")
    
  }
  
  if (class(input)=="GRangesList") {
        
          res=foreach(gr=input) %dopar% {
            
            
        if (detectGeneEnhPresence(gr)) { 
          # check whether gene~enhancer info is provided in featureType column
          
            # def min output:
            x  <- data.frame(matrix(NA,nrow=length(gr$name[-1]),ncol=4,byrow=T,
                       dimnames=list(gr$name[-1],c("n","coefs","pval","pval2"))))
            
            mat <- getGeneEnhScoresDF(gr) # get DF of gr metadata (+remove NA)
              
            if (is.matrix(mat)) {
              
                if (ncol(mat)>0) {
              # ensures that function works after stepwise filtering:
              # covers cases when all enhancers are filtered out +/- gene filter
                 if (!all(zeroVar(mat)|manyZeros(mat))){ 
                    # filtering for low gene/enh variability and filtering prob
                        ModelRes <- eval(model.func)
                        
                        x <- perGeneModelling(x,mat,ModelRes)
                  }
                }
              }  
          
            return(x)
            
        }
            
          }
            # combine perGene stats with Granges
                
                comb.res=do.call("rbind.data.frame",res)[,1:4]
                
                    rownames(comb.res)=NULL # needed for GRanges DataFrame
         
          input <- input[as.logical(sapply(res,length))] # remove GE with wrong 
                    #featureType                     
          
      }
      
      if (class(input)=="GRanges") { 
        
        if (detectGeneEnhPresence(input)) { 
          # check whether gene~enhancer info is provided in featureType column
        
        comb.res  <- data.frame(matrix(NA,nrow=length(input$name[-1]),ncol=4,byrow=T,
                    dimnames=list(input$name[-1],c("n","coefs","pval","pval2"))))
        
        mat <- getGeneEnhScoresDF(input) # get DF of gr metadata (+remove NA)
        
        if (is.matrix(mat)) {
          
          if (ncol(mat)>0) {
          
              if (!all(zeroVar(mat)|manyZeros(mat))){ 
                # filtering for low gene/enh variability and filtering prob
                ModelRes <- eval(model.func)
                
                comb.res <- perGeneModelling(comb.res,mat,ModelRes)
              
                  }
               
              }
            
          }       
      
        
        }
      }
       
    if (exists("comb.res")) {
      
      if (as.logical(nrow(comb.res))) { # not existing when featureType entered wrongly
    
       comb.res <- qvaluCal(comb.res)  # add qvalue calculations
     
       comb.res <- comb.res[,c("n","coefs","pval","qval")] # pval2&qval2 removed
       
  # a function to create GRanges object for GRangesList
  # later p-values, effect sizes etc will appended to this object
     
      
       if (asGInteractions==T){gr2 <- grlist2GI(input)}
      
       if (asGInteractions!=T){gr2 <- grlist2gr(input)}
      
      mcols(gr2)=cbind(mcols(gr2),DataFrame(comb.res))

  # return result
  return(gr2)
    
      }
  
      if (!as.logical(nrow(comb.res))) {
        
        message("Error in INPUT,featureType should contain 1 gene + min 1 enh")}
      
  }else {
     
      message("Error in INPUT,featureType should contain 1 gene + min 1 enh")}
  
}






#########-------------------------------------------------###########
########### utility functions for association prediction ############

# this part contains utility functions that will be needed for all
# prediction functions

#' convert GRangesList to GRanges with associated regulatory regions (internal)
#'
#' The function converts output of \code{\link{regActivityAroundTSS}} (a
#' GRangesList) to a GRanges object by rearranging the rows and columns so
#' that output GRanges object mainly contains regulatory region coordinates 
#' and an additional column "reg" which contains the gene info
#'
#' @param grlist a GRangesList output from \code{\link{regActivityAroundTSS}}
#'
#' @keywords internal
#' @author Altuna Akalin, Inga Patarcic
grlist2gr<-function(grlist){
  
  if (class(grlist)=="GRangesList") {
    
    return(unlist(endoapply(grlist,grReorg)))
  
    }
  
  if (class(grlist)=="GRanges") { 
    
    return(grReorg(grlist))
  
  }
}

#' Help for grlist2gr()
#' @keywords internal
#' @author Inga Patarcic
grReorg <- function(x){
  
  reg=granges(x[x$featureType != "gene",])
  gene=x[x$featureType == "gene",c("name" ,"name2")]
  
  grpairs=rep(gene,length(reg))
  reg$reg=grpairs
  reg$name=grpairs$name
  reg$name2=grpairs$name2
  
  
  reg
}



#' convert GRangesList to GInteractions 
#' 
#' The function converts output of \code{\link{regActivityAroundTSS}} (a
#' GRangesList) to a GInteractions object by rearranging the rows and columns so
#' that output GInteractions object contains regulatory region coordinates 
#' and an additional column "reg" which contains the gene info
#'
#' @param grlist a GRangesList output from \code{\link{regActivityAroundTSS}}
#'
#' @keywords internal
#' @author Altuna Akalin, Inga Patarcic

grlist2GI<-function(grlist){
  
  if (class(grlist)=="GRangesList") {
    
    gr <- unlist(endoapply(grlist,grReorg))
    
  }
  
  if (class(grlist)=="GRanges") { 
    
    gr <- grReorg(grlist)
    
  }

  # arranging as GInteractions object
  
  gi <- GInteractions(gr,gr$reg)
  
  mcols(gi) <- mcols(gi)[c("anchor1.name","anchor1.name2")]
  
  names(mcols(gi)) <- c("name","name2")
  
  return(gi)
  }





#' Estimate P-values from resampling statistics using Gamma distribution
#'
#' The function estimates gamma distribution based P-values from the null
#' distribution obtained from resampling, a gamma distribution is fit
#' to the null distribution and p-values are calculated based on that fitted
#' distribution. Details: stats::pgamma
#'
#' @param vals a vector or matrix of column vectors. Each column corresponds
#' to coefficients obtained from resampling fit
#' @param orig a single value or a vector of original coefficients obtained from the
#' original matrix
#' @param abs if TRUE, use absolute values of coefficients to fit a distribution
#' @param add an numeric value defaults to 0. if more than 0, add that value to
#' orig and vals. This way gamma distribution can model null distributions that
#' have negative values.
#'
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pgamma
#' 
#' @keywords internal
estimateGammaPval<-function(vals,orig,abs=TRUE,add=0){

  orig2=(orig) # keep original values

  # if needed take absolute values
  if(abs){
    vals=abs(vals)
    orig2=abs(orig2)
  }
  if(add){ # adding numbers might be needed, gamma can only modep positive values
    vals=vals+add
    orig=orig+add
  }

  if(class(vals) != "matrix"){

    if(any(is.na(vals))){
      result=c(coefs=NA,pval=NA,p2=NA)

    }else{
      param <-  fitdistrplus::fitdist(vals,"gamma",method="mme") # fit gamma
      # p-val based on gamma
      pval=1-pgamma(orig2,param$estimate[1],param$estimate[2])

      # if 0 p-val returned set it to minimum possible or another option is
      # 2.2e-16
      if(pval==0){pval=1-pgamma(max(vals),param$estimate[1],param$estimate[2])}

      # result vector
      result=c(coefs=orig,pval=pval,p2=sum(vals >= orig2)/length(vals))}

  }else{

    # same as above it is just when the input is a matrix
    result=mapply(function(x,orig2,orig){
      if(any(is.na(x)) | all(x==0)){
        c(coefs=NA,pval=NA,p2=NA)

      }else{
        param <-  fitdistrplus::fitdist(x,"gamma",method="mme")

        pval=1-pgamma(orig2,param$estimate[1],param$estimate[2])
        if(pval==0){pval=1-pgamma(max(x),param$estimate[1],param$estimate[2])}
        c(coefs=orig,pval=pval,p2=sum(x >= orig2)/length(x) )
      }
    },split(t(vals),1:ncol(vals)),orig2,orig)

  }
}

#' decide if the input colum matrix has zero variation
#'
#' The input column matrix represent gene expression and regulatory region
#' activities. If it has zero variation for gene expression
#'
#' @param mat a column matrix that contains gene expression values and
#'            regulatory region activities.
#' @param col  column number of response variables
#'
#' @keywords internal
#' @author Altuna Akalin
zeroVar<-function(mat,col=1){

  # flag 0 variation in gene expression
  if(sd(mat[,col],na.rm=T)==0){
    return(TRUE)
  }
  FALSE
}

#' decide if the input colum from a matrix has many zeros
#'
#' The input column matrix represent gene expression and regulatory region
#' activities. decides if the genes
#'
#' @param mat a column matrix that contains gene expression values and
#'            regulatory region activities.
#' @param col  column number of response variables
#'
#' @keywords internal
#' @author Altuna Akalin
manyZeros<-function(mat,col=1){
  # if many values are zero cross-validation or subsampling based techniques
  # will fail
  if( sum(mat[,col]==0,na.rm=T)/nrow(mat) > 0.9 ){
    return(TRUE)
  }
  FALSE
}

#' Stepwise complete cases method
#'
#' Removes NA but in e order: 1st cells that have NA for at least
#' 3/4 of enhancer regions, then the same thing is applied for cell types
#' then individual NAs are removed
#'
#' @param mat a column matrix that contains gene expression values and
#'            regulatory region activities.
#'
#' @keywords internal
#' @author Inga Patarcic
stepwise.complete.cases <- function(x){
  
  similarity=0.75
  
  remove.cc <- function(x,similarity){
    n.na <- sum(is.na(x))
    if (n.na>=similarity*length(x)){return(F)}
    if (n.na<similarity*length(x)){return(T)}
  }
  
  # remove problematic cells, then enhancers then remove remaining NAs 
  x <- x[apply(x,1,remove.cc,similarity=similarity),]
  x <- x[,apply(x,2,remove.cc,similarity=similarity)]
  
  if (is.matrix(x)){x <- x[,apply(x,2,remove.cc,similarity=0.00001)]}
  
  return(x)
  
}


#' GRanges contains gene+enhancer info
#' 
#' Function that tests each GRanges object in GRangesList that it contains 
#' necessary info to run associateReg2Gene(): info about gene expression and
#' enhancer activity
#'
#' @param x GRanges object - which contains info about gene expression and
#' enhancer activity
#'
#' @keywords internal
#' @author Inga Patarcic
detectGeneEnhPresence <- function(x){

  TestGene <- sum(x$featureType%in%"gene")==1
  TestEnh <-sum(x$featureType%in%"regulatory")>=1
 
  return(TestGene&TestEnh)
   
}  
  
  




#' internal qvalue calculation function
#'
#' private function calculates qvalue for 2 categories of pvalue
#' @importFrom qvalue qvalue
#' @keywords internal
qvaluCal <- function(comb.res){
  
  # add qvalues for pval1&2 - adjusted for filtered results
  
  qvalR <- matrix(NA,ncol=2,nrow=nrow(comb.res),
                  dimnames =  list(1:nrow(comb.res),c("qval", "qval2")))
  comb.res=cbind(comb.res,qvalR)
  
  # add qval for pval1&2  
  test <- !all(is.na(comb.res[,3])) # test if any NA present
  # adjust for stimated pi0 <= 0 
  if (test){qval <- try(qvalue::qvalue(comb.res[,3])$qvalues,silent=T)
  
  if(!inherits(qval, 'try-error')){comb.res$qval <- qval}
  
  }
  
  test <- !all(is.na(comb.res[,4]))  # test if any NA present
  if (test){qval <- try(qvalue::qvalue(comb.res[,4])$qvalues,silent=T)
  
  if(!inherits(qval, 'try-error')){comb.res$qval2 <- qval}
  
  }
  return(comb.res)
  
}


#' extract gene-enhancer scores from GRanges object
#'
#' Extract gene-enhancer scores from metadata of GRanges object and filters NA
#' values using stepwise.complete.cases()
#'
#' @param gr perGene GRanges object
#'
#' @keywords internal
#' @author Inga Patarcic
getGeneEnhScoresDF <- function(gr) {
  
  mat=t(as.matrix(mcols(gr)[,-which(names(mcols(gr)) %in% c("name","name2",
                                                            "featureType"))]))
  colnames(mat) <-  gr$name
  
  # se complete.cases (remove cells, remove enhancers, remove random NA)  
  mat <- stepwise.complete.cases(mat)
  
  return(mat)
  
}      


#' Function returns output of called modelling procedure indivudually for 
#' each GRanges object
#'
#' It extracts gene-enhancer scores from metadata of GRanges object and filters
#' NA values using stepwise.complete.cases() embedded in getGeneEnhScoresDF()
#' Performs FILTERING: drop gene and NA is returned due to low or zero variation
#' in gene expression or when no enhancers are remained after 
#' stepwise.complete.cases() 
#'
#' @param gr perGene GRanges object
#'
#' @keywords internal
#' @author Inga Patarcic
perGeneModelling <- function(x,
                             mat,
                             ModelRes){
  
  nEnhancers <- ncol(mat)-1 # how many regulatory regions
  
  if ((ncol(mat)>=2)&
      (nEnhancers==nrow(x))){
    
    x <- data.frame(cbind(n=nrow(mat),t(ModelRes)))
    
    
  }else if ((ncol(mat)>=2)&
            (nEnhancers!=nrow(x))){
    
    # adjusting for the case where just some enhancers are filtered
    x[colnames(mat)[-1],] <- cbind(n=nrow(mat), t(ModelRes))
    
  }
  
  colnames(x) <- c("n","coefs","pval","pval2")
  
  return(x)
  
}


#### END OF utility functions for association prediction ####
#########------------------------------------------------############



#########-------------------------------------------------###########
########### association prediction functions

#' Correlation between a col matrix and a vector
#'
#' A convenience function that computes the the correlation of a given vector
#' with all the other column vectors of a matrix
#' @keywords internal
corMat<-function(y,mat,method="pearson"){

  if (!is.vector(mat)){
  return(apply(mat,2,function(x,y,method){cor(x,y,method=method)},
        y=y,method=method))}
  
  if (is.vector(mat)){return(cor(mat,y,method=method))}
  
}

#' Distance correlation between a col matrix and a vector
#'
#' a convenience function that computes the the distance correlation of a given
#' vector (y) with all the other columns of a matrix (mat)
#'
#' @keywords internal
#' @importFrom energy dcor
dcorMat<-function(y,mat){
  #require(energy)
  if (!is.vector(mat)){return(apply(mat,2,
                                    function(x,y){energy::dcor(x,y)},y=y))}
  
  if (is.vector(mat)){return(dcor(mat,y))}
}

#' Correlation with resampling p-values
#'
#' private function that calculate correlation of given
#' column matrix with a column and assigns significance to
#' correlation coefficents based on resampling
#'
#' @param mat matrix of response and predictor variables
#' @param col column number of response variable
#' @param method pearson or spearman
#' @param B number of random samples (shuffling) of given column
#' @keywords internal
#'
corResample<-function(mat,
                      scaleData=scaleData,
                      method="pearson",
                      col=1,
                      B=B){


  # scaleData if necessary
  if(scaleData){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }

  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))

  # original coefs
  # if only 1 enhancer available
  # if more than one enhancer available
  orig=corMat(y=mat[,col],mat=mat[,-col],method)

  #coefs from resampling (adjusted for 1 or more enhancers)
  Nenh <- ncol(mat[,-col])
    if (is.null(Nenh)){coefs=matrix(0.0,ncol=1,nrow=(B))}
    if (!is.null(Nenh)){coefs=matrix(0.0,ncol=Nenh,nrow=(B))}

  # calculate coeff for resampled Ys
  for(i in 1:B){
    #cat(i,"\n")
    coefs[i,] = corMat(Ys[[i]],mat=mat[,-col],method)
  }

  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}


#' Distance correlation with resampling based p-values
#'
#' private function that calculates distance correlation of given
#' column matrix with a vector and assigns significance to
#' correlation coefficents based on resampling
#' @param mat matrix of response and predictor variables
#' @param col column number of response variable
#' @param B number of random samples (shuffling) of given column
#' @keywords internal
dcorResample<-function(mat,scaleData=scaleData,col=1,B=B){

  if(scaleData){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }
  
  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))

  # original coeff
  orig=dcorMat(y=mat[,col],mat=mat[,-col])

  #coefs (adjusted for 1 or more enhancers)
  Nenh <- ncol(mat[,-col])
  if (is.null(Nenh)){coefs=matrix(0.0,ncol=1,nrow=(B))}
  if (!is.null(Nenh)){coefs=matrix(0.0,ncol=Nenh,nrow=(B))}
  


  # calculate coeff for resampled Ys
  for(i in 1:B){
    #cat(i,"\n")
    coefs[i,] = dcorMat(Ys[[i]],mat=mat[,-col])
  }

  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}



#' internal function that uses glmnet for prediction
#'
#' private function that uses glmnet to predict fiven column from the rest
#' and assigns significance to coefficents based on resampling
#' @param mat matrix of response and predictor variables
#' @param col column number of response variable
#' @param B number of random samples (shuffling) of response variables
#' @keywords internal
#' @import glmnet
glmnetResample<-function(mat,scaleData=scaleData,col=1,B=B){
  #require(glmnet)

  if (ncol(mat)<=2){
    return(matrix(NA,ncol=ncol(mat)-1,nrow=3,
                  dimnames=list(c("coefs","pval","p2"),1:(ncol(mat)-1) ))
    )
  }
  
  
  if(scaleData){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }

  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))

  # calculating coef and lambda.min
      coefs <- matrix(NA,ncol=ncol(mat[,-col]),nrow=(B))
      orig <- rep(NA,ncol(mat)-1)
            names(orig) <- colnames(mat)[-col] 
  
  # original coefs
  tryCatch({
        mod<- cv.glmnet(x = mat[,-col], y = mat[,col],
                        standardize=FALSE,
                        nfolds=5,alpha=0.5)
      
        orig=coef.cv.glmnet(mod,s="lambda.min")[-1,1]
      
        #coefs from resampling
        coefs=matrix(0.0,ncol=ncol(mat[,-col]),nrow=(B))
        
  },error=function(coefs){coefs <- coefs})

  # calculate coeff for resampled Ys
  # catching errors if in resampling all 0s selected 
  
  tryCatch({
        for(i in 1:B){
      
          mod<- cv.glmnet(x = mat[,-col], y = Ys[[i]],
                          standardize=FALSE,
                          nfolds=5,alpha=0.5)
      
          coefs[i,]=glmnet::coef.cv.glmnet(mod,s="lambda.min")[-1,1]
        }
  },error=function(coefs){
    
          coefs <- matrix(NA,ncol=ncol(mat[,-col]),nrow=(B))
          
          })

  
  
  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}


#' internal RandomForest prediction function
#'
#' private function that uses randomForest to predict fiven column from the rest
#' and assigns significance to
#' correlation coefficents based on resampling
#'
#'
#' @param mat matrix of response and predictor variables
#' @param col column number of response variables
#' @param B number of random samples (shuffling) of response variables
#'
#' @keywords internal
#'
#' @importFrom ranger ranger
rfResample<-function(mat,scaleData=scaleData,col=1,B=B){
  #require(randomForest)
  #require(ranger)

  # decide if the matrix needs to be dropped and NA returned due
  # to low or zero variation in gene expression
  if (ncol(mat)<=2){
    return(matrix(NA,ncol=ncol(mat)-1,nrow=3,
                  dimnames=list(c("coefs","pval","p2"),1:(ncol(mat)-1) ))
    )
  }

  if(scaleData){
    mat=scale(mat)
    mat[is.nan(mat)]=0 # when scaled all 0 columns will be NaN, conv. to 0
  }

  # resample response variables Ys
  Ys=lapply(1:B,function(x) sample(mat[,col],nrow(mat)))
  # original coefs in this case importance values as %incMSE
  #mod<- randomForest(x = mat[,-col], y = mat[,col],
  #                   importance=TRUE,na.action=na.omit)
  #  orig=importance(mod)[,2] # gini index
  
  colnames(mat) <- c("gene",rep("enh",ncol(mat)-1))
  
  mod <- ranger( "gene ~ .", data = data.frame(mat), importance="impurity",
                num.trees = 500,write.forest = FALSE,
                num.threads=1)
  orig=mod$variable.importance # gini index
  #coefs from resampling
  coefs=matrix(0.0,ncol=ncol(mat[,-col]),nrow=(B))


  # calculate coeff/importance for resampled Ys
  # handling case where by intense sampling occasionally only zeros are selected
  tryCatch({
  for(i in 1:B){
   
    mod <-ranger( y ~ ., data = data.frame(y=Ys[[i]],mat[,-col]),
                  importance="impurity",
                  num.trees = 500,write.forest = FALSE)
    coefs[i,]=mod$variable.importance # gini index

  }
  },error=function(coefs){coefs <- matrix(NA,ncol=ncol(mat[,-col]),nrow=(B))})

  # calculate p-vals
  pvals=estimateGammaPval(coefs,orig=orig,abs=TRUE,add=0)

  # arrange output format and return
  return(pvals)
}



  


  ########### END OF association prediction functions
#########-------------------------------------------------###########
