#' Function that rearranges list of scores calculated per cell type
#' 
#' @param scores.per.cell.type list of result score matrices for each biwgig
#'  file
#' @param sampleIDs  A vector of unique sample ids/names of bigwig files 
#' @return a dataframe of quatified scores per genome window and bigwig files.
#' Each row corresponds to ranges in the genome, and columns correspond to the 
#' bigwig files. 
#' @keywords internal
#' @author Inga Patarcic
scoresAsMcols <- function(scores.per.cell.type,sampleIDs){
  
        mcols.per.cell.type <- rbind.data.frame(scores.per.cell.type,
                                                stringsAsFactors=F,
                                                make.row.names=T)
        # adjust column names - remove bw extension
        colnames(mcols.per.cell.type) <- sampleIDs
       
     return(mcols.per.cell.type)
}



#' Function that runs normalization for the quantified activity/expression level
#' 
#' @param scoresDF a dataframe of quatified scores per genome window and bigwig 
#'  file
#' @param normalizationMet  charater. which normalization procedure is wanted.
#' @return a normalized dataframe of quatified scores per genome window and 
#' bigwig files. 
#' @keywords internal
#' @author Inga Patarcic
normalizeScores <- function(scoresDF,normalizationMet){

    if (normalizationMet=="ratio") {   
            # estimate size factors and multiply gene expression using these 
            # factors
    sizeFactors <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(scoresDF))
    scoresDF.normalized <- as.matrix(scoresDF)*rep(sizeFactors,each=
                                                            nrow(scoresDF))}
                
    if (normalizationMet=="quantile") {  
 scoresDF.normalized <- preprocessCore::normalize.quantiles(as.matrix(scoresDF))
          rownames(scoresDF.normalized) <- rownames(scoresDF)
          colnames(scoresDF.normalized) <- colnames(scoresDF)}
           
 
  return(scoresDF.normalized)

  
  }
  



#' Function that quantifies exon expression based on strand specificity
#' 
#'  Function that separates stranded and unstranded libraries and quantifies 
#'  exon expression separately for stranded libraries eg expression of exons 
#'  located on the + strand is quantified using forward libraries if libraries
#'  are unstranded then exon orientations is not taken into account
#'  
#' @param exons GRanges object with exon coordinates and metadata about genes
#' @param exons.splitted GRanges object with exon coordinates splitted based on
#' the strand.
#' @param geneActSignals character, paths to bigwig files
#' @param libStrand character, library strandness ("+","-","*")
#' @param sampleIDs character, name of samples/bigwig files
#' @param summaryOperation argument for ScoreMatrixBin() how to quantify scores
#' @param mc.cores integer, number of cores for mclapply
#' @return Granges object with metadata which corresponds to quantified scores 
#' per input exon regions and bigwig files
#' @keywords internal
#' @author Inga Patarcic
exonExpressionStrandAdjusted <- function(exonsForQuant,
                                          exonsSplitted,
                                          geneActSignals,
                                          sampleIDs,
                                          libStrand,
                                          mc.cores,
                                          summaryOperation) {
  
  
  if (length(exonsSplitted$`*`)!=0){stop("Strand of exons can be only + or -")}
  
        # detect stranded tracks and unstranded tracks    
        ForwardLibraries <- stringr::str_detect(libStrand,"\\+") 
        Unstranded.libraries <- stringr::str_detect(libStrand,"\\*")
        
  if (sum(ForwardLibraries)!=(sum(!(ForwardLibraries|Unstranded.libraries)))){
    stop("Number of forward and reverse libraries do not match")}
  
  
  # separate analysis for unstranded & stranded libraries
  
  if (sum(Unstranded.libraries)!=length(libStrand)){
    # test strandness - "mixed"&"stranded"
    # Stranded libraries - positive/forward strand
    
    # forward strand
    exonscoresForwardStrand <- mclapply(geneActSignals[ForwardLibraries], 
                                        function(x) {try(
                                        genomation::ScoreMatrixBin(x,
                                              windows =  exonsSplitted$`+`,
                                              bin.num = 1,
                                              type = "bigWig",
                                              bin.op=summaryOperation, 
                                              is.noCovNA=F),
                                              silent = T)},
                                        mc.cores=mc.cores)
    
    if (max(sapply(exonscoresForwardStrand,class)=="try-error")){
      stop(paste("Error with ScoreMatrixBin and following bigwig files",
                 geneActSignals[ForwardLibraries][sapply(
                   exonscoresForwardStrand,class)=="try-error"]))}                                      
    
    exonscoresForwardStrand.df <- scoresAsMcols(scores.per.cell.type=
                                                  exonscoresForwardStrand,
                                        sampleIDs=sampleIDs[ForwardLibraries])
    
    # reverse strand               
    exonscoresReverseStrand <- mclapply(geneActSignals[!(ForwardLibraries|
                                                        Unstranded.libraries)], 
                                        function(x) {try(
                                          genomation::ScoreMatrixBin(x,
                                               windows =  exonsSplitted$`-`,
                                               bin.num = 1,
                                               type = "bigWig", 
                                               bin.op=summaryOperation,
                                               is.noCovNA=F),
                                               silent = T)},
                                        mc.cores=mc.cores)
    
    if (max(sapply(exonscoresReverseStrand,class)=="try-error")){
      stop(paste("This bigwig files do not overlap exon regions",
                 geneActSignals[!(ForwardLibraries|Unstranded.libraries)]
                 [sapply(exonscoresReverseStrand,class)=="try-error"]))}
    
    exonscoresReverseStrand.df <- abs(scoresAsMcols(scores.per.cell.type=
                                                      exonscoresReverseStrand,
                                        sampleIDs=sampleIDs[ForwardLibraries]))
    
    
    # test if exon length is equal to the length of splitted exons
    if (length(exonsForQuant)!=sum(c(nrow(exonscoresForwardStrand.df),
                             nrow(exonscoresReverseStrand.df)))) {stop("Exons
                                    and splitted exon object are not equal")}
    
    
    # pool reverse and forward strand together        
    values(exonsForQuant) <- (cbind(mcols(exonsForQuant),
                                    rbind(exonscoresForwardStrand.df,
                                    exonscoresReverseStrand.df)))
    
    }
  
  if (sum(Unstranded.libraries)!=0){
    
    
    # Unstranded
    exonscoresUnstranded <- mclapply(geneActSignals[Unstranded.libraries], 
                                     function(x){try(
                                       genomation::ScoreMatrixBin(x,
                                                      windows = exonsForQuant,
                                                      bin.num = 1,
                                                      type = "bigWig",
                                                      bin.op=summaryOperation,
                                                      is.noCovNA=F),
                                                    silent = T)},
                                     mc.cores=mc.cores)
    
    if (max(sapply(exonscoresUnstranded,class)=="try-error")){
      stop(paste("This bigwig files do not overlap exon regions",
                 geneActSignals[Unstranded.libraries]
                 [sapply(exonscoresUnstranded,class)=="try-error"]))}                                      
    
    
    exonscoresUnstranded.df <- scoresAsMcols(scores.per.cell.type=
                                               exonscoresUnstranded,
                                    sampleIDs=sampleIDs[Unstranded.libraries])                  
    
    values(exonsForQuant) <- (cbind(mcols(exonsForQuant),
                                    rbind(exonscoresUnstranded.df)))
    
    
  }
  
  return(exonsForQuant)                
}




#' Function that quantifies gene expression based on prequantified exon exp
#' 
#' Function that calculates per gene expression based on exon expression results
#' sum across exons (mean value across exon*exon width)/
#' gene_length(sum_of_all_exon_lengths)
#'  
#' @param expressionPerGene GRanges object with exon coordinates and 
#' metadata info about corresponding genes and quantifed exon expression
#' @return Granges object of genes with metadata which corresponds to quantified
#' gene expression
#' @keywords internal
#' @author Inga Patarcic
quantifyGeneExpression <- function(expressionPerGene,exonsForGene){
  
  ExonLength <- width(expressionPerGene)
  GeneLength <- sum(ExonLength,na.rm=T)
  ExpressionPerSample=as.data.frame(mcols(expressionPerGene)[-(1:exonsForGene)])
        
    GeneExpression <- apply(ExpressionPerSample*ExonLength,
                            2,sum,na.rm=T)/GeneLength
  
      return(GeneExpression)
}




#' Calculates regulatory activity over pre-defined regions
#' 
#' The function calculates regulatory activity from histone
#' modification, DNAse or methylation signals for pre-defined regulatory
#' regions and returns a GRanges object with regulatory region locations
#' and their activity over a set of samples.
#' 
#' 
#' @param regRegions a GRanges object that contains regulatory regions
#' over which the regulatory activity will be calculated. 
#' It is strongly suggested to adjust 
#' seqlengths of this object to be equal to the seqlenghts Hsapiens from the 
#' appropriate package (BSgenome.Hsapiens.UCSC.hg19 or whatever needed version).
#' @param activitySignals a named list of BigWig files. Names correspond to 
#'        unique sample ids/names.
#' @param sampleIDs NULL (default). A vector of unique sample
#' ids/names(.bw files), ordered as the bigwig files are ordered. When NULL
#' basenames of .bw files is used as a unique sample ids/names.
#' @param isCovNA (def:FALSE), if this is set to TRUE, uncovered
#' bases are set to NA, this is important when dealing with methylation
#' data, where uncovered bases on a bigWig file do not mean zero methylation.
#' @param weightCol a numeric column in meta data part used as weights. Useful
#' when genomic regions have scores other than their coverage values, 
#' such as percent methylation, conservation scores, GC content, etc.
#' @param normalize NULL(default). Optional "quantile" and "ratio"
#' If set to "quantile" activity measures are quantile normalized as 
#' implemented in \code{\link[preprocessCore]{normalize.quantiles}} and
#' returned ; if set to "ratio" then "median ratio method" implemented
#' as \code{\link[DESeq2]{estimateSizeFactorsForMatrix}} is used to normalize 
#' activity measures.
#' @param summaryOperation "mean"(default).
#' This designates which summary operation should be used over the regions.
#' Currently, only mean is available, but "median" or "sum" will be implemented
#' in the future.
#' @param mc.cores (def:1) Define the number of cores to use;
#' at most how many child processes will be run simultaneously 
#' using mclapply from parallel package. Parallelization requires at 
#' least two cores.
#' 
#' @return a GRanges object where its meta-columns correspond
#'         to calculated acvitity measures and column names 
#'         correspond to provided sample ids or names.
#' 
#' 
#' @details regulatory activity is measured by averaging logFC for
#' histone modification ChIP-seq profiles, or DNAse signal, or methylation
#' per base.Currently, relevant bigWig files are required to calculate activity       
#' activity. This function might be extended to work with BAM files
#' in the future. 
#' 
#' #@examples
#' test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
#' test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
#' regTSS_toy <- GRanges(c(rep("chr1",4),rep("chr2",2)),
#'                       IRanges(c(1,7,9,15,1,15),c(4,8,14,20,4,20)),
#'                                             c(rep("+",3),rep("-",3)))
#' regTSS_toy$reg <-  regTSS_toy[c(1,1,3:6)]
#' regTSS_toy$name2 <- regTSS_toy$name <- paste0("TEST_Reg",
#'                                         c(1,1,3:length(regTSS_toy)))
#' regActivity(regTSS_toy,c(test.bw,test2.bw))   
#' 
#' @export
regActivity <- function(regRegions,
                        activitySignals,
                        sampleIDs=NULL,
                        isCovNA=FALSE,
                        weightCol=NULL,
                        summaryOperation="mean",
                        normalize=NULL,
                        mc.cores=1){
  
  
       # test input - ranges
  if (!exists("regRegions")) {stop("regRegions object missing")}
      # test input - bigwig files
  if (!min(sapply(activitySignals,file.exists))) 
          {stop("activitySignals object missing")}
  if (is.null(sampleIDs)) {
      bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
      sampleIDs <- stringr::str_replace(basename(activitySignals),
                               paste(bw.exts,collapse="|"),"")}
          
          
  # calculating coverage
         
          scores.per.cell.type <- parallel::mclapply(activitySignals,
                                    genomation::ScoreMatrixBin,
                                           windows = regRegions,
                                           bin.num = 1,
                                           type = "bigWig", 
                                           is.noCovNA=isCovNA,
                                           weight.col=weightCol,
                                           bin.op=summaryOperation,
                                           mc.cores = mc.cores)
          
         
           # adding scores as mcols  
        scores.per.cell.type.df <- scoresAsMcols(scores.per.cell.type=
                            scores.per.cell.type,sampleIDs=sampleIDs)
                                                   
        
          
          # add normalization
      if (!is.null(normalize)){
        scores.per.cell.type.df <- normalizeScores(scores.per.cell.type.df,
                                                   normalize)}

          mcols(regRegions) <- scores.per.cell.type.df
           
          return(regRegions)
          
  
  
}

#' Identify regulatory regions around provided TSSes
#' 
#' The function identifies the regulatory regions around provided
#' TSSes over a pre-defined window. The function needs a GRanges
#' object for TSSes with meta-columns corresponding to expression
#' levels in different cell types or conditions. 
#' 
#' @param regActivity a GRanges object output from \code{regActivity}
#'        function
#' @param tss a GRanges object output from \code{bwToGeneExp} that
#' contains the TSS location and associated gene expression values 
#' per cell type or condition as meta data. Each row should have a
#' "name" and "name2" columns for unique id or name/symbol for the
#' gene which the TSS is associated with. One is Ensembl id 
#' and the other could be used for the gene symbol. Other metadata
#' column names should represent sample names/ids and should
#' match the GRanges object provided via regActivity argument.
#' @param upstream number of basepairs upstream from TSS to look for
#' regulatory regions. default 500kb
#' @param downstream number of basepairs downstream from TSS to look for 
#' regulatory regions. default 500kb
#'@param mc.cores (def:1) Define the number of cores to use;
#' at most how many child processes will be run simultaneously 
#' using mclapply from parallel package. Parallelization requires at 
#' least two cores.
#' 
#' @return a GRangesList object per gene that contain location of TSS
#' and regulatory regions around that gene. Names for the GRangesList
#' are unique gene ids/names. 
#' Metadata for a GRanges object in the list represents regulatory 
#' activity and gene expression accross the same samples. 
#' The GRanges objects have the following metadata columns:
#'  1. featureType: either "gene" or "regulatory"
#'  2. name: name/id for gene and enhancers. Gene name could be id from a 
#'  database enhancer name should be in the format as follows "chr:start-end"
#'  3. name2: a secondary name for the feature, such as gene symbol "PAX6" etc. 
#'  not necessary for enhancers could be NA
#'  4. other columns: numeric values for gene expression or regulatory actvity.
#'    Column names represent sample names/ids.
#' 
#' @examples 
#' 
#' regTSS_toy <- GRReg1_toy
#'   regTSS_toy$bw1 <- rep(1,length(GRReg1_toy))
#'   regTSS_toy$bw2 <- rep(2,length(GRReg1_toy))
#'   regTSS_toy$bw3 <- rep(3,length(GRReg1_toy))
#' regReg_toy <- GRReg2_toy
#'    regReg_toy$bw1 <- rep(3,length(regReg_toy))
#'    regReg_toy$bw2 <- rep(4,length(regReg_toy))
#' 
#' regActivityAroundTSS(regReg_toy,regTSS_toy,upstream=1,downstream=1)
#' 
#' 
#' 
#' @details only enhancers located within (+/-)upstream/downstream of TSS 
#' are identified,extracted and reported in output (together with info
#' about gene expression). Sample id's (corresponding to the cell types 
#' or conditions) are included only if both, 1) gene expression values and
#' 2) quantified regulatory activity are available in TSS and 
#' regActivity objects. Non-overlapping cell types are excluded.
#'  
#' 
#' 
#' @export
regActivityAroundTSS <- function(regActivity,
                                 tss,
                                 upstream=500000,
                                 downstream=500000,
                                 mc.cores=1){
  
  
  if (is.null(tss$name)) {stop("TSS object does not contain info about gene 
                               name. TSS GRanges object should be 2nd arg")}
  
  # extend TSS for wanted region 
  tss.prom <- promoters(tss,
                        upstream,
                        downstream)
  # split per gene 
  tss.extended <- split(tss.prom,
                        as.character(tss.prom$name))
  
  # getting EnhRegions which overlap extended TSS
  
  tssActivity <- parallel::mclapply(tss.extended,function(x){
    
    # overlap extended TSS & regRegion
    tss.regAct.overlap <- as.data.frame(findOverlaps(x,regActivity))
    
    # if there is any overlap - GO!
    if (nrow(tss.regAct.overlap)!=0) {
      regActivity <- regActivity[tss.regAct.overlap$subjectHits]
      
      # adjusting mcols for reg region
            name <- as.character(regActivity)
            featureType <- rep("regulatory",length(name))
            name2 <- rep(NA,length(name))
            
            mcols(regActivity) <- cbind(featureType,
                                        name,
                                        name2,
                                        data.frame(mcols(regActivity)))
            
      # adjusting mcols for gene expression object
            featureType <- "gene"
            mcols(x) <- cbind(featureType,
                              data.frame(mcols(x)))
            
            tss.colnames <- colnames(mcols(x))[colnames(mcols(x))%in%
                                               colnames(mcols(regActivity))]
            mcols(x) <-  mcols(x)[tss.colnames]
            mcols(regActivity) <-  mcols(regActivity)[tss.colnames]
            
            geneReg <- c(x,regActivity)
            geneReg$featureType  <- as.character(geneReg$featureType)
            geneReg$name   <- as.character(geneReg$name)
            geneReg$name2   <- as.character(geneReg$name2)
            
      return(geneReg)
      
    }
  },
  mc.cores = mc.cores)
  
  # removing NULL entries
  tssActivity <- tssActivity[!sapply(tssActivity,is.null)]
  # returning GRangesList
  tssActivity <- GRangesList(tssActivity)
  
  return(tssActivity)
  
}





#' Quantifies gene expression measured by RNA-Seq  
#' 
#' 
#' The function quantifies gene expression as a sum of exon 
#' expressions quantified over pre-defined exon regions using 
#' signal from RNA-Seq tracks (bigwig files) and returns a
#' GRanges object with TSS location and corresponding gene
#' expression levels quantified over a set of samples. 
#'  
#' 
#' 
#' @param exons A GRanges object that contains exon regions over which 
#' the expression will be calculated. It is strongly suggested to adjust 
#' seqlengths of this object to be equal to the seqlenghts Hsapiens from the 
#' appropriate package (BSgenome.Hsapiens.UCSC.hg19 or whatever needed version).
#' A meta-columns with following info are necessary: 1) name (character)- ENSEMBL 
#' or other gene identifier; 2) name2(character,optional) - 2nd gene identifier,
#' 3) reg - GRanges object - gene location, or TSS location. 
#' @param geneActSignals a named list of RNA-Seq BigWig files. Names correspond 
#' to the unique sample ids/names. Stranded and unstranded libraries allowed.
#'  BUT!!! It is crucial that forward and reverse RNA-Seq libraries are listed
#' in a row (eg one on top of each other)
#' @param sampleIDs NULL (default). A vector of unique sample
#' ids/names(.bw files), ordered as the bigwig files are ordered. When NULL
#' basenames of .bw files is used as a unique sample ids/names.
#' @param libStrand a vector of "*","+,"-" (default NULL) which needs to be 
#' entered as an argument.This vector provides info about the order of RNA-Seq 
#' libraries based on their strandness: "+" corresponds to forward/positive 
#' RNA-Seq bigwig files; "-" corresponds to reverse/negative RNA-Seq bigwig
#' files and "*" is unstranded library. When all libraries are unstranded then
#' the vector should contain a list of "*" with the lenght equal to the number
#' of analyzed RNA-Seq libraries (eg bigwig files). 
#' If libStrand=NULL than function will do it automatically, eg
#' create a vector of "*".It is crucial that stranded RNA-Seq libraries are
#' listed in a row (eg one on top of each other)
#' @param summaryOperation "mean"(default).
#' This designates which summary operation should be used over the regions.
#' Currently, only mean is available, but "median" or "sum" will be implemented
#' in the future.
#' @param mc.cores (def:1) Define the number of cores to use;
#' at most how many child processes will be run simultaneously 
#' using mclapply from parallel package. Parallelization requires at 
#' least two cores.
#' @param normalize NULL(default). Optional "quantile" and "ratio"
#' If set to "quantile" activity measures are quantile normalized as 
#' implemented in \code{\link[preprocessCore]{normalize.quantiles}} and
#' returned ; if set to "ratio" then "median ratio method" implemented
#' as \code{\link[DESeq2]{estimateSizeFactorsForMatrix}} is used to normalize 
#' activity measures.
#' 
#'         
#' @return a GRanges object where its meta-columns correspond to quantified
#' gene expression across cell type or condition. GRanges correspond
#' to the TSS location of the analyzed gene. Additionally, each TSS 
#' contains following metadata:a "name" and "name2" columns for 
#' unique id or name/symbol for the gene which the TSS is associated with. 
#' One is Ensembl id and the other could be used for the gene symbol. 
#' Other metadata column names should represent sample names/ids and should
#' match the GRanges object provided via regActivity argument.
#' 
#' @details For a gene of interest and across all samples (.bw files)
#' individaully,the expression is firstly calculated for all exon 
#' regions which correspond to the gene of interest. Then, per gene
#' exon expression scores are summed together and divided by a total exon length
#' as followed:
# #'  \deqn{Gene Expression=\sum{n=1}^{Exons{n}}
# #'  (MeanExonExpression{n}*ExonLength{n})/(GeneLength)}
#' Normalization is runned on the level of gene expression. Currently, the 
#' relevant bigWig files are required to calculate activity. This function might 
#' be extended to work with BAM files in the future. RNA-Seq .bw files can 
#' originate from stranded,unstranded or mixed libraries.
#' 
#' 
#' 
#' @examples test.bw <- system.file("extdata", "test.bw",package = "reg2gene")
#' test2.bw <- system.file("extdata", "test2.bw",package = "reg2gene")
#' regTSS_toy <- GRanges(c(rep("chr1",4),rep("chr2",2)),
#'                       IRanges(c(1,7,9,15,1,15),c(4,8,14,20,4,20)),
#'                                             c(rep("+",3),rep("-",3)))
#' regTSS_toy$reg <-  regTSS_toy[c(1,1,3:6)]
#' regTSS_toy$name2 <- regTSS_toy$name <- paste0("TEST_Reg",
#'                                         c(1,1,3:length(regTSS_toy)))
#'                                         
#'                                         
#' bwToGeneExp(exons = regTSS_toy,geneActSignals = c(test.bw,test2.bw))
#' @export
bwToGeneExp <- function(exons,
                        geneActSignals,
                        sampleIDs=NULL,
                        libStrand=NULL,
                        summaryOperation="mean",
                        mc.cores=1,
                        normalize=NULL){
  
  
  
  # test if there is any info about genes in exon granges object
  if (sum(stringr::str_detect(colnames(mcols(exons)),"name"))==0){
    stop("No info about the gene name")}
  
  
  if (is.null(sampleIDs)) {
    bw.exts = c(".bw",".bigWig",".bigwig",".BigWig", ".BIGWIG", ".BW")
    sampleIDs <- stringr::str_replace(basename(geneActSignals),
                                      paste(bw.exts,collapse="|"),"")
  }
  
  
  
  # separate genes form forward and reverse strand
  exons.splitted <- split(exons,strand(exons))
  gene.metadata <- length(mcols(exons)) # metadata info about corresponding gene
  
  # arranging a problem of strandness of RNA-Seq libraryies
  # 1. if object libStrand is empty then strands are assumed to be "*"
  if (is.null(libStrand)) {libStrand <- rep("*",length(geneActSignals))}
  
  # quantify expression per exon with adjustment to strandness of RNA-Seq 
  # libraries
  ExonExpression <- exonExpressionStrandAdjusted(exons,
                                              exons.splitted,
                                              geneActSignals,
                                              sampleIDs, 
                                              libStrand,
                                              mc.cores,
                                              summaryOperation)         
  # quantify expression per gene   
  exonExpressionPerGene <- split(ExonExpression,
                                 as.character(ExonExpression$name))
  
  GeneExpression <- list()
     for (i in 1:length(exonExpressionPerGene)){
           gene.id <- names(exonExpressionPerGene)[[i]]
          GeneExpression[[gene.id]] <- quantifyGeneExpression(
                                              exonExpressionPerGene[[gene.id]],
                                              gene.metadata)}
  
  
    GeneExpression.df <- do.call("rbind",GeneExpression)

   #############            
   # normalization
   if (!is.null(normalize)){
     GeneExpression.df <-  normalizeScores(GeneExpression.df,
                                           normalize)}


        # arranging gene coordinates as TSS
    Genes.object <- data.frame(unique(mcols(ExonExpression)[1:gene.metadata]),
                               stringsAsFactors = F)
    
    
        #extract gene location
              # TSS coord of genes
    
                 tss <- promoters(GRanges(Genes.object$reg.seqnames,
                                IRanges(Genes.object$reg.start,
                                        Genes.object$reg.end),
                                strand=Genes.object$reg.strand),
                                1,1)


        # adding meta-data
             name <- Genes.object$name
             name2 <- Genes.object$name2
                       if (length(name2)==0){
                         name2 <-  rep(NA,length(name))}# if no 2nd name
             metadata <- data.frame(name,
                                    name2,
                  GeneExpression.df[match(name,rownames(GeneExpression.df)),],
                            stringsAsFactors = F)

          mcols(tss) <- metadata

  return(tss)
  
} 







