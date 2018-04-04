#' Plots enhancer-promoter interactions
#' 
#' 
#' @param  interactionData A GInteractions object or a list of GInteractions
#' objects. Each GInteractions object has an anchor1 which corresponds to 
#' regulatory region (enahncer), and anchor2 which should correspond to the TSS 
#' of a gene (that info is used to plot promoters, thus please input TSS, not 
#' gene coordinated). A minimum necessary info that should be present as 
#' meta-data is info about gene names stored in the column "name2". Additionaly,
#' statistics (as "pval" column), and colors of the loops (as "color" column) 
#' could be added as meta-data. Read more about arguments statistics and 
#' coloring.  
#' Such object can be produced by modelling [\code{\link{associateReg2Gene}}], 
#' or meta-analysis [\code{\link{metaAssociations}}] or voting functions from 
#' reg2gene package. 
#' 
#' @param  rangeGenes (default NULL) By default, function retrieves 
#' automatically gene annotations using TxDb.Hsapiens.UCSC.hg19.knownGene 
#' package for genes present in the interactionData object. Importantly,
#' hgnc_symbol or ensembl_gene_id are stored as the name2 column in the
#' interactionData object.
#' However, by using this argument, one can define and import gene annotations 
#' in formats such as TxDb, GRanges, GRangesList, data.frame, 
#' character.In the case you want to import your own gene annotations, please
#' take a look at \code{\link[Gviz]{GeneRegionTrack}} for more details 
#' how this object should look like.
#' 
#' @param selectGene (default NULL) a character vector of gene name symbols 
#' (eg "FTO","IRX3", etc.) which user wants to plot.
#' 
#' @param selectRegulatoryRegion (default NULL) GRanges object or a 
#' character vector which stores info about the region of interest which 
#' user wants to plot (format: "chr16:53112601-53114200"). This region does not
#' necessarilly need to be equal to the regulatory regions reported in the 
#' interactionData input objects, whereas it only needs to overlap some
#' regulatory regions. In addition, to provide the context, a function plots 
#' not only regulatory regions queried, but all that other that are associated
#' with genes which are in the first place associated with queried regions. 
#' Primary interactions (those who overlap with queried region) are plotted red,
#' whereas secondary interactions (additional interactions associated with 
#' identified genes) are plotted in grey. 
#' 
#' @param  filters default "hgnc_symbol"; other option "ensembl_gene_id". 
#' One should define if the gene id's in the name2 column of the 
#' interactionData objects are ensembl gene ids or hgnc symbols.
#' In summary, filters define a restriction on the query for
#' \code{\link[biomaRt]{getBM}}. This step is necessary when gene annotations 
#' are retrieved automatically (default).
#' 
#' @param benchmarkData (default NULL) A GInteractions object or a list of
#' GInteractions objects storing info about interaction regions in the genome 
#' (for example, regions from the literature obtained using chromatin 
#' conformation capture related methods.
#' 
#' @param statistics (default NULL) Column name where info about the heights of 
#' loops is stored. If NULL, then all interaction loops have an equal height 
#' Otherwise, loops of different height are plotted based on the info stored in
#' the corresponding column, eg "pval" or "qval". Height is calculated as 
#' -log10(pval).
#' 
#' @param coloring (default NULL) Information about colors used to plot enhancer
#' promoter interactions. If it is NULL, then all interactions will be plotted
#' red. Otherwise, a character, e.g. a name of the column of the interactionData 
#' object. In selected column, info about loop color is stored. 
#' For example, interactions for different genes can be 
#' colored differently, or statistically significant/insignificant interactions
#' can be colored differently. User pre-defines this column by itself.
#' In addition, if selectRegulatoryRegion!=NULL, then primary selected 
#' interactions are colored red, wheres secondary interactions (additional
#' interactions associated with identified genes) are plotted in grey
#' 
#' @param sizes (default 0.3) Default size of tracks. In the case that many 
#' tracks are plotted this argument should be adjusted (lower height of tracks)
#' such that all tracks get plotted.
#' 
#' @param cex.title (default NULL) Numeric. This argument controls the size of 
#' text which describes each track. If, many tracks are plotted, then text might
#' disapear, but then cex.title argument should be set to smaller value, until
#' text appears again 
#' 
#' @param plotEnhancersOne (default TRUE). To plot a separate track(s)
#' for enhancer regions or not
#'
#' @import GenomicFeatures
#' @import Gviz
#' @import GenomicInteractions 
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import biomaRt
#'    
#' @return Nice plot.
#' 
#' @details  Plots InteractionTrack, GeneRegionTrack, InteractionsTracks,
#' AnnotationTrack, IdeogramTrack and GenomeAxisTrack for queried (list of) 
#' GInteraction object(s) [which contain info about enhancer-gene associations],
#' and (list of) benchmark datasets.
#' 
#' User can select genes or regions of interest, and these regions will be 
#' selected from the input GInteractions object. Check description of arguments
#' for this function for specific cases (coloring sheme).
#' Names of the tracks correspond to the name of the list (in the case that list
#' is used as an input), or one can define name for individual objects.
#' E track (stands for enhancers) indicates locations of enhancer regions, 
#' wereas P track (stands for promoters) plots TSS locations for input genes. 
#'  
#' 
#' @author Inga Patarcic
#' 
#' @details Take a look at the plot :-)
#' Gene annotations are retieved from TxDb.Hsapiens.UCSC.hg19.knownGene, and
#' exons are plotted (no info about transcript names is used). 
#' Entrez -> hgcn symbol are obtained using biomaRt. 
#' 
#' @examples
#' library(Gviz)
#' library(GenomicInteractions)
#' library(GenomicFeatures)
#' library(GenomicRanges)
#' library(InteractionSet)
#' library(reg2gene)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(biomaRt)
#' 
#' # Creating an example GInteractions set:
#' 
#'  enhancers <- GRanges(rep("chr16",6),
#'                      IRanges(c(53112601,55531601,53777201,
#'                                53778801,54084001,53946467),
#'                              c(53114200,55533250, 53778800,
#'                                53780400, 54084400 ,53947933)))
#' 
#'  genes <- GRanges(rep("chr16",6),
#'                   IRanges(c(53737874, 54964773, 54320676,
#'                             53737874, 54964773, 54320676),
#'                           c(53737874, 54964773, 54320676,
#'                             53737874, 54964773, 54320676)))
#' 
#'  GenomeInteractions <- GInteractions(enhancers,genes)
#' 
#'  GenomeInteractions$name2 <- c("FTO","IRX5","IRX3")
#' 
#'  GenomeInteractions$pval <- c(0.20857403, 0.72856090, 0.03586015,
#'                               0.32663439, 0.32534945, 0.03994488)
#' 
#'  GenomeInteractions$color <- c("red","blue","grey")
#' 
#'  # simplest example
#'  # if no info about coloring and height of loops provided, all loops are equal
#'  plotGenomeInteractions(interactionData = GenomeInteractions)
#' 
#' 
#'  # provide info about statistics * height of loops
#'  plotGenomeInteractions(interactionData = GenomeInteractions,
#'                   statistics ="pval",
#'                   coloring = "color")
#' 
#' 
#' 
#'   # Specific gene can be individually plotted
#'     plotGenomeInteractions(interactionData = GenomeInteractions,
#'                            selectGene = "FTO")
#' 
#'   # This should be equal to the example where genes are not selected
#' 
#'   # if one wants to plot regulatory region
#'   plotGenomeInteractions(interactionData = GenomeInteractions,
#'                selectRegulatoryRegion = "chr16:53112601-53114200")
#'   # NOTE: red interactions correspond to the enhancer regions that
#'   # overlap queried region, whereas grey interaction is interactions represent
#'   # in the dataset and associated with a gene with which the queried regulatory
#'   # region interacts (provides info about neighbourhood for queried interaction)
#' 
#' 
#'   # and function plots not regulatory regions queried, but all that
#'   # are present in the interactionData and overlap that region
#' 
#'   plotGenomeInteractions(interactionData = GenomeInteractions,
#'                          selectRegulatoryRegion = "chr16:53112601-53778800",
#'                          coloring = "color",
#'                          statistics = "pval")
#' 
#'   # more than one interaction track can be plotted, if input is list
#'   GenomeInteractionsList <- list(GenomeInteractions,GenomeInteractions[1:3])
#'   names(GenomeInteractionsList) <- c("EnhP1","EnhP2")
#' 
#'   # names from the list are taken as an argument for Enh&Inter plots
#'   plotGenomeInteractions(interactionData = GenomeInteractionsList,
#'                          coloring = "color",
#'                          statistics = "pval",
#'                          interactionDataName = "EP Interactions",
#'                          cex.title = 0.4,
#'                          sizes = 0.4)
#' 
#'   # additionaly, benchmark data can be plotted as well
#' 
#'   plotGenomeInteractions(interactionData = GenomeInteractions,
#'                           coloring = "color",
#'                           statistics = "pval",
#'                           benchmarkData = GenomeInteractions[1:3])
#' 
#' 
#'   # add a list of bemchmark data
#'     benchmarkData = list(GenomeInteractions[1:3],GenomeInteractions[4:5],
#'     GenomeInteractions[6])
#'   names(benchmarkData) <- c("HiC","ChIA-PET","GTEx")
#' 
#' 
#'  plotGenomeInteractions(interactionData = GenomeInteractions,
#'                         coloring = "color",
#'                         statistics = "pval",
#'                         benchmarkData = benchmarkData,
#'                         # interactionDataName = "EP Interactions",
#'                         cex.title = 0.4,
#'                         sizes = 0.4)                  
#'                                                                       
#' @export
plotGenomeInteractions <- function(interactionData,
                                   rangeGenes=NULL,
                                   selectGene=NULL,
                                   selectRegulatoryRegion=NULL,
                                   filters="hgnc_symbol",
                                   benchmarkData=NULL,
                                   statistics=NULL,
                                   coloring=NULL,
                                   interactionDataName="Enh~Prom Interactions",
                                   benchmarkName="Benchmark",
                                   sizes=0.3,
                                   cex.title=1, 
                                   plotEnhancersOne=TRUE){
  
#######################################
# setting universal
  
            require(Gviz)
            require(GenomicInteractions)
            gen <- "hg19"
            
            benchParameters <- list(col.anchors.fill ="blue",
                                    col.anchors.line = "grey",
                                    interaction.dimension="height", 
                                    interaction.measure ="counts",
                                    plot.trans=FALSE,
                                    plot.outside = TRUE, 
                                    col.outside="grey", 
                                    anchor.height = 0.1)
           enhParameters <- list(fill = "blue", 
                                  col = "blue", 
                                  just.group="above",
                                  cex.feature = 5,
                                  labelPos="above")
            grParameters <- list(fill = "grey",
                                 col = "lightgrey", 
                                 cex.feature = 5,
                                 labelPos="above")
      
# --------------------- 
            
#################################
# setting track by track:
# 1. Enhancer tracks            
  
# selection of regions
    if (class(interactionData)=="GInteractions"){
      
      # select gene
          if (!is.null(selectGene)){
        
           interactionData <- selectGeneF(x=interactionData,
                                          selectGene=selectGene)
      
      }
      
      # extracting info about regulatory regions of interest
          if (!is.null(selectRegulatoryRegion)){
        
            interactionData <- selectRegR(x=interactionData,
                                  selectRegulatoryRegion=selectRegulatoryRegion)
            coloring <- "color"
      }
      
      # extracting info about chromosome & genome set to hg19
          chr <- as.character(unique(seqnames(first(interactionData))))
      
      # add info about loop height
          interactionData <- getStatistics(x=interactionData,statistics)
         
     
      # createIneractionTrack  + adjust plotting parameters   
          interaction_track <- createIntTrack(interactionData,
                               interactionDataName=interactionDataName,
                               chr=chr,
                               coloring=coloring)
          
          # setting enhancer track
          enhTrack <- AnnotationTrack(first(interactionData), 
                                      genome=gen, 
                                      name="E",
                                      stacking="full",
                                      cex.feature = 5 )
     
          displayPars(enhTrack) <- enhParameters
    }
  
    if (class(interactionData)=="list"){
     
     # select gene
      if (!is.null(selectGene)){
        
       interactionData <-  lapply(interactionData,
                           selectGeneF,
                          selectGene=selectGene)
      
       }
     
     # extracting info about regulatory regions of interest
     if (!is.null(selectRegulatoryRegion)){
       
        interactionData <- lapply(interactionData,
                                  selectRegR,
                                  selectRegulatoryRegion=selectRegulatoryRegion)
        coloring <- "color"
        }
     
     # extracting info about chromosome & genome set to hg19
     chr <- as.character(unique(seqnames(first(interactionData[[1]]))))
     
     
     # add info about loop height removing NA values 
     interactionData <- lapply(interactionData,
                               getStatistics,
                               statistics=statistics)
     
     # add info about track names
     if (is.null(names(interactionData))) {
       names(interactionData) <- 1:length(interactionData)
     }
     interactionData <- lapply(names(interactionData),function(x){
       
                   tmp <- interactionData[[x]]
                   tmp$listName <-  x
                   
                   return(tmp)
     })
     
     
            
    # createIneractionTrack  + adjust plotting parameters   
     interaction_track <- lapply(interactionData,
                                createIntTrack,
                                interactionDataName=NULL,
                                chr=chr,
                                coloring=coloring)
     

     # setting enhancer track
     enhTrack <-  lapply(interactionData, function(x){
      
            if ("listName"%in%names(mcols(x))){
             
              enhTrack <- AnnotationTrack(first(x), 
                                 genome=gen, 
                                 name=paste0("E",unique(x$listName)),
                                 stacking="full",
                                 cex.feature = 5 )
              
            }
           
            if (!("listName"%in%names(mcols(x)))){
             
             enhTrack <- AnnotationTrack(first(x), 
                                         genome=gen, 
                                         name=NULL,
                                         stacking="full",
                                         cex.feature = 5 )
             
           }
           
               displayPars(enhTrack) <- enhParameters
               
               return(enhTrack)
           
                   })
  
  
     # arranging whether you want 1 enh region,or all
     if (plotEnhancersOne==TRUE) { enhTrack <- enhTrack[[1]]}
     
    }
  
 

                     
#################################
# 2. Gene track tracks               

# set span
# get span for plots, regardeless of how many input enhancer regions entered
   
# adjust span
     
Span <- range(GRanges(as.vector(sapply(interactionData, function(x){
                               as.character(range(c(first(x),second(x))))}))),
              ignore.strand=TRUE)
   
            
     if (!is.null(rangeGenes)){ Span <- range(c(Span,
                                             rangeGenes),
                                           ignore.strand=TRUE)}
     if (is.null(rangeGenes)){ 
    
       # obtaining info about genes of interst across list
            AllGenes <- unique(unlist(sapply(interactionData,
                                             function(x)x$name2)))
    
       # getting info about genes from Ensembl
            rangeGenes <- gettingGeneInfoEnsembl(AllGenes,filters)
        
      # Additional,setting ranges for plot based on enhancer
      # & gene span: min&max range 
              rangeGenesTmp <- rangeGenes
              mcols(rangeGenesTmp) <- NULL
            
            Span <- range(c(Span, rangeGenesTmp),ignore.strand=T)
     
  }
     
          
    tss <- getTSS(AllGenes,interactionData)

            
# extend min&max for 10000
PlotRange <- resize(Span,fix = "center",width=width(Span)+10000) 
     
      
     grtrack <- GeneRegionTrack(range=rangeGenes,
                                genome = gen,
                                chromosome = chr,
                                name = "Genes",
                                transcriptAnnotation="symbol",
                                stacking="full",
                                labelPos="above",
                                start = start(PlotRange), 
                                end = end(PlotRange)) 
     
     
#################################
# 3. Setting chromosome & location track
      
       gtrack <- GenomeAxisTrack()
          displayPars(grtrack) <- grParameters
      
       itrack <- IdeogramTrack(genome = gen,chromosome = chr)
       
#################################
# 4. setting promoter track
       
       
        promoterTrack <- AnnotationTrack(tss, 
                                         genome=gen, 
                                         name="P",
                                         #id=unique(interactionData$name2),
                                         stacking="squish")
       
       
       displayPars(promoterTrack) <- grParameters
      
#################################
# 5. Setting benchmark track(s)
  
        if (is.null(benchmarkData)){ 
    
                tracksToPlot <- unlist(list(itrack,
                                     gtrack, 
                                     grtrack,
                                     interaction_track,
                                     enhTrack,
                                     promoterTrack))
    
                trackSizes <- c(0.1,0.2,0.2,
                                rep(sizes,length(interaction_track)),
                                0.075,0.075)
                                  }
  
  # extra plot for benchmark
  
       if (!is.null(benchmarkData)){
    
            if (class(benchmarkData)=="GInteractions"){
              
              Benchmark_track <- benchPlotHelp(name=NULL,
                                               benchmarkData=benchmarkData,
                                               PlotRange=PlotRange, 
                                               chr=chr)
              
           
              displayPars(Benchmark_track) = benchParameters
              
              
              }
            
            if (class(benchmarkData)=="list"){
              
          
              Benchmark_track <- lapply(names(benchmarkData),
                                        function(x){
                                          
                                          
                      Benchmark_track <-  benchPlotHelp(name=x,
                                                benchmarkData=benchmarkData,
                                                PlotRange=PlotRange,
                                                chr=chr)
                                        
                          displayPars(Benchmark_track) = benchParameters
                       
                          return(Benchmark_track)
                                          
                                          })
              }
      
    

    
    tracksToPlot <- unlist(c(list(itrack,
                        gtrack, 
                        grtrack,
                        interaction_track,
                        enhTrack,
                        promoterTrack),
                        Benchmark_track))

    trackSizes <- c(0.1,0.2,0.2,
                    rep(sizes,length(interaction_track)),
                    0.075,0.075,
                    rep(sizes,length(Benchmark_track)))

      
                 
                 
                
    }

  #---------------------------------
  # plotting
    
  plotTracks(tracksToPlot,
             chromosome=chr, 
             from=min(start(PlotRange)), 
             to=max(end(PlotRange)), 
             sizes=trackSizes,
             background.panel = "#FFFEDB", 
             background.title = "darkblue",
             fontcolor.feature = "darkblue",
             cex.title=cex.title)
  
  }





#' Plots enhancer-promoter interactions
#' 
#' Help function to create benckmark tracks, especilly needed for a list of 
#' benchmark inputs.
#' @author IngaPa
#' @keywords internal
benchPlotHelp <- function(name,
                          benchmarkData,
                          PlotRange,
                          chr){
  
  if (!is.null(name)) {benchmarkDataOne <- benchmarkData[[name]]}
  if (is.null(name)) {benchmarkDataOne <- benchmarkData}
  # filter for regions of interest based on ranges of plot
  
  BenchFilter <- DataFrame(findOverlaps(benchmarkDataOne,PlotRange))
  benchmarkDataOne <- benchmarkDataOne[unique(BenchFilter$queryHits)]
  
  # creating GenomicInteractions
  
  if (length(benchmarkDataOne)!=0){   
    
    benchmarkDataOne <- GenomicInteractions(first(benchmarkDataOne), 
                                            second(benchmarkDataOne),
                                            counts=1)
    
    
    Benchmark_track <- InteractionTrack(benchmarkDataOne, 
                                        name = name, 
                                        chromosome = chr)
    if (is.null(name)) {
      
      Benchmark_track <- InteractionTrack(benchmarkDataOne, 
                                          name="Benchmark", 
                                          chromosome = chr)
    }
    
    displayPars(Benchmark_track) = list(col.anchors.fill ="blue",
                                        col.anchors.line = "grey",
                                        interaction.dimension="height", 
                                        interaction.measure ="counts",
                                        plot.trans=FALSE,
                                        plot.outside = TRUE, 
                                        col.outside="grey", 
                                        anchor.height = 0.1)
    
    
    return(Benchmark_track)
  }}

#' Selects genes of interest and associated reg regions from GInteractions obj
#' 
#' @author IngaPa
#' @keywords internal
selectGeneF <- function(x,
                        selectGene){

  # eliminating other genes if selectGene chosen
    x <- x[x$name2%in%selectGene]

  return(x)

}


#' Selects reg regions of interest from GInteractions obj,and associated genes
#' however for associated genes it reports as well all other interactions 
#' that are associated with genes associated with queried reg region
#' 
#' @author IngaPa
#' @keywords internal
selectRegR <- function(x,
                       selectRegulatoryRegion){

  #x <- interactionData[[1]]
  regRegion <- GRanges(selectRegulatoryRegion)
  SelectedEPpairs <- DataFrame(findOverlaps(regRegion,x))
  
  # save original x such that you can retrieve interactions with associated gene
    x.original <- x
  
  # select genes associated with subselected region  
    x <- x[SelectedEPpairs$subjectHits]
    x$color <- "red"
    
  # retrieve all interactions for given genes such that context can be plotted  
    x.original <- x.original[x.original$name2%in%unique(x$name2)]
    x.original$color <- "grey"
  
  # get all interactions associated with genes associated with gene of interest  
    x.context <- unique(c(x,x.original))
    
    }


#' Selects column that will be used as an input to plot heigth of the loops,
#' eg "pval". However -log10(pvalue) is calculated, not original value
#' 
#' @author IngaPa
#' @keywords internal
getStatistics <- function(x,
                          statistics){

          if (!is.null(statistics)){
            
                  x <- x[complete.cases(mcols(x)[statistics])]
                  
                  counts <- -log10(unlist(mcols(x)[statistics]))
            
          }

          if (is.null(statistics)){counts <- rep(1,length(x))}
  
    x$counts <- counts
  
    return(x)
  
}

#' Adjusts for colors of plot: three options:
#' 1) all red, in not info provided
#' 2) uses user-provided color vector (eg "color" column from GInteractions
#' object)
#' 3) If user selects region of interest, this function still plots other 
#' interactions around, but only in grey, whereas identified interaction
#' supposed to be red
#' @author IngaPa
#' @keywords internal
setColors <- function(x,
                      coloring) {
      # select colors from added column in interaction data input object      
      if (!is.null(coloring)) { 
        coloring <- unlist(mcols(x)[coloring])
      }             
      # interaction colors      
      if (is.null(coloring)) { 
        coloring <-  rep("red", length(x))
      }
    x$coloring <- coloring
}


#' Help function to create InteractionTrack
#' 
#' @author IngaPa
#' @keywords internal
createIntTrack <- function(x,
                           interactionDataName,
                           chr,
                           coloring){
  
  
      # creating GenomicInteractions
      EnhPromGInteractions <- GenomicInteractions(first(x), 
                                                  second(x),
                                                  x$counts)      

      # creating interaction track
      if (!is.null(interactionDataName)){
        
      interaction_track <- InteractionTrack(EnhPromGInteractions, 
                                            name =  interactionDataName, 
                                            chromosome = chr)
      }
      
      
      if ("listName"%in%names(mcols(x))){
        
        interaction_track <- InteractionTrack(EnhPromGInteractions, 
                                              name =  unique(x$listName), 
                                              chromosome = chr)
      }
      
      # get info about loop color     
      color <- setColors(x=x,
                         coloring=coloring)
      
      # set display parameters for track
      displayPars(interaction_track) = list(col.interactions=color, 
                                            col.anchors.fill ="blue",
                                            col.anchors.line = "black",
                                            interaction.dimension="height", 
                                            interaction.measure ="counts",
                                            plot.trans=FALSE,
                                            plot.outside = TRUE, 
                                            col.outside="grey", 
                                            anchor.height = 0.1)
      
      return(interaction_track)
      
      
}

#' Help function to get info about exon location for genes of interest from
#' Ensembl
#' 
#' filters <- c("hgnc_symbol","ensembl_gene_id")
#' 
#' @author IngaPa
#' @keywords internal
gettingGeneInfoEnsembl <- function(GeneList,
                                   filters){
  
      require(GenomicFeatures)
      require(TxDb.Hsapiens.UCSC.hg19.knownGene)
      require(biomaRt)
      
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
      
      GRList <- exonsBy(txdb, by = "gene") #exon coordinates
      
      # getting entrez IDs -> common names (since entrez is input)
      
      mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      dataset="hsapiens_gene_ensembl",
                      host="grch37.ensembl.org")
      
      
      EntrezIDs <- getBM(attributes=c("hgnc_symbol","entrezgene",
                                      "ensembl_gene_id"),
                         filters=filters,
                         values=GeneList,mart=mart)
  
      # selecting exons of interest      
      GRList <- GRList[names(GRList)%in%EntrezIDs$entrezgene]
      
      # adding corresponding hgnc_symbol
      names(GRList) <-  EntrezIDs$hgnc_symbol[match(names(GRList),
                                                    EntrezIDs$entrezgene)]
      
      # I will use GRanges as an input - easiest
      
      exonGeneObj <- unlist(GRList)
      # adding necassary meta-data
      exonGeneObj$symbol <- exonGeneObj$gene <- 
        exonGeneObj$transcript <-   names(exonGeneObj)
      exonGeneObj$feature= "protein coding"

      return(exonGeneObj)
      
      
}  





#' My stupid help function to get TSS for lists of enhancer~gene interactions 
#' sets
#' 
#' @author IngaPa
#' @keywords internal
getTSS <- function(AllGenes,interactionData) {
  
  return(do.call(c,do.call(c,lapply(AllGenes,function(x){
    
    return(unique(sapply(interactionData,function(y){
      
      return(unique(second(y)[y$name2==x]))})))}))))
  
}