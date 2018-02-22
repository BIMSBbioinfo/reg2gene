#' Plots enhancer-promoter interactions
#' 
#' 
#' @param  interactionData A GInteractions object with anchor1 corresponding to 
#' regulatory region, and anchor2 to the TSS of a gene. Minimum meta-data is 
#' info about gene names store in column "name2".
#' This object can be produced by modelling [ \code{\link{associateReg2Gene}}], 
#' or meta-analysis [\code{\link{metaAssociations}}] or voting functions from 
#' reg2gene package. 
#' 
#' @param  range (default NULL) this object stores info about gene annotations.
#' If default, function retrieves automatically gene annotations based on
#' TxDb.Hsapiens.UCSC.hg19.knownGene package.
#' Otherwise, this argument corresponds to the range argument from
#' \code{\link[Gviz]{GeneRegionTrack}} which handles many different data input
#' types:  TxDb, GRanges, GRangesList,data.frame, character scala. In the case
#' you want to write your own documentation, for more details about this 
#' argument take a look at \code{\link[Gviz]{GeneRegionTrack}}.
#' 
#' @param selectGene (default NULL) a character vector of gene name symbols 
#' (eg "FTO","IRX3", etc.) which user wants to plot.
#' 
#' @param selectRegulatoryRegion (default NULL) GRanges object or a 
#' character vector which stores info about the region of interest which 
#' user wants to plot (format: "chr16:53112601-53114200"). This region does not
#' necessarilly need to be equal to the regulatory regions reported in the 
#' interactionData input objects, whereas it only needs to overlap some
#' regulatory regions
#' @param benchmarkData (default NULL) A GInteractions object or a list of
#' GInteractions objects storing info about interaction regions in the genome 
#' (for example, regions from the literature obtained using chromatin 
#' conformation capture related methods.
#' 
#' @param statistics (default NULL) Column name where info about the heights of 
#' loops is stored. If NULL, then all interaction loops have an equal height 
#' Otherwise, loops of different height are plotted based on the info stored in
#' the corresponding column, eg "pval" or "qval" 
#' 
#' @param coloring (default NULL) Information about colors used to plot enhancer
#' promoter interactions. If it is NULL, then all interactions will be plotted
#' red. Otherwise, an ID of the column from interactionData object where info 
#' about color is stored. For example, interactions for different genes can be 
#' colored differently, or statistically significant/insignificant interactions
#' can be colored differently. User pre-defines this column by itself.
#'
#' @import GenomicFeatures
#' @import Gviz
#' @import GenomicInteractions 
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import biomaRt
#'    
#' @return Nice plot
#' 
#' @details  Plots InteractionTrack, GeneRegionTrack, AnnotationTrack, 
#' IdeogramTrack and GenomeAxisTrack for requested GInteractions object of
#' enhancer-gene associations 
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
#' 
#' # Creating an example GInteractions set:
#' 
#' enhancers <- GRanges(rep("chr16",6),
#'                      IRanges(c(53112601,55531601,53777201,
#'                                53778801,54084001,53946467),
#'                              c(53114200,55533250, 53778800,
#'                                53780400, 54084400 ,53947933)))
#' 
#' genes <- GRanges(rep("chr16",6),
#'                  IRanges(c(53737874, 54964773, 54320676,
#'                            53737874, 54964773, 54320676),
#'                          c(53737874, 54964773, 54320676,
#'                            53737874, 54964773, 54320676)))
#' 
#' GenomeInteractions <- GInteractions(enhancers,genes)
#' 
#' GenomeInteractions$name2 <- c("FTO","IRX5","IRX3")
#' 
#' GenomeInteractions$pval <- c(0.20857403, 0.72856090, 0.03586015,
#'                              0.32663439, 0.32534945, 0.03994488)
#' 
#' GenomeInteractions$color <- c("red","blue","grey")
#' 
#' plotGenomeInteractions(interactionData = GenomeInteractions,
#'                  statistics ="pval",
#'                  coloring = "color")
#' 
#' # if no info about coloring and height of loops provided, all loops are equal
#' plotGenomeInteractions(interactionData = GenomeInteractions)
#'  
#'  # Specific gene can be individually plotted
#'  
#'  plotGenomeInteractions(interactionData = GenomeInteractions,
#'  selectGene = "FTO")
#'  
#'  # This should be equal to the example where genes are not selected
#'  
#'  plotGenomeInteractions(interactionData = GenomeInteractions,
#'                         selectGene = c("FTO","IRX3","IRX5"),
#'                         coloring = "color",
#'                         statistics = "pval")
#'                         
#'  # if one wants to plot regulatory region                       
#'  plotGenomeInteractions(interactionData = GenomeInteractions,
#'               selectRegulatoryRegion = "chr16:53112601-53114200")
#'  
#'  # and function plots not regulatory regions queried, but all that
#'  # are present in the interactionData and overlap that region
#'                             
#'  plotGenomeInteractions(interactionData = GenomeInteractions,
#'                         selectRegulatoryRegion = "chr16:53112601-53778800",
#'                         coloring = "color",
#'                         statistics = "pval")
#'                         
#'  # add bemchmark data
#'  
#'  plotGenomeInteractions(interactionData = GenomeInteractions, 
#'  coloring = "color", statistics = "pval",
#'  benchmarkData = GenomeInteractions[1:3])
#'  
#'   # add a list of bemchmark data  
#'    benchmarkData = list(GenomeInteractions[1:3],GenomeInteractions[4:5],
#'    GenomeInteractions[6])
#'    names(benchmarkData) <- c("Bench1","Bench2","Bench3")

#'   plotGenomeInteractions(interactionData = GenomeInteractions,
#'   coloring = "color",
#'   statistics = "pval",
#'   benchmarkData = benchmarkData)                     
#'                                                                       
#' @export
plotGenomeInteractions <- function(interactionData,
                                   range=NULL,
                                   selectGene=NULL,
                                   selectRegulatoryRegion=NULL,
                                   benchmarkData=NULL,
                                   statistics=NULL,
                                   coloring=NULL){
  
  require(Gviz)
  require(GenomicInteractions)
  
  # extracting info about chromosome & genome set to hg19
  chr <- as.character(unique(seqnames(first(interactionData))))
  gen <- "hg19"
  
  # extracting info about genes
  if (!is.null(selectGene)){
    GeneList <- selectGene
    # eliminating other genes if selectGene chosen
    interactionData <- interactionData[interactionData$name2%in%
                                         selectGene]
  }
  if (is.null(selectGene)){GeneList <- unique(interactionData$name2)}
  
  # extracting info about regulatory regions of interest
  
  if (!is.null(selectRegulatoryRegion)){
    
    regRegion <- GRanges(selectRegulatoryRegion)
    SelectedEPpairs <- DataFrame(findOverlaps(regRegion,
                                              interactionData))
    
    interactionData <- interactionData[SelectedEPpairs$subjectHits]
    GeneList <- unique(interactionData$name2)
    
    
  }
  
  
  
  
  # getting height of loops removing NA values 
  
  if (!is.null(statistics)){
    
    interactionData <- interactionData[
      complete.cases(mcols(interactionData)[statistics])]
    counts <- -log10(unlist(mcols(interactionData)[statistics]))
    
  }
  
  if (is.null(statistics)){ counts <- rep(1,length(interactionData))}
  
  
  # extracting info about enhancers&genes
  Enhancers <- first(interactionData)
  Promoters <- second(interactionData)
  
  
  # creating GenomicInteractions
  EnhPromGInteractions <- GenomicInteractions(Enhancers, 
                                              Promoters,
                                              counts)
  
  
  # SETTING PLOTTING OPTIONS:
  
  # select colors from added column in interaction data input object      
  if (!is.null(coloring)) { 
    coloring <- unlist(mcols(interactionData)[coloring])
  }             
  # interaction colors      
  if (is.null(coloring)) { 
    coloring <-  rep("red", length(EnhPromGInteractions))
  }
  
  
  
  # creating interaction track
  interaction_track <- InteractionTrack(EnhPromGInteractions, 
                                        name = "Enh~Prom Interactions", 
                                        chromosome = chr)
  
  
  # creating benchmark track for given region of interest
  
  
  
  # if no info provided about genes - retrieve info automatically
  
  if (is.null(range)){ 
    
    require(GenomicFeatures)
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    require(biomaRt)
    
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
    
    GRList <- exonsBy(txdb, by = "gene") #exon coordinates
    
    # getting entrez IDs -> common names (since entrez is input)
    
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    dataset="hsapiens_gene_ensembl",
                    host="grch37.ensembl.org")
    
    EntrezIDs <- getBM(attributes=c("hgnc_symbol","entrezgene"),
                       filters="hgnc_symbol",
                       values=GeneList,mart=mart)
    
    # selecting exons of interes      
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
    
    
    #   setting ranges for plot based on enhancer & gene span 
    # identify min and
    # max range
    
    Span <- range(c(start(c(first(interactionData),
                            second(interactionData))),
                    start(exonGeneObj),
                    end(c(first(interactionData),
                          second(interactionData))),
                    end(exonGeneObj)))
    
    # extend min&max for 10000
    PlotRange <- GRanges(chr,IRanges(Span[1],Span[2]))
    PlotRange <- resize(PlotRange,
                        fix = "center",
                        width=width(PlotRange)+10000)    
    
    
    
    grtrack <- GeneRegionTrack(exonGeneObj,
                               genome = gen,
                               chromosome = chr,
                               name = "Gene symbol",
                               transcriptAnnotation="symbol",
                               # feature = "protein coding",
                               stacking="full",
                               labelPos="above",
                               start = start(PlotRange), 
                               end = end(PlotRange)) 
    
    
    
    
    
    
  }
  if (!is.null(range)){ 
    
    grtrack <- GeneRegionTrack(range=range,
                               genome = gen,
                               chromosome = chr,
                               name = "Genes",
                               transcriptAnnotation="symbol",
                               stacking="full",
                               labelPos="above")
    
    
    #   setting ranges for plot based on enhancer span
    
    Span <- range(c(start(c(first(interactionData),
                            second(interactionData))),
                    end(c(first(interactionData),
                          second(interactionData)))))
    
    # extend min&max for 10000
    PlotRange <- GRanges(chr,IRanges(Span[1],Span[2]))
    PlotRange <- resize(PlotRange,
                        fix = "center",
                        width=width(PlotRange)+10000)   
    
  }
  
  # setting promoter track
  # promoterTrack <- AnnotationTrack(unique(second(interactionData)), 
  #                                  genome=gen, 
  #                                  name="Promoters",
  #                                  id=unique(interactionData$name2),
  #                                  stacking="full")
  
  
  
  # setting enhancer track
  
  enhTrack <- AnnotationTrack(Enhancers, 
                              genome=gen, 
                              name="Assoc enhancers",
                              stacking="full")
  
  # setting chromosome & location track
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = gen, 
                          chromosome = chr)
  
  # setting display parameters
  
  displayPars(enhTrack) <- list(fill = "blue", 
                                col = "blue", 
                                just.group="above",
                                cex.feature = 0.7,
                                labelPos="above")
  
  
  
  displayPars(interaction_track) = list(
    col.interactions=coloring, 
    col.anchors.fill ="blue",
    col.anchors.line = "black",
    interaction.dimension="height", 
    interaction.measure ="counts",
    plot.trans=FALSE,
    plot.outside = TRUE, 
    col.outside="grey", 
    anchor.height = 0.1)
  
  
  
  displayPars(grtrack) <- list(fill = "grey",
                               col = "lightgrey", 
                               cex.feature = 10,
                               labelPos="above")
  
  if (is.null(benchmarkData)){ 
    
    plotTracks(list(itrack,
                    gtrack, 
                    grtrack,
                    #promoterTrack, 
                    enhTrack,
                    interaction_track),
               chromosome=chr, 
               from=min(start(PlotRange)), 
               to=max(end(PlotRange)), 
               sizes=c(0.1,0.2 ,0.2, 0.2,0.3),
               background.panel = "#FFFEDB", 
               background.title = "darkblue",
               fontcolor.feature = "darkblue")
  }
  
  # extra plot for benchmark
  
  if (!is.null(benchmarkData)){
    
    if (class(benchmarkData)=="GInteractions"){
      
      Benchmark_track <- benchPlotHelp(name=NULL,
                                       benchmarkDataOne=benchmarkData,
                                       PlotRange=PlotRange)
      
      
      
      displayPars(Benchmark_track) = list(col.anchors.fill ="blue",
                                          col.anchors.line = "grey",
                                          interaction.dimension="height", 
                                          interaction.measure ="counts",
                                          plot.trans=FALSE,
                                          plot.outside = TRUE, 
                                          col.outside="grey", 
                                          anchor.height = 0.1)
      
      plotTracks(list(itrack,
                      gtrack, 
                      grtrack,
                      enhTrack,
                      interaction_track,
                      Benchmark_track),
                 chromosome=chr, 
                 from=min(start(PlotRange)), 
                 to=max(end(PlotRange)), 
                 sizes=c(0.1,0.2 ,0.2, 0.2,0.3,0.3),
                 background.panel = "#FFFEDB", 
                 background.title = "darkblue",
                 fontcolor.feature = "darkblue")
    }
    
    if (class(benchmarkData)=="list"){
      
      
      Benchmark_track <- lapply(names(benchmarkData),
                                function(x){
                                  benchPlotHelp(name=x,
                                                benchmarkDataOne = benchmarkData,
                                                PlotRange=PlotRange)})
      
      
      plotTracks(c(list(itrack,
                        gtrack, 
                        grtrack,
                        enhTrack,
                        interaction_track),
                   Benchmark_track),
                 chromosome=chr, 
                 from=min(start(PlotRange)), 
                 to=max(end(PlotRange)), 
                 sizes=c(0.1,0.2 ,0.2, 0.2,0.3,
                         rep(0.2,length(Benchmark_track))),
                 background.panel = "#FFFEDB", 
                 background.title = "darkblue",
                 fontcolor.feature = "darkblue")
    }
    
    
  }
}




#' Plots enhancer-promoter interactions
#' 
#' Help function to create benckmark tracks, especilly needed for a list of 
#' benchmark inputs.
#' @author IngaPa
#' @keywords internal
benchPlotHelp <- function(name=NULL,
                          benchmarkDataOne,
                          PlotRange){
  
  if (!is.null(name)) {benchmarkDataOne <- benchmarkData[[name]]}
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






