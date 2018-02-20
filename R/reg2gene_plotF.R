#' Plots enhancer-promoter interactions
#' 
#' 
#' @param  interactionData A GInteractions object. It can be produced by 
#' modelling [ \code{\link{associateReg2Gene}}], or
#' meta-analysis [\code{\link{metaAssociations}}] or voting functions from 
#' reg2gene. It can be GInteractions without any additional meta-data, but the
#' first element of input associating pair corresponds to regulatory region, 
#' whereas the second memeber corresponds to TSS/gene region. 
#' 
#' @param  geneInfo a data-table with info about genes stored in the following 
#' format: "seqnames" ,"start","end","width","strand", "symbol", "transcript"
#' are required. "symbol" corresponds to gene name, whereas "transcript" to info
#' about transcripts for given gene.
#' 
#' @param benchmarkData A GInteractions object with regions from the literature
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
#' @return Nice plot
#' 
#' @details  Plots InteractionTrack, GeneRegionTrack, AnnotationTrack, 
#' IdeogramTrack and GenomeAxisTrack for requested GInteractions object of
#' enhancer-gene associations 
#' 
#' @author Inga Patarcic
#' 
#' @examples
#' library(Gviz)
#' library(GenomicInteractions)
#' 
#'  \dontrun{ 
#' 
#' interactionData.ss <- interactionData[interactionData$name2=="FTO"]
#' geneInfo.ss <- geneInfo[geneInfo$gene.name=="FTO",]
#' plotAssociations(interactionData.ss,
#'                  geneInfo.ss,
#'                  statistics ="pval")
#' 
#' 
#' interactionData.ss <- interactionData[interactionData$name2%in%c("FTO","IRX3")]
#' interactionData.ss$color <- rep(c("red","blue"),table(interactionData.ss$name2))  
#' 
#' 
#' interactionData.ss <- interactionData[interactionData$name2%in%c("FTO","IRX3","IRX5")]
#' interactionData.ss$color <- rep(c("red","blue","grey"),
#'                                 table(interactionData.ss$name2)[
#'                                   unique(interactionData.ss$name2)]) 
#' geneInfo.ss <- geneInfo[geneInfo$gene.name==c("FTO","IRX3","IRX5"),]
#' 
#' 
#' plotAssociations(interactionData.ss,
#'                  geneInfo.ss,
#'                  statistics ="pval",
#'                  coloring = "color")
#' 
#' 
#' interactionData.ss <- interactionData.ss[which(interactionData.ss$pval<0.05)]  
#' 
#' plotAssociations(interactionData.ss,
#'                  geneInfo.ss,
#'                  statistics ="pval",
#'                  coloring = "color")
#' } 
#' 
#' @export
plotAssociations <- function(interactionData,
                             geneInfo,
                             benchmarkData=NULL,
                             statistics=NULL,
                             coloring=NULL){
                               
           
           
           # extracting info about chromosome & genome set to hg19
           chr <- as.character(unique(seqnames(first(interactionData))))
           gen <- "hg19"
           
           
           #   setting ranges for plot based on enhancer span
           
           Span <- range(c(start(c(first(interactionData),
                                   second(interactionData))),
                           end(c(first(interactionData),
                                 second(interactionData)))))
           
           
           PlotRange <- GRanges(chr,IRanges(Span[1],Span[2]))
           PlotRange <- resize(PlotRange,
                               fix = "center",
                               width=width(PlotRange)+10000)
           
           
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
           
         
           
           # gene track
           grtrack <- GeneRegionTrack(unique(as.data.frame(geneInfo)),
                                      genome = gen,
                                      chromosome = chr,
                                      name = "Genes",
                                      transcriptAnnotation="symbol",
                                      stacking="full",
                                      labelPos="above")
           
           
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
             
             # filter for regions of interest based on ranges of plot
             
             BenchFilter <- DataFrame(findOverlaps(benchmarkData,PlotRange))
             benchmarkData <- benchmarkData[unique(BenchFilter$queryHits)]
             
             # creating GenomicInteractions
             
             if (length(benchmarkData)!=0){   
               
               benchmarkData <- GenomicInteractions(first(benchmarkData), 
                                                    second(benchmarkData),
                                                    counts=1)
               
               
               Benchmark_track <- InteractionTrack(benchmarkData, 
                                               name = "Benchmark Interactions", 
                                                   chromosome = chr)
               
               
               displayPars(Benchmark_track) = list(
                 col.interactions=coloring, 
                 col.anchors.fill ="blue",
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
           }
           
           
           
         }

         
         
#' Plots enhancer-promoter interactions for genes of interest
#' 
#' Wrapper for plotAssociations(), such that individual genes can be requested
#' 
#' @param gene character vector of gene names.
#' 
#' @param  interactionData A GInteractions object. It can be produced by 
#' modelling [ \code{\link{associateReg2Gene}}], or
#' meta-analysis [\code{\link{metaAssociations}}] or voting functions from 
#' reg2gene. It can be GInteractions without any additional meta-data, but the
#' first element of input associating pair corresponds to regulatory region, 
#' whereas the second memeber corresponds to TSS/gene region. 
#' 
#' @param  geneInfo a data-table with info about genes stored in the following 
#' format: "seqnames" ,"start","end","width","strand", "symbol", "transcript"
#' are required. "symbol" corresponds to gene name, whereas "transcript" to info
#' about transcripts for given gene.
#' 
#' @param benchmarkData A GInteractions object with regions from the literature
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
#' @return Nice plot for genes of interest
#' 
#' @details  Plots InteractionTrack, GeneRegionTrack, AnnotationTrack, 
#' IdeogramTrack and GenomeAxisTrack for requested GInteractions object of
#' enhancer-gene associations 
#' 
#' @author Inga Patarcic
#' 
#' @examples
#' \dontrun{
#' library(Gviz)
#' library(GenomicInteractions)
#' 
#' plotGene(gene=c("FTO"),
#' interactionData,
#' geneInfo,
#' statistics="pval",
#' coloring=NULL)
#' 
#' plotGene(gene=c("HBA1","HBA2"),
#'          interactionData,
#'          geneInfo,
#'          statistics="pval",
#'          coloring=NULL)
#' }
#' @export
plotGene <- function(gene,
                     interactionData,
                     geneInfo,
                     benchmarkData,
                     statistics,
                     coloring){
  
  require(stringr)
  # selecting genes
    # allow ensembl ID's and common gene ID's
        if (!any(str_detect(gene,"ENSG"))) {
        
          interactionData <- interactionData[interactionData$name2%in%gene]
            geneInfo <- geneInfo[geneInfo$symbol%in%gene,]
         
                   }  
        
        if (any(str_detect(gene,"ENSG"))) {
          
          interactionData <- interactionData[interactionData$name%in%gene]
          geneInfo <- geneInfo[geneInfo$symbol%in%gene,]
          
        }  
  
  # plot genes
  plotAssociations(interactionData,
              geneInfo,
              benchmarkData,
              statistics,
              coloring)

}



#' Plots enhancer-promoter interactions for reagulatory region of interest
#' 
#' Wrapper for plotAssociations(), such that individual regulatory regions
#' can be requested
#' 
#' @param regRegion GRanges object or a character vector which stores info about
#' the region of interest (format: "chr16:53112601-53114200").
#' 
#' @param  interactionData A GInteractions object. It can be produced by 
#' modelling [ \code{\link{associateReg2Gene}}], or
#' meta-analysis [\code{\link{metaAssociations}}] or voting functions from 
#' reg2gene. It can be GInteractions without any additional meta-data, but the
#' first element of input associating pair corresponds to regulatory region, 
#' whereas the second memeber corresponds to TSS/gene region.
#' 
#' @param  geneInfo a data-table with info about genes stored in the following 
#' format: "seqnames" ,"start","end","width","strand", "symbol", "transcript"
#' are required. "symbol" corresponds to gene name, whereas "transcript" to info
#' about transcripts for given gene.
#' 
#' @param benchmarkData A GInteractions object with regions from the literature 
#' 
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
#' 
#' @return Nice plot for genes of interest
#' 
#' @details  Plots InteractionTrack, GeneRegionTrack, AnnotationTrack, 
#' IdeogramTrack and GenomeAxisTrack for requested GInteractions object of
#' enhancer-gene associations. 
#' IMPORTANT! It is good to reduce the number of plotted interactions a priori, 
#' since many interactions can be tested for regRegion of interest, eg filter
#' based on p-value results.
#' 
#' @author Inga Patarcic
#' 
#' @examples
#' 
#' \dontrun{
#' library(Gviz)
#' library(GenomicInteractions)
#' interactionData <- readRDS("/data/akalin/Projects/AAkalin_reg2gene/Results/associateReg2Gene/PearsonRoadmaDNAme_forPlottingExample.rds")
#' geneInfo <- readRDS("/data/akalin/Base/Annotation/hg19/GENCODE/v24/gencode.v24lift37.basicannotationAndNoncodingGRanges170202_adjustedPlotting.rds")
#' regRegion <- "chr16:53112601-53114200"
#' interactionData <- interactionData[interactionData$pval<0.1]
#' 
#' plotRegulatoryRegion(regRegion,
#'                     interactionData,
#'                     geneInfo)
#' }
#' 
#' @export
plotRegulatoryRegion <- function(regRegion,
                                 interactionData,
                                 geneInfo,
                                 benchmarkData,
                                 statistics,
                                 coloring){
  
  regRegion <- GRanges(regRegion)
  SelectedEPpairs <- DataFrame(findOverlaps(regRegion,interactionData))
  
  interactionData <- interactionData[SelectedEPpairs$subjectHits]
  geneInfo <- geneInfo[geneInfo$symbol%in%interactionData$name2,]
  
  # plot genes
  plotAssociations(interactionData,
                   geneInfo,
                   benchmarkData)
  
}





