#' Hierarchically annotates input GRanges object with corresponding genes
#'
#' This function either annotates input regions to their nearby genes or it 
#' hierarchical runs association procedure (when GIneraction object is used as
#' an input) as follows: promoters,enhancers, nearby genes.
#'
#' @param windows GRanges object that contains the windows (genomic regions) 
#' of interest:enhancers, promoters, CpG islands, ChIP-Seq or DNase-Seq peaks...
#' @param interactions (default NULL) GInteractions object which stores 
#' interactions of interest, either as regulatory regions - gene associations or 
#' as interacting locations obtained by chromatin conformation capture and
#' related methods. Info about gene names should be stored as a meta-data. 
#' Anchor1 needs to be regulatory region, whereas anchor2 is location of 
#' gene/TSS.
#' @param geneAnnotations GRanges which contains info about TSS for genes of 
#' interest, or ideally TSS for all genes in the genome, since nearest gene 
#' annotation methods (step 1 TSS+/-1000bp & step 3 TSS +/- 1Mb) is based
#' on the gene TSS coordinates from this object. 
#' @param annotateInteractions (default FALSE) This argument is useful in the 
#' case when interaction dataset does not store info about genes, whereas it 
#' provides only an info about interacting locations in the genome (not
#' necessarilly regulatory region - promoter interactions). If FALSE, no 
#' additional action done. If TRUE, interactions are firstly annotated
#' to the genes (unannotated interactions are removed from the dataset),
#' and then windows are annotated using annotated interaction object.
#' @param upstream number of basepairs upstream from TSS to look for
#' input windows. default 1000
#' @param downstream number of basepairs downstream from TSS to look for 
#' input windows. default 1000
#' @param distance (default 1Mb). Maximal allowed distance between genes TSS &
#' input peak. Used in the 3rd step of the association procedure (genomic 
#' regions is associated to the closest gene if the distance between these 
#' two locations is smaller than prediefined distance threshold.)
#' @param identified (default TRUE). 
#' If TRUE, report genomic regions AND corresponding genes;
#' If FALSE, report regions for which info about gene is missing.
#' NOTE! When this argument set to TRUE, regions are sometimes annotated based 
#' on the nearest gene procedure instead of enhancer. However, gene annotation
#' is correct.
#' 
#' @param ... Any parameter in the future
#' 
#' @details This function annotates input windows (genomicRegions) to the promoter 
#' regions of genes (and correspondingly to that gene) if only location of genes
#' (as GRanges object) and of windows (genomicRegions) is provided. 
#' Meaning that, if the input windows (genomicRegions) is located within +/- 
#' upstream/downstream distance from TSS of a gene, then this gene is annotated
#' to the queried region.
#' When GIneraction object is used as an input, then hierarchical association
#' procedure is runned as follows: promoters,enhancers, nearby genes, eg:
#' 1) genomic regions of interest are first considered to be 
#' promoters and annotated with nearby genes if they are located within a 
#' certain distance from TSS of nerbay gene (default +/-1000bp); otherwise
#' 2) remaning genes are overlapped with enhancer regions, and genes annotated
#' to that enhancer regions are reported, 
#' 3) if no overlap with either promoters nor enhancers is identified, then 
#' closest gene is reported if it is located within 1Mb
#' 4) if no gene located within 1Mb has been identified then, this region is 
#' filtered out.
#' IMPORTANT! 
#' Anchor1 in GInteractions object needs to be regulatory region, whereas
#' anchor2 is the location of gene/TSS.
#' If GInteractions is missing as an input, then input windows are 
#' tested ONLY whether they are located in the vicinity of TSS of genes stored
#' in this object, and if yes the corresponding gene is reported. 
#' 
#' @return A \code{\link[InteractionSet]{GInteractions}} object that contains
#' info about genes (location+meta-data) annotated with regions of interest.
#' Anchor1 corresponds to the queried windows (genomicRegions location),
#' whereas anchor2 corresponds to the extended (+/-upstream/downstream) 
#' TSS location. It additionaly reports whether genomic region was annotated 
#' based on the overlap with promoter or enhancer region of the annotated gene,
#' or an annotation was assessed based on the proximity to the gene (possible 
#' values "enhancer","promoter","nearestGene").
#' 
#' @import GenomicRanges
#' @import InteractionSet
#' @import genomation
#' 
#' @examples # CREATE genomicRegions test object - windows argument
#' 
#' library(GenomicRanges)
#' library(InteractionSet)
#' library(GenomicInteractions)
#' library(Gviz)
#' library(biomaRt)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(GenomicFeatures)
#' library(reg2gene)
#' 
#' genomicRegions <- GRanges(c("chr1:1-2", # 1. overlap prom
#'                             "chr2:1-2",  # 2. overlap enh
#'                             "chr3:1-2", # 3. overlap tss +/- 1,000,000
#'                             "chr4:1-2")) # 4. do not overlap tss +/- 1,000,000
#' 
#' # CREATE GInteractions test object
#' annotationsEnh <- GRanges(c("chr1:1-2",
#'                             "chr2:1-2",
#'                             "chr3:100000-100002",
#'                             "chr4:10000001-10000002"))
#' 
#' annotationsGenes <- GRanges(c("chr1:1-2",
#'                               "chr2:100000-100002",
#'                               "chr3:99999-100002",
#'                               "chr4:10000001-10000002"))
#' annotationsGenes$name=c("gen1","gen2","gen3","gen4")
#' seqlengths(annotationsEnh) <- seqlengths(annotationsGenes) <- rep(10000002,
#'                                                                   4)
#' 
#' 
#' annotations = GInteractions(annotationsEnh,annotationsGenes)
#' 
#' 
#' reg2gene(windows=genomicRegions,
#'                     annotations =annotations,
#'                     identified=TRUE)
#' reg2gene(windows=genomicRegions,
#'                     annotations,
#'                     identified=FALSE)
#' 
#' @export
reg2gene <- function(windows,
                    geneAnnotations,
                    interactions=NULL,
                    annotateInteractions=FALSE,
                    upstream=1000,
                    downstream=1000,
                    distance=1000000,
                    identified=TRUE,
                    ...){
  
              # extending genes & enhancers +/- upstream/downstream
          
            TSSextended <- trim(suppressWarnings(
                                promoters(geneAnnotations,upstream,downstream)))
                          
              # step 0
              # assign to closest TSS when geneAnnotations are GRanges
            
            PromotersAnnotated <- annotateToPromoter(windows,
                                                     TSSextended=TSSextended,
                                                geneAnnotations=geneAnnotations,
                                                     identified = TRUE)
            
            # UNassigned to the closest TSS based on geneAnnotations are GRanges
            UnAnnotatedWindows <- annotateToPromoter(windows,
                                                     TSSextended=TSSextended,
                                                geneAnnotations=geneAnnotations,
                                                     identified = FALSE)
              
                    # anntotated only to promoters
                     if (is.null(interactions)&(identified==TRUE)){
                       return(PromotersAnnotated)}
                    
                    # notanntotated but based only on promoters
                      if (is.null(interactions)&(identified==FALSE)){
                        return(UnAnnotatedWindows)}
            
              
        # for remaining windows:
        #-----------------
        # step2: annotate to enhancers
       
           if (!is.null(interactions)&(length(UnAnnotatedWindows)!=0)){
            
             # if only locations reported in interactions data, then 1st 
             # annotate them to genes, and then proceed with annotating windows
             
             if (annotateInteractions==TRUE) {
                  
                interactions <- annotateInteractionsToGenes(
                                                TSSextended=TSSextended,
                                                interactions = interactions)
                                             }
             
             # extending genes & enhancers +/- upstream/downstream
              
                  rRegionextended <- trim(suppressWarnings(
                                              promoters(first(interactions),
                                              upstream,
                                              downstream)))
                 
                        
                 EnhAss <- DataFrame(findOverlaps(UnAnnotatedWindows,
                                              rRegionextended))
                 
                 #  identify remaining unexplained windows
                 windowsClosestG <- UnAnnotatedWindows
                 
                       if (nrow(EnhAss)!=0){
                           
                        # identify corresponding genes
                           EnhancerPeak <- GInteractions(
                                          UnAnnotatedWindows[EnhAss$queryHits],
                                      second(interactions)[EnhAss$subjectHits])
                        
                                mcols(EnhancerPeak) <- mcols(interactions
                                                          [EnhAss$subjectHits])
                        
                                EnhancerPeak$annotatedAs <- "enhancer"
                                
                                windowsClosestG <- UnAnnotatedWindows[-(unique(
                                                            EnhAss$queryHits))]
                          }
                        
                        
      # step 3. or assign to the closest gene within +/- 1000000
                    
       
                    if (length(windowsClosestG)!=0){
                          # is there anything that can be associated with 
                          # nearest gene
                           
                      nearestGenes <- DataFrame(distanceToNearest(
                                                  windowsClosestG,TSSextended))
                      # if identified location is below threshold
                      identifiedGen <-  which(nearestGenes$distance<distance)
                      
                      # if certain locations were deleted due to the 
                      # distanceToNearest not taking into account different gene  
                      missingLocations <- which(!1:length(windowsClosestG)%in%
                                                    nearestGenes$queryHits)
                        
                      nearestGenesAss <- nearestGenes[identifiedGen,]
                            
                            if (nrow(nearestGenesAss)!=0){        
                        # threshold nearest genes by distance
                           
                            NearestGenePeak <- GInteractions(
                                  windowsClosestG[nearestGenesAss$queryHits],
                                  geneAnnotations[nearestGenesAss$subjectHits])
                        
                                colnames(mcols(NearestGenePeak)) <- "name"
                              NearestGenePeak$annotatedAs <- "nearestGene"
                        }  
                             
                # 3. unassign if none of these categories is matched
                
                    RemainGenes <- nearestGenes[which(nearestGenes$distance>
                                                        distance),]
                    # pooling locations that dropped and once filetered by 
                    # distance
                    
                    UnAnnotatedWindows <- c(windowsClosestG[missingLocations],
                                      windowsClosestG[RemainGenes$queryHits])
               
                        }
                  # return windows that remained unexplained
                     
          
      
      # step 4 pooling together info
            if (!identified) {return(UnAnnotatedWindows)}
                   
            if (identified) {
       
               obj <- list(
                    tryCatch(get("PromotersAnnotated"),error=function(e){NULL}),
                    tryCatch(get("NearestGenePeak"),error = function(e){NULL}),
                    tryCatch(get("EnhancerPeak"),error = function(e){NULL}))
        
        if (!is.null(unlist(obj))) {return(
                                        tryCatch(
                                          do.call("c",unlist(obj)),
                                            error=function(e) 
                      print("ERROR! Try setting annotateInteractions = TRUE")))}

         if (is.null(unlist(obj))) {print("No overlap identified!")}
                        
                                 }
                
                     }
              
              
            }
  

#' Annotate genomic location(windows) to genes reported in geneAnnotation 
#' 
#' Help f() to annotate loc1 (windows) to genes with which they potentially 
#' interact
#' 
#' @author Inga Patarcic
#' @keywords internal 

annotateToPromoter <- function(windows,
                               TSSextended,
                               geneAnnotations,
                               identified){

  
  # identify promoter-TSS combinations
  PromotersAnnotated <- DataFrame(findOverlaps(windows,TSSextended))
  
  if (nrow(PromotersAnnotated)!=0) {
  
    PromotersPeak <- GInteractions(windows[PromotersAnnotated$queryHits],
                                geneAnnotations[PromotersAnnotated$subjectHits])
  
  # adding meta-data
    mcols(PromotersPeak) <- mcols(geneAnnotations[PromotersAnnotated$subjectHits])
    
    PromotersPeak$annotatedAs <- "promoter"
      
  # return results: explained vs unexplained windows     
      if(!identified){return(windows[-(unique(PromotersAnnotated$queryHits))])}
      if(identified){return(PromotersPeak)}
  }

  if (nrow(PromotersAnnotated)==0) {
    
    if(identified){}
    if(!identified){return(windows)}
  
  
}



}


#' Annotate interaction set to genes reported in geneAnnotation 
#' 
#' Help f() to annotate loc2 | loc1 to genes with which they potentially 
#' interact
#' 
#' @author Inga Patarcic
#' @keywords internal 
annotateInteractionsToGenes <- function(TSSextended,
                                        interactions){
  
  
    return(unique(annotateToPromoter(windows = c(first(interactions),
                                   second(interactions)),
                       TSSextended=TSSextended,
                       geneAnnotations = geneAnnotations,
                       identified = TRUE)))
    
  
  
  }
