#' Hierarchically annotates input GRanges object with corresponding genes
#'
#' This function either annotates input regions to their nearby genes or it 
#' hierarchical runs association procedure (when GIneraction object is used as
#' an input) as follows: promoters,enhancers, nearby genes.
#'
#' @param genomicRegions GRanges object that contains the windows of interest. 
#' It could be enhancers, promoters, CpG islands, ChIP-Seq or DNase-Seq peaks...
#' @param annotations GRanges which contains info about genes of interest, or
#' GInteractions object which stores regulatory regions - gene associations.
#' Info about gene names should be stored as a meta-data. 
#' IMPORTANT! 
#' If GInteractions object is used as an input then, anchor1 needs 
#' to be regulatory region, whereas anchor2 is location of gene/TSS.
#' If GRanges object is used as an input, then input genomicRegions are 
#' tested whether they are 
#' located in the vicinity of TSS of genes stored in this object, and if yes 
#' the corresponding gene is reported. Otherwise, hierarchical association
#' of input GRanges object with corresponding genes is performed.
#' @param upstream number of basepairs upstream from TSS to look for
#' input genomicRegions. default 1000
#' @param downstream number of basepairs downstream from TSS to look for 
#' input genomicRegions. default 1000
#' @param distance (default 1Mb). Maximal allowed distance between genes TSS and
#' input peak. Used in the 3rd step of the association procedure (genomic 
#' regions is associated to the closest gene if the distance between these 
#' two locations is smaller than prediefined distance threshold.)
#' @param identified (default TRUE). 
#' If TRUE, report genomic regions AND corresponding genes;
#' If FALSE, report regions for which info about gene is missing.
#' 
#' 
#' @details This function annotates input genomicRegions to the promoter 
#' regions of genes (and correspondingly to that gene) if only location of genes
#' (as GRanges object) and of genomicRegions is provided. 
#' Meaning that, if the input genomicRegions is located within +/- 
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
#'
#'
#' @return A \code{\link[InteractionSet]{GInteractions}} object that contains
#' info about genes (location+meta-data) annotated with regions of interest.
#' Anchor1 corresponds to the queried genomicRegions location, whereas anchor2
#' corresponds to gene location.  It additionaly reports whether genomic region
#' overlapped with promoter or enhancer region of the annotated gene, or an
#' association was assessed based on the proximity to the gene (possible values
#' "enhancer","promoter","nearestGene").
#' 
#' @import GenomicRanges
#' @import InteractionSet
#' @import genomation
#' 
#' @examples # CREATE genomicRegions test object
#' 
#' library(GenomicRanges)
#' library(InteractionSet)
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
#' annotateGenomicRegions(genomicRegions,
#'                     annotations =annotations,
#'                     identified=T)
#' annotateGenomicRegions(genomicRegions,
#'                     annotations,
#'                     identified=F)
#' 
#' @export
annotateGenomicRegions <- function(genomicRegions,
                      annotations,
                      upstream=1000,
                      downstream=1000,
                      distance=1000000,
                      identified=TRUE,
                      ...){
              
              # step 0
              # assign to closest TSS when annotations are GRanges
              
          if (class(annotations)=="GRanges"){
                
                # extending genes & enhancers +/- upstream/downstream
                
                      TSSextended <- trim(suppressWarnings(promoters(
                                            annotations,
                                            upstream,
                                            downstream)))
                
                # identify promoter-TSS combinations
                      PromoterAssociated <- DataFrame(
                                        findOverlaps(genomicRegions,TSSextended)
                                                      )
                
                      PromotersPeak <- GInteractions(genomicRegions[
                                        PromoterAssociated$queryHits],
                                    TSSextended[PromoterAssociated$subjectHits])
                # adding meta-data
                      mcols(PromotersPeak) <- mcols(annotations[
                                                PromoterAssociated$subjectHits])
                      PromotersPeak$annotatedAs <- "promoter"
                # return results: explained vs unexplained genomicRegions     
                        if (!identified) {
                          return(genomicRegions[-(unique(
                                                PromoterAssociated$queryHits))])
                        }
                        
                        if (identified) {
                          return(PromotersPeak)
                        }
                        
                        
                         }
              
          if (class(annotations)=="GInteractions"){
                
                # extending genes & enhancers +/- upstream/downstream
              
                      TSSextended <- trim(suppressWarnings(promoters(
                                                      second(annotations),
                                                      upstream,
                                                      downstream)))
                      rRegionextended <- trim(suppressWarnings(
                                                  promoters(first(annotations),
                                                  upstream,
                                                  downstream)))
                        
                 # identify promoter-TSS combinations
                      PromoterAssociated <- DataFrame(
                                      findOverlaps(genomicRegions, TSSextended)
                                                      )
                      if (nrow(PromoterAssociated)!=0){ 
                        # if something is associated with promoter - proceed
                      
                      PromotersPeak <- GInteractions(
                                  genomicRegions[PromoterAssociated$queryHits],
                                  TSSextended[PromoterAssociated$subjectHits])
                      
                      mcols(PromotersPeak) <- mcols(annotations[
                                                PromoterAssociated$subjectHits])
                        
                          PromotersPeak$annotatedAs <- "promoter"
                      }
                
                # for remaining genomicRegions:
                #-----------------
                
                # 1. try to explain as an enhancer
                        # identify remaining unexplained genomicRegions
                        
                       
                      genomicRegionsEnhancer <- genomicRegions[
                                       -(unique(PromoterAssociated$queryHits))]
                        
                        if (length(genomicRegionsEnhancer)!=0){ 
                          
                          # is there anything that can be associated with enh
                        
                          # overlap with regulatory regions
                          EnhancerAssociated <- DataFrame(
                            findOverlaps(genomicRegionsEnhancer,rRegionextended)
                                                    )
                          if (nrow(EnhancerAssociated)!=0){
                            
                            # is there anything that is associated with enh
                            
                        # identify corresponding genes
                          EnhancerPeak <- GInteractions(
                           genomicRegionsEnhancer[EnhancerAssociated$queryHits],
                                    TSSextended[EnhancerAssociated$subjectHits])
                        
                          mcols(EnhancerPeak) <- mcols(annotations[
                                              EnhancerAssociated$subjectHits])
                        
                            EnhancerPeak$annotatedAs <- "enhancer"
                        
                          }
                        
                        }
                      
                      
                # 2. or assign to the closest gene within +/- 1000000
                        
                        #  identify remaining unexplained genomicRegions
                      genomicRegionsClosestG <- genomicRegionsEnhancer[-(unique(
                                                 EnhancerAssociated$queryHits))]
                        
                      
                        if (length(genomicRegionsClosestG)!=0){
                          # is there anything that can be associated with 
                          # nearest gene
                           
                            nearestGenes <- DataFrame(
                           distanceToNearest(genomicRegionsClosestG,TSSextended)
                                  
                                         )
                      
                            nearestGenesAss <- nearestGenes[
                              which(nearestGenes$distance<distance),]
                            
                        if (nrow(nearestGenesAss)!=0){        
                        # threshold nearest genes by distance
                           
                          # genes 
                            NearestGenePeak <- GInteractions(
                              genomicRegionsClosestG[nearestGenesAss$queryHits],
                                       TSSextended[nearestGenesAss$subjectHits])
                        
                            mcols(NearestGenePeak) <- mcols(annotations[
                                                  nearestGenesAss$subjectHits])  
                      
                            NearestGenePeak$annotatedAs <- "nearestGene"
                        }  
                             
                # 3. unassign if none of these categories is matched
                
                    RemainGenes <- nearestGenes[which(nearestGenes$distance>
                                                        distance),]
                  
                        }
                  # return genomicRegions that remained unexplained
                  
                       
                      
                      # what to report
                        if (!identified) {
                          
                          return(genomicRegionsClosestG[RemainGenes$queryHits])
                              
                          }
                  
                        if (identified) {
                        
                              return(c(PromotersPeak,
                                       EnhancerPeak,
                                       NearestGenePeak))
                              }
                
                }
              
              
}
  
  
  