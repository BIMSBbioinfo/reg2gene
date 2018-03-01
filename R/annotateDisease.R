#' Add associated traits to the interaction object
#' 
#' This function adds associated traits to the input interaction object. 
#' Traits are extracted from the object which stores either variant-gene-disease
#' associations of gene-disease associations. Traits are idenfied based on 
#' based on gene symbol overlap or or based on variant-genomic region overlap in
#' addition to the gene symbol overlap.
#' 
#' @param  interactionData A GInteractions object. It can be produced by 
#' modelling [ \code{\link{associateReg2Gene}}], or
#' meta-analysis [\code{\link{metaAssociations}}] or voting functions from 
#' reg2gene. It can be GInteractions object with only gene symbols as an 
#' additional meta-data. Important is that the
#' first element of input associating pair corresponds to regulatory region, 
#' whereas the second memeber corresponds to TSS/gene region. 
#' 
#' @param  DisGeneDB a data-table or GRanges object which stores info about 
#' gene-disease associations (in data table format) or variant-gene-disease 
#' associations (as a GRanges object where ranges correspond to variant 
#' location). For example DISGENET db or GWAS Catalog db can be used as an 
#' input for this function. 
#' 
#' @param geneIDs character vector of length 2 which stores information
#' about column names where gene symbols are stored for both input objects. 
#' By default c("name2","GeneName")
#' 
#' 
#' @return A GInteractions object (interactionData) with added meta-data
#' from the DisGeneDB object. Added info is presumably info about traits
#' associated with given genes or both variant-genes combination.
#' 
#' @details  Function works with DisGeneDB object that can be either a
#' Granges object or dataframe storing info about gene~trait or
#' variant~gene~trait associations.
#' If input object is GRanges, then only variant (not gene) location should
#' be stored as a genomic ranges. In that case, function finds overlaps 
#' between anchor1 from GInteraction object and DisGeneDB genomic ranges
#' (presumably overlap between variant and regulatory regions), and in 
#' addition it identifies overlaps between genes from these two objects
#' based on gene symbols (stored in columns "name2","GeneName"). Thus,
#' only if gene symbols and genomic locations from both objects overlap
#' the function will return these entries.
#' 
#' If DisGeneDB stores info about gene~trait associations as a dataframe
#' (no info about variant location), then info about associated trait is
#' fetched only based on gene name overlap.
#' 
#' @author Inga Patarcic
#' 
#' @examples #[to do]
#' 
#' @import GenomicRanges
#' @import InteractionSet
#' 
#' @keywords internal
annotateRegionToTrait <- function(interactionData,
                                  DisGeneDB,
                                  geneIDs=c("name2","GeneName"),
                                  ...){
  
  # Finding info about genes when variant-gene-disease
  # association info available
  
  if (class(DisGeneDB)=="GRanges"){
    # expect GRanges or data-frame
    
    OverlapGI <-  DataFrame(findOverlaps(first(interactionData),
                                         DisGeneDB))
    
    # interactions for which overlap in disease-gene-variant 
    # database is identified
    interactionData <- interactionData[OverlapGI$queryHits]
    # interactions for which overlap in disease-gene-variant 
    # database is identified
    DisGeneDB <- DisGeneDB[OverlapGI$subjectHits]
    
    mcols(interactionData) <- DataFrame(mcols(interactionData),
                                        DisGeneDB)
    
    equalGenes <-  which(unlist(mcols(interactionData)[geneIDs[1]]==
                                  mcols(interactionData)[geneIDs[2]]))
    
    interactionData <- interactionData[equalGenes] 
    
    return(interactionData)
    
    
    
  }
  
  if (class(DisGeneDB)!="GRanges"){
    # expect GRanges or data-frame
    
    NameMatch <- unlist(match(mcols(interactionData)[geneIDs[1]],
                              DisGeneDB[geneIDs[2]]))
    
    if (!all(is.na(NameMatch))) {  # if any match identified
      
      # identify appropriate diseases 
      DisGeneDB <- DisGeneDB[NameMatch[complete.cases(NameMatch)],]
      
      # identify appropriate interactions objects which have available 
      # info about traits
      IndexForInterObj <- (1:length(interactionData))[complete.cases(NameMatch)]
      
      interactionData <- interactionData[IndexForInterObj]
      mcols(interactionData) <- c(mcols(interactionData),
                                  DisGeneDB )
    }
    
    
    return(interactionData)
    
    
    
  }
  
}


#' Compare two objects of genomic region~gene~trait associations
#'   
#' 
#' This function reports per trait number of reported 
#' genomic region~gene~trait from two different objects.
#' 
#' @param  reg2Trait A GInteractions object, likely produced by 
#' [\code{\link{annotateRegionToTrait}}]. This object stores info about
#' genomic region~gene~trait associations, as follows: anchor1 corresponds
#' to genomic region (likely SNP location), anchor2 corresponds to gene 
#' location, and gene symbol and associated disease is stored as a meta-data.
#' 
#' @param   reg2Trait2 A GInteractions object, likely produced by 
#' [\code{\link{annotateRegionToTrait}}]. This object stores info about
#' genomic region~gene~trait associations, as follows: anchor1 corresponds
#' to genomic region (likely SNP location), anchor2 corresponds to gene 
#' 
#' @param traitID character vector of length 2 which stores information
#' about column names where trait info are stored for both input objects. 
#' By default c("Trait","Trait")
#' 
#' 
#' @return a dataframe with rows corresponding to traits, and two columns with
#' info about gene counts for these traits. Each column gathers info from one 
#' reg2Trait GInteractions object (variant~gene interactions + trait 
#' association stored as meta-data).
#'  
#' @details For each trait that is present in at least one of the reg2Trait 
#' GInteractions objects report how many associated genes is present in the 
#' first GInteractions object, and how many is present in the 2nd 
#' GInteractions object. It allows an easy identification of traits for which
#' [\code{\link{annotateGenomicRegions}}] analysis improved number of disease
#' associated genes.
#' 
#' @author Inga Patarcic
#' 
#' @examples #[TO do]
#' 
#' @import GenomicRanges
#' @import InteractionSet
#' 
#' @keywords internal
compareReg2Traits <- function(reg2Trait,
                              reg2Trait2,
                              traitID=c("Trait","Trait")){
  
  # obtain per disease statistics
  reg2t.table <- table(mcols(reg2Trait)[traitID[1]])
  reg2t2.table <- table(mcols(reg2Trait2)[traitID[2]])
  
  # select diseases present in both datasets  
  Diseases <- names(reg2t.table)[names(reg2t.table)%in%names(reg2t2.table)]
  
  Compare2datasets <- cbind(reg2t.table[Diseases],reg2t2.table[Diseases])
  
  colnames(Compare2datasets) <- c("reg2Trait","reg2Trait2")
  
  # filtering when both entries have diseases with 0 gene counts
  
  Compare2datasets <- Compare2datasets[apply(Compare2datasets,1,
                                             function(x){return(!all(x==0))}),]
  
  return(Compare2datasets)
}




# reg2Trait <- annotateRegionToTrait(interactionData,DisGeneDB) 
# reg2Trait2 <-annotateRegionToTrait(interactionData,data.frame(DisGeneDB))
# 
# compareReg2Traits(reg2Trait,reg2Trait2)
# 
# 
# 
# 
# 
# gda <- read.delim("D:/Projects_Helping/PhD/reg2gene/dATA/curated_variant_disease_associations.tsv")
# 
# 
# 
# 
# # annotateRegionToDisease
# 
# 
# "diseaseName"
# "geneSymbol"
# 
# interactionData <- readRDS("D:/Projects_Helping/PhD/reg2gene/ModelResult.rds")
# 
# DisGeneDB <- read.delim("D:/Projects_Helping/PhD/reg2gene/dATA/curated_gene_disease_associations.tsv")
# 
# #1. GI-GI OVERLAP
# library(traseR)
# data("taSNP")
# DisGeneDB <- taSNP
# DisGeneDB$GeneName <- taSNP$GENE_NAME
# 
# 
# geneIDs=c("name2","GeneName")





# Finding info about genes when variant-gene-disease
# association info available

# if (class(DisGeneDB)=="GInteractions"){
# 
#  
#     OverlapGI <-  linkOverlaps(interactionData,
#                                first(DisGeneDB),
#                                second(DisGeneDB))
#     
#     # if subject1 overlaps subject1 then interaction is confirmed
#     benchRow <- which(OverlapGI$subject1==OverlapGI$subject2)
#      
#     # interactions for which overlap in disease-gene-variant 
#     # database is identified
#     interactionData <- interactionData[OverlapGI$query[benchRow]]
#     # interactions for which overlap in disease-gene-variant 
#     # database is identified
#     DisGeneDB <- DisGeneDB[OverlapGI$subject1[benchRow]]
#    
#    mcols(interactionData) <- DataFrame(mcols(interactionData),
#                                           DisGeneDB)
#  
#      return(interactionData)
#      
# }