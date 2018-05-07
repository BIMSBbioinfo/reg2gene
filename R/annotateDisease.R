#' Add associated traits to the interaction object
#' 
#' This function adds associated traits to the input interaction object. 
#' Traits are extracted from the object which stores either variant-gene-disease
#' associations of gene-disease associations. Traits are idenfied based on 
#' based on gene symbol overlap or or based on variant-genomic region overlap in
#' addition to the gene symbol overlap.
#' 
#' @param  interactions A GInteractions or GRanges object with gene names as a
#' minimal metadata info. If GInteractions is inputed, the first
#' element corresponds to regulatory region, whereas the 
#' anchor anchor corresponds to TSS/gene region. 
#' 
#' @param  disGeneDB A GInteractions or GRanges object which stores info about 
#' gene-disease associations (GRanges object type) or variant-gene-disease 
#' associations (as a GInteractions object where anchor1 correspond to the 
#' variant location and anchor2 to gene location). 
#' For example, DISGENET db or GWAS Catalog db can be used as an 
#' input for this function. 
#' 
#' @param geneIDs character vector of length 2 which stores information
#' about column names where gene symbols are stored for both input objects. 
#' By default c("name2","GeneName")
#' 
#' 
#' @return A GInteractions object (interactions) with added meta-data
#' from the disGeneDB object. Added info is presumably info about traits
#' associated with provided gene names or variant-gene combination.
#' 
#' @details  Function works with disGeneDB object that can be either a
#' Granges or GInteractions object which stores info about gene~trait or
#' variant~gene~trait associations.
#' If input object is GInteractions, then function finds overlaps 
#' between anchor1 from GInteraction object and anchor1 of disGeneDB 
#' (presumably overlap between variant and regulatory regions), and in 
#' addition it identifies overlaps between genes from these two objects
#' based on gene symbols (stored in columns "name2","GeneName"). Thus,
#' only if gene symbols and genomic locations from both objects overlap
#' the function will return these entries.
#' 
#' If disGeneDB stores info about gene~trait associations as a GRanges object
#' (no info about variant location), then info about associated trait is
#' obtained only based on gene name overlap.
#' 
#' @author Inga Patarcic
#' 
#' @examples 
#' library(InteractionSet)
#' library(GenomicRanges)
#' library(reg2gene)
#' 
#' # 1. get disease-gene DB example:
#' 
#' GRReg2_toyDisease <- GRReg2_toy
#' GRReg2_toyDisease$disease <- paste0("dis",1:length(GRReg2_toy))
#' 
#' # 1.a get disease-gene DB example as GInteractions
#'  disGeneDB <- GInteractions(GRReg2_toyDisease,GRReg2_toyDisease$reg)
#'  mcols(disGeneDB) <- mcols(GRReg2_toyDisease)
#' 
#' 
#' # 1. get enhancer-gene:GInteractions 
#' 
#' GRReg1_toyGI <- GInteractions(anchor1 = GRReg1_toy,
#'                               anchor2 = GRReg1_toy$reg,
#'                               gene=GRReg1_toy$name)
#' 
#' mcols(GRReg1_toyGI) <- mcols(GRReg1_toyGI)[c("gene")]
#' 
#' 
#' 
#' # run reg2trait; version geneDisease
#' 
#'  reg2trait(interactions = GRReg1_toyGI,
#'                  disGeneDB = GRReg2_toyDisease,
#'                  geneIDs = c("gene","name"))
#'                  
#'                  
#'  # run reg2trait; version SNPgeneDisease       
#'           
#'  reg2trait(interactions = GRReg1_toyGI,
#'                  disGeneDB = disGeneDB,
#'                  geneIDs = c("gene","name"))                
#' 
#'  
#' @import GenomicRanges
#' @import InteractionSet
#' 
#' @keywords internal
reg2trait <- function(interactions,
                      disGeneDB,
                      geneIDs=c("name2","GeneName"),
                                  ...){
  
  # Finding info about genes when variant-gene-disease
  # association info available
  
  if (class(disGeneDB)=="GInteractions"){
    # expect GRanges or GInteractions
    
    if (class(interactions)=="GInteractions"){
      
    OverlapGI <-  DataFrame(findOverlaps(first(interactions),
                                         first(disGeneDB)))}
    
    if (class(interactions)=="GRanges"){
      
      OverlapGI <-  DataFrame(findOverlaps(interactions,
                                           first(disGeneDB)))
    }
    
    
    # interactions for which overlap in disease-gene-variant 
    # database is identified
    interactions <- interactions[OverlapGI$queryHits]
    # interactions for which overlap in disease-gene-variant 
    # database is identified
    disGeneDB <- disGeneDB[OverlapGI$subjectHits]
    
    }
    
  
  if (class(disGeneDB)=="GRanges"){
    # expect GRanges or GInteractions
    
    NameMatch <- unlist(match(mcols(interactions)[geneIDs[1]],
                              mcols(disGeneDB)[geneIDs[2]]))
    
    if (!all(is.na(NameMatch))) {  # if any match identified
      
      # identify appropriate diseases 
      disGeneDB <- disGeneDB[NameMatch[complete.cases(NameMatch)]]
      
      # identify appropriate interactions objects which have available 
      # info about traits
      IndexForInterObj <- (1:length(interactions))[complete.cases(NameMatch)]
      
      interactions <- interactions[IndexForInterObj]
      mcols(interactions) <- c(mcols(interactions),
                               mcols(disGeneDB) )
    }
    
    
    return(interactions)
    
    
    
  }
     equalGenes <-  which(unlist(mcols(interactions)[geneIDs[1]]==
                                  mcols(disGeneDB)[geneIDs[2]]))
    
    
    mcols(interactions) <- DataFrame(mcols(interactions),
                                        disGeneDB)
    
    interactions <- interactions[equalGenes] 
    # remove X column
    mcols(interactions) <- mcols(interactions)[
      colnames(mcols(interactions))!="X"]
    return(interactions)
    
    
  }

  


#' Compare two objects of genomic region~gene~trait associations
#'   
#' 
#' This function reports per trait number of reported 
#' genomic region~gene~trait from two different objects.
#' 
#' @param  interactionsTraits A GInteractions object, likely produced by 
#' [\code{\link{reg2trait}}]. This object stores info about
#' genomic region~gene~trait associations, as follows: anchor1 corresponds
#' to genomic region (likely SNP location), anchor2 corresponds to gene 
#' location, and gene symbol and associated disease is stored as a meta-data.
#' 
#' @param   interactionsTraits2 A GInteractions object, likely produced by 
#' [\code{\link{reg2trait}}]. This object stores info about
#' genomic region~gene~trait associations, as follows: anchor1 corresponds
#' to genomic region (likely SNP location), anchor2 corresponds to gene 
#' 
#' @param traitID character vector of length 2 which stores information
#' about column names where trait info are stored for both input objects. 
#' By default c("Trait","Trait")
#' 
#' @param naming a character vector (default: c("Traits1,"Trait2").
#' Column names in the output object.
#' 
#' @return a matrix with rows corresponding to traits, and two columns with
#' info about gene counts for these traits. Each column gathers info from one 
#' interactionsTraits GInteractions object (variant~gene interactions + trait 
#' association stored as meta-data).
#'  
#' @details For each trait that is present in at least one of the interactionsTraits 
#' GInteractions objects report how many associated genes is present in the 
#' first GInteractions object, and how many is present in the 2nd 
#' GInteractions object. It allows an easy identification of traits for which
#' [\code{\link{annotateGenomicRegions}}] analysis improved number of disease
#' associated genes.
#' 
#' @author Inga Patarcic
#' 
#' @examples library(GenomicRanges)
#' library(InteractionSet)
#' 
#' # 1. creating interaction object 1 with added diseases 
#' 
#' GRReg1_toyDisease <- GRReg1_toy
#' GRReg1_toyDisease$disease <- paste0("dis",c(1,1,3:(length(GRReg1_toy))))
#' 
#'  
#' # 2. creating interaction object 2 with added diseases 
#' GRReg2_toyDisease <- GRReg2_toy
#' GRReg2_toyDisease$disease <- paste0("dis",1:length(GRReg2_toy))
#' 
#' 
#' # comparing them with compareReg2Traits
#' compareReg2Traits(interactionsTraits=GRReg1_toyDisease,
#'                   interactionsTraits2=GRReg2_toyDisease,
#'                   traitID=c("disease","disease"),
#'                   naming=c("Trait1","Trait2"))
#'  
#' @import GenomicRanges
#' @import InteractionSet
#' 
#' @keywords internal
compareReg2Traits <- function(interactionsTraits,
                              interactionsTraits2,
                              traitID=c("Trait","Trait"),
                              naming=c("Trait1","Trait2")){
  
  # obtain per disease statistics
  reg2t.table <- table(mcols(interactionsTraits)[traitID[1]])
  reg2t2.table <- table(mcols(interactionsTraits2)[traitID[2]])
  
  # select diseases present in both datasets  
  Diseases <- names(reg2t.table)[names(reg2t.table)%in%names(reg2t2.table)]
  
  Compare2datasets <- cbind(reg2t.table[Diseases],reg2t2.table[Diseases])
  
  colnames(Compare2datasets) <- c("interactionsTraits","interactionsTraits2")
  
  # filtering when both entries have diseases with 0 gene counts
  
  Compare2datasets <- Compare2datasets[apply(Compare2datasets,1,
                                             function(x){return(!all(x==0))}),]
  
  colnames(Compare2datasets) <- naming
  
  return(Compare2datasets)
}


