# create example GRangesList from Inga's raw matrices

# read in data 
m1=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000140718_FTO.rds")
m2=readRDS("data/sample_rawActivityMatrices/Roadmap/H3K27ac/ENSG00000177508_IRX3.rds")

m1=apply(m1, 2,as.numeric) # change to numeric matrix
m2=apply(m2, 2,as.numeric) # change to numeric matrix

loc1=do.call("rbind",strsplit(colnames(m1),"_"))
loc2=do.call("rbind",strsplit(colnames(m2),"_"))

#' object in the list represents regulatory activity and gene expression accross
#' the same samples. The GRanges objects have the following 
#' metadata columns:
#' featureType: either "gene" or "regulatory"
#' name: name/id for gene and enhancers. Gene name could be id from a database
#' enhancer name should be in the format as follows "chr:start-end"
#' name2: a secondary name for the feature, such as gene symbol "PAX6" etc. not
#' necessary for enhancers could be NA
#' other columns: numeric values for gene expression or regulatory actvity.
#' Column names represent sample names/ids.
#' 

# conver to GRanges and final objects
g1=GRanges(seqnames=loc1[,1],IRanges(as.numeric(loc1[,2])-1,
        as.numeric(loc1[,3])))
rownames(m1)=paste0("E0",1:nrow(m1))
mcols(g1)=DataFrame(featureType=c("gene",rep("enhancer",99)),
                    name=c("ENSG00000140718",paste(paste(loc1[-1,1],loc1[-1,2],sep=c(":")),loc1[-1,3],sep="-")),
                    name2=c("FTO",rep(NA,99)),
                    t(m1))


g2=GRanges(seqnames=loc2[,1],IRanges(as.numeric(loc2[,2])-1,
                                     as.numeric(loc2[,3])))
rownames(m2)=paste0("E0",1:nrow(m2))
mcols(g2)=DataFrame(featureType=c("gene",rep("enhancer",69)),
                    name=c("ENSG00000177508",
                           paste(paste(loc2[-1,1],loc2[-1,2],sep=c(":")),loc2[-1,3],sep="-")),
                    name2=c("IRX3",rep(NA,69)),
                    t(m2))

output=GRangesList(ENSG00000140718=g1,ENSG00000177508=g2)
saveRDS(output,"data/sampleRegActivityAroundTSS.rds")
output=readRDS("data/sampleRegActivityAroundTSS.rds")

allpairs=unlist(endoapply(output, function(x){
  
  reg=GRanges(x[x$featureType != "gene",]$name)
  gene=x[x$featureType == "gene",c("name" ,"name2")]
  
  grpairs=rep(gene,length(reg))
  grpairs$reg=reg
  
  grpairs
  
}))

allpairs$LiebermanHiC=as.integer(rnorm(length(allpairs),4))
allpairs$GTEx=as.integer(rnorm(length(allpairs),30))
allpairs$LiebermanHiC[allpairs$LiebermanHiC>6]=6
allpairs$LiebermanHiC[allpairs$LiebermanHiC < 0]=0

allpairs$LiebermanHiC[4:5]=0
allpairs$GTEx[4:5]=0

saveRDS(allpairs,"targetPrediction/data/sampleAllPairsAndBencmark.rds")
readRDS("data/sampleAllPairsAndBencmark.rds")
