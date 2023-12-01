library(tidyverse)
options(stringsAsFactors = F)

setwd("<SET_YOUR_PATH>")

##-------------------------------------------------------------##
##                        package prepare                      ##
##-------------------------------------------------------------##
#### ChIPseeker ####
library(ChIPseeker)
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)

#### clusterProfiler ####
library(clusterProfiler)
library(AnnotationHub)
## create background gene
require(AnnotationHub)
hub <- AnnotationHub()
query(hub,"Ovis aries") # determine your species
# AH94094 | org.Ovis_aries.eg.sqlite 
sheep_orgdb <- hub[["AH101101"]]

# create txDb.org
library(GenomicFeatures)
gtf <- c("<SET_YOUR_PATH/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.101.gtf>")
species <- c("Ovis aries")
txdb <- makeTxDbFromGFF(file = gtf, format = "gtf",organism = species)
detach("package:GenomicFeatures")

# convert geneID
library(biomaRt)
library(org.Hs.eg.db)
# load database
ensembl <- useEnsembl(biomart = "genes",dataset ="oarambouillet_gene_ensembl",mirror = "useast") # mirror= useast, uswest, asia
# orthology human gene to GO enrichment
human_go <- function(x1) {
  return(enrichGO(x1$ENTREZID, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, maxGSSize = 800,minGSSize = 15,
                  ont = "BP",pvalueCutoff = 0.05,  pAdjustMethod = "BH",readable = TRUE))
}
# sheep gene to GO enrichment
sheep_go <- function(x2) {
  return(enrichGO(x2, keyType = "ENTREZID", OrgDb = sheep_orgdb, maxGSSize = 800,minGSSize = 15,
                  ont = "BP",pvalueCutoff = 0.05,  pAdjustMethod = "BH",readable = TRUE))
}
go_simplify <- function(y) {
  return(simplify(y, cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang"))
}
GENE_CONVERT <- function(x){
  return(biomaRt::getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name"),
                        filters = "ensembl_gene_id",values = x ,mart = ensembl))
}


##-------------------------------------------------------------##
##                           analysis                          ##
##-------------------------------------------------------------##
#### use peaks in one tissues ####
ts_peak <- readPeakFile("ts_peak_1tissue.bed") # bed file (without header)
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
peakAnno <- annotatePeak(ts_peak,TxDb=txdb, tssRegion=c(-2000,2000)) %>% as.GRanges() %>% as.data.frame()
head(peakAnno)
str(peakAnno)
# saveRDS(peakAnno,"ts_peak_oneTissue.rds")
peakAnno$category <- str_split(string = peakAnno$annotation,pattern = "\\(",simplify = T)[,1] 
peakAnno <- rename(peakAnno,"tissue" = "V4")
length(unique(peakAnno$geneId))
table(peakAnno$category)
saveRDS(peakAnno,"ts_peak1_chipseek_res.rds")

# overlap gene
Tissue <- unique(peakAnno$tissue)
common_gene1 <- list()
for (i in 1:length(Tissue)) {
  gene1 <- TSGs[[paste(Tissue[i],"_specific",sep="")]]
  gene2 <- filter(peakAnno,tissue == Tissue[i]) %>% pull(geneId)
  common_gene1[[paste(Tissue[i])]] <- intersect(gene1,gene2)
  # VennDiagram::venn.diagram(filename =  paste(mytissue[i],"_overlap.tiff",sep=""),
  # list("TSG"=gene1,"TSPG"=gene2), fill=)
}
saveRDS(common_gene1,"tsg_tspg1_ovlap_gene.rds")

# GO analysis
result_BM <- lapply(common_gene, GENE_CONVERT)
sheep_gene <- lapply(result_BM, function(x){unique(na.omit(x$entrezgene_id))}) # remove missing entrezgene_id
orthology_gene <- lapply(result_BM, function(x){unique(na.omit(x$external_gene_name))})
human_gene <- plyr::ldply(orthology_gene, function(x){bitr(x,fromType = "SYMBOL",toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)})
human_golist1 <- plyr::dlply(.data =  human_gene, .variables = ".id", .fun = human_go)
human_golist2 <- lapply(human_golist1, go_simplify) 

go1 <- human_golist2$Adipose@result
go2 <- human_golist2$Duodenum@result
go3 <- human_golist2$Heart@result
go4 <- human_golist2$Hypothalamus@result
go5 <- human_golist2$Liver@result
go6 <- human_golist2$Lung@result
go7 <- human_golist2$Rumen@result
go8 <- human_golist2$Spleen@result


