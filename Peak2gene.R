library(tidyverse)
library(rtracklayer) 
library(InTAD) 

setwd("<SET_YOUR_PATH>")

# ------------------------------------------------------------------------ #
#                       input files (as GRange format)                     #   
# ------------------------------------------------------------------------ #
## input peak coordinate
peak_GR <- import.bed("./data/all_peak_merged.bed")
head(peak_GR)
head(as.data.frame(peak_GR))

## input TAD file (*.bed)
tad_GR <- import.bed("./data/tad_40k.bed")
head(tad_GR)
head(as.data.frame(tad_GR))

## input gene coordinate
extraCols_gene_coordinate <- c(gene_id=c("character"),gene_name=("character"))
import_bed6 <- function(x) {
  import(x, format="BED", extraCols = extraCols_gene_coordinate)
}
gene_GR <- import_bed6("./data/inTAD_gene_coordinate.bed")
head(gene_GR)
head(as.data.frame(gene_GR))


## input peak and gene expresseion  files
# expression files
tpm <- get(load("./data/all_tpm_matrix.RData"))
tpm[1:5,1:5]
gene_sample <- readRDS("./data/all_sample_phenotype.rds")
dim(gene_sample)
head(gene_sample)
# peak files
rpkm <- readRDS("./data/peakRPKM.rds")
rpkm[1:5,1:5]
rpkm_mat <- rpkm[,4:ncol(rpkm)] %>% as.matrix()
rownames(rpkm_mat) <- str_c(rpkm$chr,":",rpkm$start,"-",rpkm$end)
rpkm_mat[1:5,1:5]
# peak sample
peak_sample <- readRDS("./data/peak_sample_info.rds")
head(peak_sample)
rm(rpkm)

# check gene_coordinate whether match to expression matrix
gene_coordinate <- as.data.frame(gene_GR)
head(gene_coordinate)
dim(gene_coordinate) 
dim(tpm) 




# ------------------------------------------------------------------------ #
#                              link peak to gene                           #
# ------------------------------------------------------------------------ #
mytissue1 <- c("Liver", "Spleen", "Duodenum", "Adipose", "Rumen")  
Group1 <- c("D0_Hu","M08_Hu","Tibetan")
mytissue2 <- c("Heart","Lung","Hypothalamus")
Group2 <- c("D0_Hu","M08_Hu","Tibetan","Hu_plain_lamb","Hu_plateau_lamb","Tibetan_lamb") 
# set above parameters base on your data

#### for group1 and tissue1 ####
list1 <- list()
NAME1 <- NA
for (i in 1:length(mytissue1)) {
  ## expression file
  df2 <- filter(gene_sample,tissue == mytissue1[i], generation=="F0")
  exp_sample <- filter(df2,group %in% Group1) %>% pull(sampleid)
  gene_mat1 <- t(tpm[exp_sample, gene_coordinate$gene_id]) # extract corresponding sample and gene
  # compute median and mean tpm for gene expression
  list2 <- list()
  for (j in 1:length(Group1)) {
    df3 <- filter(df2, group == Group1[j]) %>% pull(sampleid)
    gene_median <- apply(gene_mat1, 1, median)
    gene_mean <- apply(gene_mat1, 1,mean)
    list2[[Group1[j]]] <- data.frame(gene_median,gene_mean) %>% as.matrix()
  }
  list1[[mytissue1[i]]] <- do.call(cbind,list2)
  NAMES1 <- c(paste("A_",mytissue1[i],"_1",sep=""),
              paste("A_",mytissue1[i],"_2",sep=""),
              paste("B_",mytissue1[i],"_1",sep=""),
              paste("B_",mytissue1[i],"_2",sep=""),
              paste("C_",mytissue1[i],"_1",sep=""),
              paste("C_",mytissue1[i],"_2",sep=""))
  NAME1 <- c(NAME1,NAMES1)
}
rm(gene_mat1)

gene_mat1 <- do.call(cbind,list1)
dim(gene_mat1)
gene_mat1[1:5,1:5]
rm(df2,df3,list1,list2)
colnames(gene_mat1) <- NAME1[2:length(NAME1)]

# extract corresponding samples
peak_mat1 <- rpkm_mat[,colnames(gene_mat1)]
peak_mat1 <- log2(peak_mat1+1)
peak_mat1[1:5,1:5]


#### for group2 and tissue2 ####
list1 <- list()
NAME2 <- NA
for (i in 1:length(mytissue2)) {
  # expression file
  DF2 <- filter(gene_sample,tissue == mytissue2[i])
  exp_sample <- filter(DF2,group %in% Group2) %>% pull(sampleid)
  GENE_MAT1 <- t(tpm[exp_sample, gene_coordinate$gene_id])
  # compute median and mean tpm for gene expression
  list2 <- list()
  for (j in 1:length(Group2)) {
    DF3 <- filter(DF2, group == Group2[j]) %>% pull(sampleid)
    gene_median <- apply(GENE_MAT1, 1, median)
    gene_mean <- apply(GENE_MAT1, 1,mean)
    list2[[Group2[j]]] <- data.frame(gene_median,gene_mean) %>% as.matrix()
  }
  list1[[mytissue2[i]]] <- do.call(cbind,list2)
  NAMES2 <- c(paste("A_",mytissue2[i],"_1",sep=""),
             paste("A_",mytissue2[i],"_2",sep=""),
             paste("B_",mytissue2[i],"_1",sep=""),
             paste("B_",mytissue2[i],"_2",sep=""),
             paste("C_",mytissue2[i],"_1",sep=""),
             paste("C_",mytissue2[i],"_2",sep=""),
             paste("a_",mytissue2[i],"_1",sep=""),
             paste("a_",mytissue2[i],"_2",sep=""),
             paste("b_",mytissue2[i],"_1",sep=""),
             paste("b_",mytissue2[i],"_2",sep=""),
             paste("c_",mytissue2[i],"_1",sep=""),
             paste("c_",mytissue2[i],"_2",sep=""))
  NAME2 <- c(NAME2,NAMES2)
}
rm(GENE_MAT1)

GENE_MAT1 <- do.call(cbind,list1)
dim(GENE_MAT1)
GENE_MAT1[1:5,1:5]
rm(DF2,DF3)
colnames(GENE_MAT1) <- NAME2[2:length(NAME2)]

# extract corresponding samples
peak_mat2 <- rpkm_mat[,colnames(GENE_MAT1)]
peak_mat2 <- log2(peak_mat2+1)
peak_mat2[1:5,1:5]
gene_mat <- cbind(gene_mat1,GENE_MAT1)
peak_mat <- cbind(peak_mat1,peak_mat2)
# colnames(gene_mat) <- c(df2$SampleID,DF2$SampleID)


# - ----------------------------------------------------------------------- #
#                               analysis -- InTAD                           #
# - ----------------------------------------------------------------------- #
# check data
summary(colnames(gene_mat) == colnames(peak_mat)) # The names of ATAC-seq samples should match with gene expression datasets
length(peak_GR) == nrow(peak_mat) # compare number of signal regions and in the input table
# matched sample info.
df2 <- filter(peak_sample, Generation == "F0" ,Tissue %in% mytissue1)
DF2 <- filter(peak_sample,Tissue %in% mytissue2)

link_sample <- rbind(df2,DF2)
rownames(link_sample) <- link_sample$SampleID
# Created signals and genes object for all samples
inTadSig <- newSigInTAD(peak_mat,peak_GR,gene_mat,gene_GR,link_sample) # Default: log2(TPM+1)
inTadSig
inTadSig_filter <- filterGeneExpr(inTadSig,cutVal=1) # default:cutVal=0
rm(inTadSig)
# combine signals and genes in TADs
combine_inTadSig <- combineInTAD(inTadSig_filter,tad_GR)
# compute correlation
corData <- findCorrelation(combine_inTadSig,adj.pval = TRUE) # default: pearson cor
# filter result
top_cor_0.05 <- filter(corData, cor>=0.25, qvalue<=0.05) # corRank == 1,
top_cor_0.1 <- filter(corData, cor>=0.25, qvalue<=0.1) # corRank == 1,

# ## output
# # output raw result
saveRDS(corData, "./result/all_link_FDR_res.rds")
write.table(corData,"./result/all_link_FDR_res.txt",row.names = F,quote = F,sep="\t")
write.table(top_cor_0.05[,-4], "./result/all_FDR0.05_linkGene.txt",row.names = F,quote = F,sep="\t")
write.table(top_cor_0.1[,-4], "./result/all_FDR0.1_linkGene.txt",row.names = F,quote = F,sep="\t")
# output bed file
df2 <- str_split(top_cor_0.05$peakid,":",simplify = TRUE) %>% as.data.frame()
df3 <- str_split(df2$V2,"-",simplify = TRUE) %>% as.data.frame()
df4 <- data.frame(df2$V1,df3,top_cor_0.05$gene)
head(df4)
write.table(df4,"./result/all_0.05_FDR_linkGene.bed",row.names = F,col.names = F,quote = F,sep="\t")
rm(df2,df3,df4)

df2 <- str_split(top_cor_0.1$peakid,":",simplify = TRUE) %>% as.data.frame()
df3 <- str_split(df2$V2,"-",simplify = TRUE) %>% as.data.frame()
df4 <- data.frame(df2$V1,df3,top_cor_0.1$gene)
head(df4)
write.table(df4,"./result/all_0.1_FDR_linkGene.bed",row.names = F,col.names = F,quote = F,sep="\t")

