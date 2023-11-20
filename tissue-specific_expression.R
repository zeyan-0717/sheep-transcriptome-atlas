library(tidyverse)

# set dir
setwd("<SET_YOUR_PATH>")
# input file
load("<SET_YOUR_PATH/_raw_count_matrix.RData>")
info <- readRDS("<SET_YOUR_PATH/all_sample_phenotype.rds>")
head(info)

## step 1: prepare phenotype table
meta <- data.frame(sampleid=rownames(tpm))
head(meta)
head(info)
meta <- inner_join(meta,info, by= "sampleid")
rownames(meta) <- meta$sampleid
head(meta) # contain all sample
meta1 <- distinct(meta[,c("tissueID","tissue","type")])
rownames(meta1) <- meta1$tissue
head(meta1)

## step 2: combat calibration
library(sva)
# filter low expression gene
total_count <- apply(t(raw_count),1,sum) %>% as.data.frame()
names(total_count) <- c("raw_count")
total <- data.frame(gene_id = rownames(total_count), count = total_count$raw_count) %>% filter(count>1) %>% pull(gene_id) 
raw_count <- raw_count[,total] # raw_count > 1
dim(raw_count)
# combat input file
combat_count <- ComBat_seq(counts = t(raw_count),batch = meta$group, group = meta$tissue)
 # combat_tpm <- ComBat_seq(tpm_normal, batch = meta$group, group = NULL)
combat_count[1:10,1:8]
t(raw_count)[1:10,1:8]
save(combat_count,file = "combat_all_count.RData")
rm(raw_count)

detach("package:sva")


## step 3: identify tissue-specific genes
library(limma) # t-test, require data obey normal distribution
# normal transform
combat_count_normal <- log2(combat_count+1)
Tissue <- unique(meta$tissue)
# limma R package
TSGs <- list()
TSG_res <- list()
df3 <- NA
for (i in 1:length(Tissue)) {
  # prepare meta table
  df1 <- filter(meta, tissue == Tissue[i])
  df2 <- filter(meta, (tissue != Tissue[i]) & (type != unique(df1$type)))
  df2$tissue <- c("others")
  df <- rbind(df1,df2)
  design <- model.matrix(~0+factor(df$tissue)) 
  colnames(design) <- levels(factor(df$tissue))
  rownames(design) <- rownames(df)
  # prepare expression matrix
  matrix1 <- combat_count_normal[,rownames(df1)] # expression of target tissue (exclude sub-groups)
  matrix2 <- combat_count_normal[,rownames(df2)] # expression of other tissues
  matrix3 <- cbind(matrix1,matrix2)
  # analysis
  # step 1
  contrast_matrix <- makeContrasts(paste(colnames(design),collapse = "-"),levels = design)
  # step 2
  fit1 <- lmFit(matrix3, design)
  fit2 <- contrasts.fit(fit1,contrast_matrix)
  fit2 <- eBayes(fit2)
  temp_res <- topTable(fit2, sort.by = "t", n = Inf) %>% filter(abs(logFC) > 1 & adj.P.Val < 0.01)
  n <- ceiling(0.05*nrow(matrix3)) # top 5%
  top0.05 <- temp_res[1:n,] 
  # n <- ceiling(0.1*nrow(matrix3)) # top 10%
  # top0.1 <- temp_res[1:n,] 
  res <- data.frame(gene_id=rownames(top0.05),top0.05)
  res$tissue <- Tissue[i]
  TSG_res[[paste(Tissue[i],"_TSG_res",sep="")]] <- res
  TSGs[[paste(Tissue[i],"_specific",sep="")]] <- res$gene_id # input file for fst enrichment (list)
  df1 <- res[,c("geneid","tissue")] # input file for annotation (data.frame) 
  df3 <- rbind(df3,df1) # input file for annotation (data.frame)
  
}
head(df3)
df3 <- df3[-1,]

top_gene <- plyr::ddply(.data = df3,.variables = "tissue",.fun = function(x) x[1,])
save(TSG_res,TSGs,top_gene,file = "count_tissue_specific_gene.RData")
saveRDS(df3,"TSG_df.rds")
rm(fit1,fit2,matrix1,matrix2,matrix3,res,temp_res,top0.05,meta1,df,df1,df2,design)
detach("package:limma")
