library(tidyverse)
library(maSigPro)
library(DESeq2)
library(Mfuzz)

# set dir
setwd("<SET_YOUR_PATH>")

# input file
MATRIX <- get(load("<SET_YOUR_PATH/raw_count_matrix.RData>")) 
MATRIX[1:5,1:5]
# load("<SET_YOUR_PATH/tpm_matrix.RData>")
# raw_count[1:5,1:5]

sample_info <- readRDS("<SET_YOUR_PATH/sample_phenotype.rds>")
head(sample_info)
table(sample_info$group)
tissue_info <- readRDS("<SET_YOUR_PATH/tissue_info.rds>")
head(tissue_info)

Tissue <- tissue_info$tissue
Tissue

DCG <- list()
for (i in 1:length(Tissue)) {
  #### step 1: time serise analysis (maSigPro) ####
  ## normalize matrix
  df1 <- filter(sample_info,tissue==Tissue[i] & breed=="Hu")
  df1$Time <- str_replace_all(df1$group,c("D0_Hu"="1","D07_Hu"="7","D14_Hu"="14", "D21_Hu"="21","M08_Hu"="295")) %>% as.numeric()
  df1$group <- as.factor(df1$group)
  
  Mat <- t(MATRIX[df1$sampleid,])
  dds <- DESeqDataSetFromMatrix(countData = Mat, colData = df1[,c("sampleid","group")], design = ~ group)
  dds <- dds[rowSums(counts(dds)) > 10,]
  dim(dds)
  Mat_normal <- vst(dds,blind = F) %>% assay()

  ## prepare meta sample info
  meta <- data.frame(Time=df1$Time, Replicate=df1$Time,Group=1)
  head(meta)
  rownames(meta) <- df1$sampleid
  head(meta)
  table(meta$Time)
  
  ## build model
  design <- make.design.matrix(meta, degree = length(unique(meta$Time)) - 1)
  ## p.vector
  fit <- p.vector(Mat_normal, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
  # fit$i # significant genes
  # head(fit$SELEC) # expression matrix with sig. genes
  ## T.fit
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
  ## lists of sig. genes
  sigs <- get.siggenes(tstep, rsq = 0.3, vars = "all")
  # names(sigs)
  dim(sigs$sig.genes$sig.profiles) # check gene number
  DCG[[paste(Tissue[i],"_DCG",sep="")]] <- sigs$sig.genes$sig.profiles

}
rm(raw_count)


for (j in 1:length(Tissue[j])) {
  #### step 2: c-means clustering ####
  # calculate mean value of expression
  df1 <- filter(sample_info,tissue==Tissue[j] & breed=="Hu")
  head(df1)
  scg_mat <- DCG[[paste(Tissue[j],"_DCG",sep="")]]
  # scg_mat[1:5,1:5]
  D0 <- filter(df1,group == "D0_Hu") %>% pull(sampleid)
  D0_avg <- scg_mat[,D0] %>% apply(1,FUN = mean)
  D7 <- filter(df1,group == "D07_Hu") %>% pull(sampleid)
  D7_avg <- scg_mat[,D7] %>% apply(1,FUN = mean)
  D14 <- filter(df1,group == "D14_Hu") %>% pull(sampleid)
  D14_avg <- scg_mat[,D14] %>% apply(1,FUN = mean)
  D21 <- filter(df1,group == "D21_Hu") %>% pull(sampleid)
  D21_avg <- scg_mat[,D21] %>% apply(1,FUN = mean)
  M8 <- filter(df1,group == "M08_Hu") %>% pull(sampleid)
  M8_avg <- scg_mat[,M8] %>% apply(1,FUN = mean)
  avg_exp <- data.frame(D0 = D0_avg, D7=D7_avg,D14=D14_avg,D21=D21_avg,M08=M8_avg) %>% as.matrix()
  head(avg_exp)
  rm(D0,D0_avg,D6,D6_avg,D13,D13_avg,D20,D20_avg,M8_avg,M8)
  
  ## analysis
  mfuzz_class <- new("ExpressionSet",exprs=avg_exp)
  mfuzz_class <- standardise(mfuzz_class)
  
  set.seed(4321)
  cluster_num <- 4
  mfuzz_cluster <- mfuzz(mfuzz_class, c=cluster_num, m=mestimate(mfuzz_class))
  mfuzz.plot2(mfuzz_class, mfuzz_cluster, centre.col = "black",
              mfrow=c(3,3), x11= FALSE,centre=TRUE, min.mem = 0.3,
              time.labels=colnames(avg_exp),colo = rev(mycolor1))
  ## extract clustering results
  names(mfuzz_cluster$size) <- 1:cluster_num
  mfuzz_cluster$size # check the number of each cluster
  
  ## visualize
  # mycolor <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
  mycolor1 <- colorspace::sequential_hcl(5000,"BuPu")
  pdf(paste(Tissue[j],"_c_means.pdf",sep=""),width = 9,height = 9)
  mfuzz.plot2(mfuzz_class, mfuzz_cluster, centre.col = "black",
              mfrow=c(3,3), x11= FALSE,centre=TRUE, min.mem = 0.3,
              time.labels=colnames(avg_exp),colo = rev(mycolor1))
  dev.off()

  # output original expression matrix (VST normalized)
  gene_cluster <- data.frame(avg_exp[names(mfuzz_cluster$cluster),], 
                             Cluster= str_c("Cluster","_",mfuzz_cluster$cluster))
  # output  standardlized matrix (mFuzz normalized)
  gene_standard <- mfuzz_class@assayData$exprs
  gene_cluster_standard <- data.frame(gene_standard[names(mfuzz_cluster$cluster),],
                                      Cluster= str_c("Cluster","_",mfuzz_cluster$cluster))
  
  write.table(gene_cluster,paste(Tissue[j],"_c_means_vst.txt",sep=""),row.names = F,quote = F,sep="\t")
  write.table(gene_cluster_standard,paste(Tissue[j],"_c_means_standard.txt",sep=""),row.names = F,quote = F,sep="\t")
  
}


  
