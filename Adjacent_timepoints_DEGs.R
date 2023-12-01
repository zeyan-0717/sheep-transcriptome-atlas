library(tidyverse)
library(DESeq2)
options(stringsAsFactors = T)

## set dir
setwd("<SET_YOUR_PATH/DEG/raw_DEG_res>")
## input files
# load("<SET_YOUR_PATH/raw_count_matrix.RData>")
# load("<SET_YOUR_PATH/tpm_matrix.RData>")
raw_count <- readRDS("F:/sheep_adaptation/data/mature_raw_count_matrix.rds")
# check data
raw_count[1:5,1:5] 
dim(raw_count)

sample_info <- readRDS("<SET_YOUR_PATH/sample_phenotype.rds>")
dim(sample_info)
head(sample_info)
table(sample_info$group)


## define group
group1 <- c("D0_Hu","D06_Hu") 
group2 <- c("D06_Hu","D13_Hu")
group3 <- c("D13_Hu","D20_Hu")
group4 <- c("D20_Hu","M08_Hu")
group5 <- c("M08_Hu","Tibetan")
Group <- list(group1,group2,group3,group4,group5)

# Adjacent time-points
pairwise1 <- list(compare1 = c("group","D07_Hu","D0_Hu"),
                  compare2 = c("group","D14_Hu","D07_Hu"),
                  compare3 = c("group","D21_Hu","D14_Hu"),
                  compare4 = c("group","M08_Hu","D21_Hu"),
                  compare5 = c("group","Tibetan","M08_Hu")) # case VS. control

Tissue <- unique(sample_info$tissue)
Tissue

TB <- NA
for (i in 1:length(Tissue)) {
  meta <- dplyr::filter(sample_info,tissue==Tissue[i])
  df <- NA
  DF <- data.frame(time=NA,up =NA, down=NA)
  
  for (j in 1:length(Group)) {
    meta1 <- dplyr::filter(meta,group %in% unlist(Group[[j]]))
    raw_count <- as.matrix(raw_count)
    count_mat <- t(raw_count[rownames(meta1),])
    dim(count_mat)
    meta1$group <- factor(meta1$group,levels = unlist(Group[[j]][2:3]))
    
    ## check meta_tb & count matrix
    countMatrix <- as.matrix(count_mat)
    all(colnames(countMatrix) %in% rownames(meta1))
    all(colnames(countMatrix) == rownames(meta1))
    rm(count_mat)
    ## normalize matrix
    dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = meta1[,c("sampleid","group")], design = ~ group)
    dds <- dds[rowSums(counts(dds)) > 4,] # filter sum(counts) < 1
    dds

    #### analysis  ####
    ## step 1 :
    p_value <- <SET_YOUR_PARAMETER>
    FC <- <SET_YOUR_PARAMETER>
    
    dds_wd <- DESeq(dds) # DEG analysis
    res_wd <- results(dds_wd) # build the results table
    resultsNames(dds_wd) # see which object to compare
    # rld_mat_wd <- assay(vsd) # convert results to data frame
    res1 <- results(dds_wd,contrast = pairwise1[[j]])
    summary(res1)
    table(res1$padj < p_value & abs(res1$log2FoldChange) >= FC)
    
    # step 2: select DEG
    df1 <- data.frame(res1, stringsAsFactors = FALSE, check.names = FALSE) 
    df1 <- df1[order(df1$padj, df1$log2FoldChange, decreasing = c(FALSE,TRUE)),] 
    df1[which(df1$log2FoldChange > FC & df1$padj < p_value),'change'] <- 'up'
    df1[which(df1$log2FoldChange < -FC & df1$padj < p_value),'change'] <- 'down'
    df1[which(abs(df1$log2FoldChange) <= FC | df1$padj >= p_value),'change'] <- 'none'
    head(df1)
    table(df1$change)
    
    df2 <- filter(df1, change %in% c('up', 'down'))
    df2 <- data.frame(gene_id=rownames(df2),df2)
    table(df2$change)
    df2$group <- paste("adjacent_",j,sep="")
    df <- rbind(df,df2)
    
    # summarize up and dw gene number
    DF[j,c("down")] <- table(df2$change)[1]
    DF[j,c("up")] <- table(df2$change)[2]
    DF[j,c("time")] <- paste("adjacent_",j,sep="")
    DF$up <- str_replace_na(DF$up, replacement = 0) %>% as.numeric()
    DF$down <- str_replace_na(DF$down, replacement = 0) %>% as.numeric()
    DF$total <- apply(DF[,c("up","down")],1,sum)
    DF$tissue <- Tissue[i]
  }

  TB <- rbind(TB,DF)
  dim(df)
  df <- df[-1,]
  write.table(df,paste(Tissue[i],"_adjacent.txt",sep=""),row.names = F,quote = F,sep="\t")
}

head(TB)
TB <- TB[-1,]
dim(TB)
table(TB$tissue)
write.table(TB,"DEG_stat_all.txt",row.names = F,quote = F,sep="\t")



#### visualize change tendency ####
df <- read.table("<SET_YOUR_PATH/DEG/DEG_stat_all.txt>",header = T)
table(df$tissue)
length(unique(df$tissue))

df <- filter(TB[,-4],time != "adjacent_5")
df <- gather(df[,-4], key = "change",value = "number",2:3)
df$change <- factor(df$change, ordered = T, unique(df$change))
head(df)

color1 <- c("Pituitary"="darkorange1","Hypothalamus"="#ef8737","Cerebellum"="#ffb242","Cerebrum"="#ffd353",  
            "Abomasum"="#275024", "Duodenum"="#47632a", "Jejunum"="#748f46", "Ileum"="#98ab76","Colon"="#ced1af",
            "Leukocyte"="#dc322f", "Spleen"="#f05b43","Heart" = "#d33682", "Muscle"="#9f2d55","Artery" = "#CD5B45",
            "Adipose"="darkcyan","Lung"="#FDAF9199", "Liver"="#6c71c4", "Kidney"="#016392","Rumen" = "#B4EEB4")

pdf("<YOUR_PATH>/DEG_change.pdf",width = 12, height = 16)
ggplot(df,aes(time,number,linetype=change,color=tissue)) +
  geom_line(aes(time,number,group=change),size=0.8)+
  geom_point(size=2) +
  labs(x=NULL, y=" Number of DEGs") +
  scale_color_manual(values = color1) +
  scale_x_discrete(labels = c("0d VS. 6d","6d VS. 13d","13d VS.20d","20d VS. 8mon")) +
  theme_bw() +
  theme(axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12,face = "bold",angle = 30,vjust = 1,hjust = 1),
        axis.text.y=element_text(size = 12,face = "bold"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size=13,face = "bold"),
        panel.grid=element_blank()) +
  facet_wrap(~tissue,ncol=4, scales = "free_y")
dev.off()

rm(df1,df2,df3,dds,dds_wd,countMatrix,res_wd,res1,raw_count,pairwise1,meta,meta1,TB,df,DF);gc()

