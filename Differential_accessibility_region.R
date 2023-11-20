library(DiffBind)
library(parallel)
library(tidyverse)

# setdir
setwd("<SET_YOUR_PATH/diffbind>")

# define variables
filter <- dplyr::filter # define conflict functions
arrange <- dplyr::arrange
mytissue <- c("Adipose","Rumen") # "Heart", "Liver", "Spleen", "Lung", "Hypothalamus", "Duodenum",


# ------------------------------------------------------------------------ #
#                              diffbind                                    #
# ------------------------------------------------------------------------ #
# set threshold
p_value <- 0.05
logFC <- 0.5 # fold: absolute log Fold value in dba.report()
# p_value <- 0.05
# logFC <- 1 # fold: absolute log Fold value in dba.report()

for (i in 1:length(mytissue)) 
{
  # input data info.
  df <- read.csv("diffbind_sample_meta.csv",header = T) %>% filter(Tissue == mytissue[i] & Condition == "F0" )
  head(df)
  dbj <- dba(sampleSheet = df)
  # pca plot
  # pdf("test_sample.pdf")
  # dba.plotPCA(DBA = dbj,attributes=DBA_FACTOR, dotSize = 2.5, vColors = c("darkred","steelblue","gray50"),alpha=0.7) 
  # dev.off()
  
  ## within tissues, between groups
  # get peak count
  dbObj <- dba.count(dbj, bUseSummarizeOverlaps=TRUE) 
  saveRDS(dbObj,paste0(mytissue[i],"_F0_PeakCount.rds"))
  # basic statistic
  # info <- dba.show(dbObj)
  # libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads*info$FRiP))
  # rownames(libsizes) <- info$ID
  # write.table(libsizes,paste(tissue[i],"_library_F0.txt",sep=""), row.names = F,quote = F,sep="\t")
  
  # build contrast
  dbObj1 <-  dba.contrast(dbObj,contrast=c("Factor","B_Hu_sheep","A_Hu_sheep")) # case/control - B/A
  dbObj2 <-  dba.contrast(dbObj,contrast=c("Factor","C_tibetan_sheep","A_Hu_sheep")) # C/A
  dbObj3 <-  dba.contrast(dbObj,contrast=c("Factor","B_Hu_sheep","C_tibetan_sheep")) # B/C
  
  # differential peak analysis
  dbObj1_res <- dba.analyze(dbObj1, method=DBA_ALL_METHODS, bBlacklist = FALSE) 
  dbObj2_res <- dba.analyze(dbObj2, method=DBA_ALL_METHODS, bBlacklist = FALSE)
  dbObj3_res <- dba.analyze(dbObj3, method=DBA_ALL_METHODS, bBlacklist = FALSE)
  # DBA result summary
  # dba.show(dbObj, bContrasts=TRUE)
  # res_summary <- dba.show(dbObj, bContrasts=TRUE) %>% as.data.frame() # summary for result
  # write.table(res_summary,paste(tissue[i],"_within_F0_summary.txt",sep=""),row.names = F, quote = F,sep="\t")
  
  # extract result
  # deseq_res <- dba.report(dbObj,contrast=1,method = DBA_DESEQ2, bUsePval = 0.01, fold = 1) %>% as.data.frame()
  # output result
  res1 <- dba.report(dbObj1_res,th = p_value, bUsePval = T, fold = logFC) %>% as.data.frame() %>% filter(seqnames %in% c(1:26,"X")) # group A VS. group B
  res1 <- res1[order(res1$FDR, res1$Fold, decreasing = c(FALSE, TRUE)),]
  res1 <- mutate(res1,sig = case_when(Fold > logFC ~ "up",Fold < logFC*-1 ~ "down"))
  write.table(res1,paste("./parent/",mytissue[i],"_BA_F0.txt",sep=""),row.names = F, quote = F,sep="\t")
  write.table(res1[,c(1:3,12)],paste("./parent/",mytissue[i],"_BA_F0_tmp.bed",sep=""),sep = "\t",quote = F,row.names = F,col.names = F) # *.bed
  
  res2 <- dba.report(dbObj2_res,th = p_value, bUsePval = T, fold = logFC) %>% as.data.frame() %>% filter(seqnames %in% c(1:26,"X")) # group A VS. group C
  res2 <- res2[order(res2$FDR, res2$Fold, decreasing = c(FALSE, TRUE)),]
  res2 <- mutate(res2,sig = case_when(Fold > logFC ~ "up",Fold < logFC*-1 ~ "down"))
  write.table(res2,paste("./parent/",mytissue[i],"_CA_F0.txt",sep=""),row.names = F, quote = F,sep="\t")
  write.table(res2[,c(1:3,12)],paste("./parent/",mytissue[i],"_CA_F0_tmp.bed",sep=""),sep = "\t",quote = F,row.names = F,col.names = F) # *.bed
  
  res3 <- dba.report(dbObj3_res,th = p_value, bUsePval = T, fold = logFC) %>% as.data.frame() %>% filter(seqnames %in% c(1:26,"X")) # group C VS. group B
  res3 <- res3[order(res3$FDR, res3$Fold, decreasing = c(FALSE, TRUE)),]
  res3 <- mutate(res3,sig = case_when(Fold > logFC ~ "up",Fold < logFC*-1 ~ "down"))
  write.table(res3,paste("./parent/",mytissue[i],"_BC_F0.txt",sep=""),row.names = F, quote = F,sep="\t")
  write.table(res3[,c(1:3,12)],paste("./parent/",mytissue[i],"_BC_F0_tmp.bed",sep=""),sep = "\t",quote = F,row.names = F,col.names = F) # *.bed
  
}
