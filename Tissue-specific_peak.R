library(tidyverse)
options(stringsAsFactors = F)


# set dir
setwd("<SET_YOUR_DIR>")
# input peak sample info
peak_info <- readRDS("<SET_YOUR_DIR>/atac/peak_sample_info.rds")
str(peak_info)

 # input consenus peak count matrix
df <- data.table::fread("<SET_YOUR_DIR>/readRPKM_atac_sorted.tab", sep = "\t",stringsAsFactors = F,header = T)
df[1:5,1:5] # check data
names(df) <- c("chr","start","end",peak_info$SampleID)
df$peakID <- str_c(df$chr,"_",df$start,"_",df$end)
dim(df)

# get RPKM matrix
peak_mat <- df[,4:69] %>% as.matrix()
class(peak_mat)
rownames(peak_mat) <- df$peakID
peak_mat[1:5,1:5]

# normalize
peak_mat_normal <- preprocessCore::normalize.quantiles(peak_mat)
row.names(peak_mat_normal) <- df$peakID
colnames(peak_mat_normal) <- colnames(peak_mat)
peak_mat_normal[1:5,1:5]

##-------------------------------------------------------------##
##              Normalize raw readscount by RPM                ##
##-------------------------------------------------------------##
# data=read.table(filecount,header = T)
# row.names(data)=data[,1]
# data=data[,-1]
# colsum=colSums(data)
# RPM<-function(x){return(x/(colsum/1000000))}
# data_RPM=as.data.frame(t(apply(data,1,RPM)))


##-------------------------------------------------------------##
##                       t-sne clustering                      ##
##-------------------------------------------------------------##
# calculation
set.seed(12324)
tsne_res <- Rtsne(t(peak_mat_normal), dims = 2, pca = T,
                  max_iter = 3000,theta = 0.3,perplexity = 20,verbose = F) # row:sampleID, col:peakID
head(tsne_res)
# extract t-SNE result
tsne <- as.data.frame(tsne_res$Y)
colnames(tsne) <- c("tSNE1","tSNE2")
tsne$tissue <- peak_info$Tissue 
tsne$group <- peak_info$Group
class(tsne)
saveRDS(tsne,"<SET_YOUR_DIR>/atac/atac_tsne.rds")

# visualization
mycolor <- c("Hypothalamus"="#ef8737","Duodenum"="#47632a", "Spleen"="#f05b43","Heart" = "#d33682",
             "Adipose"="darkcyan","Lung"="#FDAF9199", "Liver"="#6c71c4","Rumen" = "#B4EEB4")
p1 <- ggplot(tsne,aes(tSNE1,tSNE2,color=tissue,shape=group)) + 
  geom_point(size=2.8,alpha=0.8) + 
  labs(x="t-SNE1", y="t-SNE2") +
  scale_color_manual(values = mycolor) +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        legend.title = element_text(size=13, face="bold"), 
        legend.text = element_text(size=13, face="bold"),  
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13,face = "bold"),
        axis.text.y=element_text(size = 13,face = "bold"),
        panel.grid=element_blank(),
        plot.margin = margin(t=2,r=2,b=2,l=2,unit = "cm"))
p1

# save as pdf
pdf("F:/sheep_adaptation/atac/atac_tsne.pdf",height = 5 ,width = 7.2)
p1
dev.off()
rm(tsne,tsne_res)


##-------------------------------------------------------------##
##               Identify tissue specific peaks                ##
##-------------------------------------------------------------##
# calculate median/average RPKM 
tissue <- unique(peak_info$Tissue)
tissue
MEDIAN_mat <- NA
for (i in 1:length(tissue)) {
  df1 <- filter(peak_info,Tissue==tissue[i]) %>% pull(SampleID)
  mat1 <- peak_mat_normal[,df1] # use normalized peak set
  MEDIAN_peak <- apply(mat1, 1,median)
  MEDIAN_mat <- cbind(MEDIAN_mat,MEDIAN_peak)
  
}
MEDIAN_mat[1:5,1:5]
MEDIAN_mat <- MEDIAN_mat[,-1]
colnames(MEDIAN_mat) <- tissue
MEDIAN_mat[1:5,1:ncol(MEDIAN_mat)]
MEDIAN_mat_bed <- cbind(df[,1:3],MEDIAN_mat) %>% as.data.frame()
head(MEDIAN_mat_bed)
str(MEDIAN_mat_bed)
write.table(MEDIAN_mat_bed,"../ts_normalized_median_RPKM.bed",row.names = F, quote = F,col.names = F,sep = "\t")
rm(mat1)
saveRDS(MEDIAN_mat,"../ts_normalized_median_RPKM.rds")

## calculate shannon entropy
MEDIAN_mat <- readRDS("<SET_YOUR_DIR>/ts_normalized_median_RPKM.rds")
ratio <- as.data.frame(t(apply(MEDIAN_mat,1,function(x) x/sum(x)))) # calculate Ri
# define function
entropy <- function(a){
  b <- -sum(a*log2(a))
  return(b)
}
entro <- as.data.frame(apply(ratio,1,entropy))
head(entro)
# plot curve pics
sum(is.na(entro$`apply(ratio, 1, entropy)`)) # no missing value
# omit_entro <- as.data.frame(entro[-which(entro[,1]=="NaN"),])
# colnames(omit_entro) <- c("entropy")
sta <- as.data.frame(seq(from=0,to=3.2,by=0.01))
# summarise distribution based on entropy score
for(i in 1:nrow(sta)){
  sub_entropy <- subset(entro, entro$`apply(ratio, 1, entropy)`< sta[i,1])
  sta[i,2] <- nrow(sub_entropy)
}
colnames(sta)=c("entropy","peak_num")
pdf("shannon_entropy.pdf",width = 3.4,height = 3)
ggplot(data = sta, mapping = aes(x = entropy, y = peak_num)) + 
  geom_point(size=0.3)+
  geom_vline(xintercept = 2.5, size=1,color = "tomato3",linetype="dashed") +
  theme_bw()+
  xlab("Entropy score")+
  ylab("Cumulative number of peak") +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  theme(panel.grid=element_blank(),
        plot.margin = margin(t = 1,r = 1,b = 1,l = 1,unit = "cm")) 
dev.off()

## select tissue-specific peak
plot <- cbind(MEDIAN_mat,entro)
sub_plot<- subset(plot,plot$`apply(ratio, 1, entropy)`<2.5) # Peaks with entropy less than 2.5 were selected as tissue-restricted peaks
order_sub_plot <- sub_plot[order(sub_plot$`apply(ratio, 1, entropy)`),]
saveRDS(order_sub_plot,"order_sub_plot.rds")

order_sub_plot <- readRDS("<SET_YOUR_DIR>/order_sub_plot.rds")
tissue_name=colnames(sub_plot[,-9])
for(i in seq(nrow(order_sub_plot))){
  a <- as.matrix(order_sub_plot[i,-9])
  order_sub_plot[i,10] <- tissue_name[which(a==a[which.max(a)],arr.ind=T)[2]]
}
re_order_sub_plot <- order_sub_plot[order(order_sub_plot$V10),]
table(re_order_sub_plot$V10)
re_order_sub_plot <- re_order_sub_plot[,c("Adipose", "Duodenum", "Heart", "Hypothalamus", "Liver", "Lung","Rumen", "Spleen")] 
                               
## visualize 
mycolor2 <- list(tissue=c("Hypothalamus"="#ef8737","Duodenum"="#47632a", "Spleen"="#f05b43",
                          "Heart" = "#d33682", "Adipose"="darkcyan","Lung"="#FDAF9199", 
                          "Liver"="#6c71c4", "Rumen" = "#B4EEB4")) 
# peakID x tissue
pdf("tissue_specific_peak_heatmap.pdf",width = 4.5,height = 4.5)
Heatmap(log2(re_order_sub_plot+1), name="Normailzed peak", 
              viridis::inferno(500), 
        cluster_rows = F,cluster_columns = F,show_row_names = FALSE) 
        # heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(3, "cm"), title_position = "lefttop"))
        # top_annotation = HeatmapAnnotation(tissue = colnames(re_order_sub_plot),col = mycolor2))
# p1 <- draw(ht1, heatmap_legend_side = "bottom")
dev.off()

# tissue x peakID 
pdf("tissue_specific_peak_heatmap_v2.pdf",width = 4.8,height = 4)
ht <- Heatmap(log2(t(re_order_sub_plot+1)), name="Normailzed peak", 
              col= viridis::inferno(500),
              # col = colorspace::sequential_hcl(100,"Purples 3"),# viridis::inferno(100), 
              cluster_rows = F,cluster_columns = F,show_column_names = FALSE,
              heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(3, "cm"), title_position = "lefttop"))
p2 <- draw(ht, heatmap_legend_side = "bottom")
dev.off()

## extract tissue specific peaks
forsub <- cbind(ratio,entro)
sub <- subset(forsub,forsub$`apply(ratio, 1, entropy)`<2.5)
tissue_name=colnames(sub[,-9])
for(i in seq(nrow(order_sub_plot))){
  a=as.matrix(sub[i,-9])
  b=i
  for(i in 1:8){
    if(sort(a,decreasing = T)[i] > 0.1){sub[b,9+i] = tissue_name[which(a==sort(a,decreasing = T)[i])] }
    else{sub[b,9+i] = "FALSE"}
  }
}
nth_1=subset(sub[10],sub[10] != "FALSE")#subset peaks in one tissue
nth_2=subset(sub[11],sub[11] != "FALSE")#subset peaks in two tissues
nth_3=subset(sub[12],sub[12] != "FALSE")#subset peaks in three tissues
nth_4=subset(sub[13],sub[13] != "FALSE")#subset peaks in four tissues
nth_5=subset(sub[14],sub[14] != "FALSE")#subset peaks in five tissues
nth_6=subset(sub[15],sub[15] != "FALSE")#subset peaks in six tissues
colnames(nth_1)=c("tissue")
colnames(nth_2)=c("tissue")
colnames(nth_3)=c("tissue")
colnames(nth_4)=c("tissue")
colnames(nth_5)=c("tissue")
colnames(nth_6)=c("tissue")
nth_1$id=row.names(nth_1)
nth_2$id=row.names(nth_2)
nth_3$id=row.names(nth_3)
nth_4$id=row.names(nth_4)
nth_5$id=row.names(nth_5)
nth_6$id=row.names(nth_6)

ALL_peaks=rbind(nth_1,nth_2)
head(ALL_peaks)
names(ALL_peaks) <- c("tissue","peakID")
table(ALL_peaks$tissue)
## subset peaks in two tissues
for(i in 1:8){
  test <- subset(ALL_peaks,ALL_peaks$tissue==tissue_name[i])[2]
  df1 <- str_split(test$peakID,"_",simplify = T)
  write.table(df1, paste(tissue_name[i],"_ts_peaks.bed",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
}

## subset peaks in one tissues (may cause unsufficient peak)
table(nth_1$tissue)
names(nth_1) <- c("tissue","peakID")
nth_1_peak <- left_join(nth_1,df[,c("chr","start","end","peakID")] ,by="peakID")
head(nth_1_peak)
str(nth_1_peak)
write.table(nth_1_peak[,c("chr","start","end","tissue")],"ts_peak_1tissue.bed",
            sep="\t",row.names = F,quote = F,col.names = F)


##-------------------------------------------------------------##
##                    identify conservative peak               ##
##-------------------------------------------------------------##
##
plot <- cbind(MEDIAN_mat,entro)
sub_plot<- arrange(plot, desc(plot$`apply(ratio, 1, entropy)`))
dim(sub_plot)
order_sub_plot <- sub_plot[1:500,] # or shannon_index > 2.98

df1 <- data.frame(peakID = rownames(order_sub_plot))
df1 <- str_split(df1$peakID,"_",simplify = T)
write.table(df1,"conservative_peak.bed",row.names = F,col.names = F,quote = F,sep="\t")

