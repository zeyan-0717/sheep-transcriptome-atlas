library(tidyverse)

# set dir
setwd("<YOUR_WORK_DIR>")

## step 1: generate expression matrix (raw_count & tpm)
# extract file path
myfile1 <- list.files("<YOUR_PATH_TO_DATA/raw_count_matrix>")
head(myfile1)
myfile2 <- list.files("<YOUR_PATH_TO_DATA/tpm_matrix>")
head(myfile2)
# extract tissue 
tissue <- substr(myfile1,1,nchar(myfile1)-nchar("_raw_count_matrix.txt"))

# merge all expression files
DF <- NA # raw count
df <- NA # tpm
for (i in 1:length(tissue)) {
  # raw count
  DF1 <- data.table::fread(paste("<YOUR_PATH_TO_DATA/raw_count_matrix>",myfile1[i],sep=""),header = T)
  dim(DF1)
  DF2 <- DF1[,-1] %>% as.matrix()
  DF3 <- t(DF2)
  DF <- rbind(DF,DF3)
  # tpm
  df1 <- data.table::fread(paste("<YOUR_PATH_TO_DATA/tpm_matrix>/",myfile2[i],sep=""),header = T)
  dim(df1)
  df2 <- df1[,-1] %>% as.matrix()
  df3 <- t(df2)
  df <- rbind(df,df3)
  
}
dim(DF)
DF[1:10,1:5]
class(DF)
DF <- DF[-1,]
DF[1:10,1:5]
colnames(DF) <- DF1$Geneid
DF[1:10,1:5]
raw_count <- DF

dim(df)
df <- df[-1,]
df[1:5,1:5]
dim(df)
colnames(df) <- df1$geneid
class(df)
tpm <- df
# save expression matrix
save(raw_count,file = "<YOUR_PATH>/all_raw_count_matrix.RData")
save(tpm,file = "<YOUR_PATH>/all_tpm_matrix.RData")
rm(DF,DF1,DF2,DF3,df,df1,df2,df3,raw_count);gc()


## step 2: generate sample phenotype table (option)
info <- data.frame(sampleid=rownames(tpm))
# df1 <- readxl::read_excel("<YOUR_PATH>/tissue_abrev.xlsx")
df1 <- str_split(info$sampleid,"_",simplify = T) %>% as.data.frame()
head(df1)
names(df1) <- c("id","tissueID")
info <- cbind(info,df1)
head(info)
rm(df1)
 # group info
# mature individual
g1 <- c("S01","S02","S03","S04","S05", "S06","S07","S08","S09","S10") # 0 day
g2 <- c("S11","S12","S13","S14","S15", "S16","S17","S18","S19","S20") # 7 day
g3 <- c("S21","S22","S23","S24","S25", "S26","S27","S28","S29","S30") # 14 day
g4 <- c("S31","S32","S33","S34","S35", "S36","S37","S38","S39","S40") # 21 day
g5 <- c("S41","S42","S43","S44","S45", "S46","S47","S48","S49","S50") # tibetan
g6 <- c("S81","S82","S83","S84","S85","S86","S87","S88","S89","S90") # 8 month
# offspring
g7 <- c("S51","S52","S53","S54","S55","S56") # hu sheep born in plain
g8 <- c("S61","S62","S63","S64","S65","S66") # hu sheep born in plateau
g9 <- c("S71","S72","S73","S74","S75","S76") # tibetan

tb <- c(g5,g9)
hu <- c(g1,g2,g3,g4,g6,g7,g8)
f0 <- c(g1,g2,g3,g4,g5,g6)
f1 <- c(g7,g8,g9)
female <- c(f0,"S51","S55","S56","S61","S62","S63","S64","S65","S71","S72","S73","S74")
male <- c("S52","S53","S54","S66","S75","S76")


# prepare meta table (option)
df4 <- info %>% mutate(group = case_when(id %in% g1 ~ "D0_Hu",id %in% g2 ~ "D06_Hu",
                                         id %in% g3 ~ "D13_Hu",id %in% g4 ~ "D20_Hu",
                                         id %in% g5 ~ "Tibetan",id %in% g6 ~ "M08_Hu",
                                         id %in% g7 ~ "Hu_plain_lamb",
                                         id %in% g8 ~ "Hu_plateau_lamb",
                                         id %in% g9 ~ "Tibetan_lamb"),
                       tissue = case_when(tissueID == "BJ" ~ "Muscle",
                                          tissueID == "CT" ~ "Pituitary",
                                          tissueID == "DN" ~ "Cerebrum",
                                          tissueID == "heart" ~ "Heart",
                                          tissueID == "JC" ~ "Colon",
                                          tissueID == "KC" ~ "Jejunum",
                                          tissueID == "kidney" ~ "Kidney",
                                          tissueID == "liver" ~ "Liver",
                                          tissueID == "lung" ~ "Lung",
                                          tissueID == "LW" ~ "Rumen",
                                          tissueID == "QN" ~ "Hypothalamus",
                                          tissueID == "spleen" ~ "Spleen",
                                          tissueID == "WBC" ~ "Leukocyte",
                                          tissueID == "XN" ~ "Cerebellum",
                                          tissueID == "ZC" ~ "Duodenum",
                                          tissueID == "ZF" ~ "Adipose",
                                          tissueID == "ZW" ~ "Abomasum",
                                          tissueID == "XDM" ~ "Artery",
                                          tissueID == "HC" ~ "Ileum"),
                       generation = case_when(id %in% f0 ~ "F0",
                                              id %in% f1 ~ "F1"),
                       breed = case_when(id %in% hu ~ "Hu",id %in% tb ~ "Tibetan"),
                       type = case_when(tissue %in% c("Pituitary") ~ "Pituitary",
                                        tissue %in% c("Cerebellum","Cerebrum","Hypothalamus") ~ "CNS",
                                        tissue %in% c("Leukocyte","Spleen") ~ "Blood/Immune",
                                        tissue %in% c("Adipose") ~ "Adipose",
                                        tissue %in% c("Muscle","Heart") ~ "Muscle/Heart",
                                        tissue %in% c("Kidney") ~ "Kidney",
                                        tissue %in% c("Liver") ~ "Liver",
                                        tissue %in% c("Lung") ~ "Lung",
                                        tissue %in% c("Artery") ~ "Artery",
                                        tissue %in% c("Rumen") ~ "Rumen",
                                        tissue %in% c("Abomasum","Ileum","Colon","Jejunum","Duodenum") ~ "Intestine"),
                       sex = case_when(id %in% female ~ "Female",
                                       id %in% male ~ "Male"),
                       age = case_when(id %in% f0 ~ "1.5year",
                                       id %in% g7 ~ "1month",
                                       id %in% c(g8,g9) ~ "2month"))

head(df4)
dim(df4)
rownames(df4) <- df4$sampleid
saveRDS(df4, file = "<YOUR_PATH>/all_sample_phenotype.rds")
# tissue info
tissue_info <- distinct(meta[,c("tissueID","tissue","type")])
saveRDS(tissue_info, file = "<YOUR_PATH>/tissue_info.rds")


## step 3: t-SNE clustering all samples
library(Rtsne)

set.seed(1234)
tsne_res = Rtsne(tpm, dims = 2, pca = T, max_iter = 2000,theta = 0.4,perplexity = 30,verbose = F) # iteration number can change(max_iter)
tsne_res$itercosts

## extract t-SNE result
info <- df4
tsne <- as.data.frame(tsne_res$Y)
colnames(tsne) <- c("tSNE1","tSNE2")
tsne$Tissue <- info$tissue 
tsne$Group <- info$group
tsne$Type <- info$type
tsne$Generation <- info$generation
class(tsne)
saveRDS(tsne,file = "<YOUR_PATH>/all_sample_tsne.rds")

## visualize
# method 1: based on tissue type
color1 <- c("Pituitary"="darkorange1","CNS"="#EEAD0E", "Intestine"="#228B22",
            "Blood/Immune"="#dc322f","Artery" = "#CD5B45","Muscle/Heart"="#d33682",
            "Adipose"="darkcyan","Kidney"="#016392",
            "Liver"="#6c71c4","Lung"="#FDAF9199","Rumen" = "#B4EEB4")

color2 <- c("Pituitary"="darkorange1","Hypothalamus"="#ef8737","Cerebellum"="#ffb242","Cerebrum"="#ffd353",  
                        "Abomasum"="#275024", "Duodenum"="#47632a", "Jejunum"="#748f46", "Ileum"="#98ab76","Colon"="#ced1af",
                        "Leukocyte"="#dc322f", "Spleen"="#f05b43","Heart" = "#d33682", "Muscle"="#9f2d55","Artery" = "#CD5B45",
                        "Adipose"="darkcyan","Lung"="#FDAF9199", "Liver"="#6c71c4", "Kidney"="#016392","Rumen" = "#B4EEB4")
# ggplot2
p1 <- ggplot(tsne,aes(tSNE1,tSNE2,color=Type,shape=Generation)) + 
  geom_point(size=2.7) + 
  labs(x="t-SNE1", y="t-SNE2") +
  scale_color_manual(values = color1) +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        legend.title = element_text(size=13, face="bold"), 
        legend.text = element_text(size=13, face="bold"),  
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13,face = "bold"),
        axis.text.y=element_text(size = 13,face = "bold"),
        panel.grid=element_blank(),
        plot.margin = margin(t=3,r=3,b=3,l=3,unit = "cm"))
p1

pdf("<YOUR_PATH/all_sample_tissuetype_tsne1.pdf>",height = 6 ,width = 8)
p1
dev.off()

# method 2: based on tissue
p2 <- ggplot(tsne,aes(tSNE1,tSNE2,color= Tissue,shape=Generateion)) + 
  geom_point(size=2.7,alpha= 0.7) + 
  labs(x="t-SNE1", y="t-SNE2") +
  scale_color_manual(values = color2) +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        legend.title = element_text(size=13, face="bold"),
        legend.text = element_text(size=13, face="bold"),  
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13,face = "bold"),
        axis.text.y=element_text(size = 13,face = "bold"),
        panel.grid=element_blank(),
        plot.margin = margin(t=4,r=4,b=4,l=4,unit = "cm")) # legend.position = "bottom"
p2

pdf("<YOUR_PATH/all_sample_tissue_tsne2.pdf>",height = 7.5 ,width = 9.5)
p2
dev.off()
