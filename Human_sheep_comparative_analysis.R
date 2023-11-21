library(Seurat)
library(Rtsne)
library(tidyverse)

setwd("<SET_YOUR_PATH/>")

# ------------------------------------------------------------------------ #
#                          seurat remove batch effect                      #
# ------------------------------------------------------------------------ #
load("combine_meta_table.Rdata") # sample_id, gene_id
load("combine_tpm_matrix.Rdata") # x=sample_id, y= gene_id

expression <- CreateSeuratObject(t(MAT)) # target matrix format: row=gene_id,column=sample_id
expression@meta.data[["species"]] <-c (c(rep("Human",6792)),c(rep("Sheep",1277)))
expression@meta.data[["tissue_new"]] <- meta$tissue
expression <- SplitObject(expression,split.by = "species")
expression<-expression[c("Human","Sheep")]


for (i in 1:length(expression)) {
  expression[[i]] <- NormalizeData(expression[[i]], verbose = FALSE)
  expression[[i]] <- FindVariableFeatures(expression[[i]], selection.method = "vst",
                                          nfeatures = 17282, verbose = FALSE) #  nfeatures: orthologous gene number
}

expression <- FindIntegrationAnchors(object.list =expression, dims = 1:30,anchor.features = 17282)
expression <- IntegrateData(anchorset = expression, dims = 1:30)
## extract new expression matrix  ##
expression_inter <- expression@assays[["integrated"]]@data
expression_inter <- t(as.matrix(expression_inter))
dim(expression_inter)
expression_inter[1:5,1:5]
## output matrix  ##
save(expression_inter,file="integratedallgenes.Rdata")



# ------------------------------------------------------------------------ #
#                              t-sne clustering                            #
# ------------------------------------------------------------------------ #
# study: https://cloud.tencent.com/developer/article/1556202
# study:https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf

load("combine_meta_table.Rdata")
load("integratedallgenes.Rdata")

# meta <- mutate(meta,species=case_when(grepl("^GTEX",sampleid) ~ "Human",
#                                       grepl("^S",sampleid) ~ "Sheep"))
table(meta$species)
expression_tsne <- Rtsne(expression_inter,dims = 2, perplexity=30, theta=0.4,verbose=TRUE,
                         max_iter = 2000,check_duplicates = FALSE,partial_pca = T,num_threads=20)
# expresssion_inter: row=sample_id, column=gene_id

tsne <- as.data.frame(expression_tsne$Y)
colnames(tsne) <- c("tSNE1","tSNE2")
tsne$tissue <- meta$tissue 
tsne$species <- meta$species
class(tsne)
save(tsne,file = "human_sheep_tsne.RData")

# ------------------------------------------------------------------------ #
#                               visualize t-sne                            #
# ------------------------------------------------------------------------ #
load("human_sheep_tsne.RData")
## based on species
p1 <- ggplot(tsne,aes(tSNE1,tSNE2,color=species)) + 
  geom_point(size=1.8, alpha= 0.4) + 
  labs(x="t-SNE1", y="t-SNE2") +
  scale_color_manual(values = c("red","blue")) +
  theme_bw()+
  theme(plot.title = element_text(size = 10, face = "bold"),
        text = element_text(size = 10),
        legend.title = element_text(size=9, face="bold"), 
        legend.text = element_text(size=8, face="bold"), 
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11,face = "bold"),
        axis.text.y=element_text(size = 11,face = "bold"),
        panel.grid=element_blank(),
        plot.margin = margin(t=1,r=1,b=1,l=1,unit = "cm"))
p1
pdf("human_sheep_species_clustering.pdf",width = 5,height = 4)
p1
dev.off()


mycolor <- c("Pituitary"="darkorange1","Hypothalamus"="#ef8737","Cerebellum"="#ffb242","Cerebrum"="#ffd353","Stomach"="yellowgreen",  
             "Abomasum"="#275024", "Duodenum"="#47632a", "Jejunum"="#748f46", "Ileum"="#98ab76","Colon"="#ced1af",
             "Leukocyte"="#dc322f", "Spleen"="#f05b43","Heart" = "#d33682", "Muscle"="#9f2d55","Artery" = "#CD5B45",
             "Adipose"="darkcyan","Lung"="#FDAF9199", "Liver"="#6c71c4", "Kidney"="#016392","Rumen" = "#B4EEB4")

Tissue <- c("Adipose", "Artery","Cerebellum", "Colon", 
            "Heart", "Hypothalamus","Ileum", "Kidney", 
            "Leukocyte", "Liver", "Lung", "Muscle", 
            "Pituitary", "Spleen")

p4 <- tsne %>% filter(tissue %in% Tissue) %>% 
  ggplot(aes(tSNE1,tSNE2,color=tissue)) + 
  geom_point(size=1.8, alpha= 0.5) + 
  labs(x=expression(paste(italic(t),"-SNE1")),
       y=expression(paste(italic(t),"-SNE2"))) +
  scale_color_manual(values = mycolor) +
  theme_bw()+
  theme(plot.title = element_text(size = 10, face = "bold"),
        text = element_text(size = 10),
        legend.title = element_text(size=9, face="bold"), 
        legend.text = element_text(size=8, face="bold"),  
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold"),
        axis.text.y=element_text(size = 10,face = "bold"),
        panel.grid=element_blank(),
        plot.margin = margin(t=1,r=1,b=1,l=1,unit = "cm"))
p4

pdf("human_sheep_tissue_clustering_v2.pdf",width = 5.5,height = 4)
p4
dev.off()
