library(tidyverse)
library(DESeq2)
library(WGCNA)

options(stringsAsFactors=F)

setwd("<YOUR_WORK_DIR>")

#***********************************************************************************************#
#                               Data input, cleaning and pre-processing                         #
#***********************************************************************************************#
Tissue <- c("blood)

for (i in 1:length(Tissue)) {
  raw_count <- readRDS("F:/sheep_adaptation/data/mature_raw_count_matrix.rds")
  raw_count[1:5,1:5]
  sample_info <- readRDS("F:/sheep_adaptation/data/all_sample_phenotype.rds")
  head(sample_info)
  table(sample_info$group)
  tissue_info <- readRDS("F:/sheep_adaptation/data/tissue_info.rds")
  head(tissue_info)
  Tissue <- unique(tissue_info$tissue)
  
  df <- filter(sample_info,tissue == Tissue[i])
  head(df)
  ## input expression matrix (raw counts from featureCount)
  traitData <- df[,c("sampleid","group")]
  head(traitData)
  traitData$group <- as.factor(traitData$group)
  countMatrix <- t(raw_count[rownames(df),])
  countMatrix[1:10,1:5]
  
  ## DEseq2 filter
  dds <- DESeqDataSetFromMatrix(countMatrix, colData = traitData, design = ~group)
  dds <- dds[rowSums(counts(dds)) >=10,] # filter sum(counts) < 1
  dds # summary of data
  vsd <- varianceStabilizingTransformation(dds)
  exprSet <- assay(vsd)  
  # head(exprSet)
  # detach("package:DESeq2")
  # output
  # saveRDS(exprSet,paste(tissue,"_WGCNA_input.csv",sep = ""))
  # write.csv(exprSet,paste(tissue,"_WGCNA_input.csv",sep = ""),row.names = T,quote = F)
  
  
  #***********************************************************************************************#
  #                                             WGCNA                                             #
  #***********************************************************************************************#
  ## option
  type = "signed" # (unsigned, signed, signed hybird)
  corType = "pearson"
  maxPOutliers = ifelse(corType=="pearson",1,0.05) 
  robustY = ifelse(corType=="pearson",T,F) 
  
  #***********************************************************************************************#
  #                                          clean data                                           #
  #***********************************************************************************************#
  ## step 1.1: filtration 
  m.mad <- apply(exprSet,1,mad) 
  dataExprVar <- exprSet[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
  dim(dataExprVar)
  dataExpr <- as.data.frame(t(dataExprVar)) 
  dim(dataExpr) 
  dataExpr[1:5,1:5]
  
  ## step 1.2 detect missing value
  gsg = goodSamplesGenes(dataExpr, verbose = 3) 
  gsg$allOK # TRUE: all genes have passed the cuts
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", 
                       paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", 
                       paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  rm(countMatrix,vsd,dds,dataExprVar)
  
  ## step 1.3: sample clustering and detect outliers
  #tiff("sampleTree.tiff",units = "in",res = 400 ,width = 18, height = 9)
  sampleTree = hclust(dist(dataExpr), method = "average")
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  abline(h=150,col="red",lwd=1.5)
  #dev.off()
  # cut outliers
  clust <- cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
  table(clust)
  keepSamples <- (clust==1)
  dataExpr <- dataExpr[keepSamples,]
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  dim(dataExpr)
  # head(dataExpr)[,1:8]
  # new traitData
  sample <- rownames(dataExpr)
  traitData <- traitData[sample,]
  save(dataExpr,traitData,file=paste("./data/",Tissue[i],"_WGCNA_input.RData",sep="")) 
  
  #***********************************************************************************************#
  #                                              WGCNA                                            #
  #***********************************************************************************************#
  ## calculate soft threshold
  powers <- c(c(1:10), seq(12,30, by=2))
  sft <- pickSoftThreshold(dataExpr, powerVector = powers, verbose=5, networkType=type)
  sft$powerEstimate  
  ## plot the results
  # pdf(paste(tissue,"_soft_value.pdf",sep=""),width = 10, height = 6)
  sizeGrWindow(9,5)
  par(mfrow = c(1,2)) 
  # par(mfrow = c(2,1)) # par():combine multiple pics
  cex1 = 0.8
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.85,col="steelblue",lwd=1.5)  # R^2 = 0.9
  # abline(h=0.9,col="steelblue",lwd=1.5)  # R^2 = 0.9
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  # dev.off()


#***********************************************************************************************#
#         Step-by-step construction of the gene network and identification of modules           #
#***********************************************************************************************#
#################################### Linux environment ####################################
cor <- WGCNA::cor
# enableWGCNAThreads(10) # multi-threads (useless under windows)
## step 2.1: calculate adjacency matrix 
adjacency = adjacency(dataExpr, power = sft$powerEstimate) # sft$powerEstimate
# adjacency = adjacency(dataExpr, power = 17) (if β=NA)
TOM = TOMsimilarity(adjacency) 
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average") 
save(geneTree,dissTOM,file=paste("./data/",Tissue[i],"_WGCNA_matrix.RData",sep="")) 



#***********************************************************************************************#
#                                   Associate with phenotype                                    #
#***********************************************************************************************#
cor <-stats::cor 
design1 <- model.matrix(~ traitData$group + 0) 
head(design1)

MEs0 <- moduleEigengenes(dataExpr,mergedColors)$eigengenes  
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs,design1, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

## plot module-trait relationship (correlations and their p-values)
modNames <- substring(names(MEs),3)
df1 <- table(mergedColors) %>% as.data.frame() 
rownames(df1) <- df1$mergedColors
df2 <- df1[paste(dput(modNames)),]

sizeGrWindow(12,12)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file = paste(Tissue[i],"_module_trait_cor.pdf",sep=""),width = 6.5, height = 9.5)
par(oma=c(1,1,1,1),mar = c(6, 10, 3, 1.2))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               # xLabels = unique(traitData$condition),
               yLabels = names(MEs), # color-blocks
               xLabels = paste(c("Day0","Day6","Day13","Day20","Mon8","Tibetan")),
               ySymbols = paste(modNames," (",df2$Freq,")",sep=""),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-Trait relationships"))
dev.off()
rm(dissTOM)

p.adj <- p.adjust(moduleTraitPvalue,method = "BH") 
padj_matrix <- matrix(data = p.adj, ncol = 6) 
colnames(padj_matrix) <- c("Day0","Day6","Day13","Day20","Mon8","Tibetan")
padj_matrix <- data.frame(Module=modNames,padj_matrix)
# write.table(MEs,paste(Tissue[i],"_moduleTrait_cor.txt",sep=""),sep = "\t",quote = F)
write.csv(padj_matrix,paste(Tissue[i],"_moduleTrait_FDR.csv",sep=""),quote = F,row.names = F)
rm(df1,df2)


# calculate ModuleMembership (MM)
geneModuleMembership <- as.data.frame(cor(dataExpr, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="")
# calculate geneTraitSignificance (GS)
geneTraitSignificance <- as.data.frame(cor(dataExpr, design1, use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(design1), sep="")
names(GSPvalue) <- paste("p.GS.", names(design1), sep="")
geneInfo <- cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
write.table(geneInfo, file=paste("./result/",Tissue[i],"_geneInfo.txt",sep=""), sep="\t", quote=F)


#***********************************************************************************************#
#                                  Select significant modules                                   #
#***********************************************************************************************#
module_gene <- list()
# Select module (here to output all modules)
for(j in 1:length(modNames)) {
  gene_cluster <- colnames(dataExpr)
  inModule <- (moduleColors == modNames[j])
  module_gene[[paste(modNames[j])]] <- gene_cluster[inModule] # store in list
  # write.csv(modgene,paste(tissue,"_",modNames[i],".csv",sep = ""),quote = F,row.names = F)
}
summary(module_gene)
save(module_gene,file =  paste("./result/",Tissue[i],"_cluster_gene.RData",sep=""))

rm(list=ls())

}



