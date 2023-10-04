library(tximport)
library(sva)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(WGCNA)
library(rrcov)
library(gplots)
library(stringr)
library(ggprism)

# Import metadata
metadata <- read.csv("metadata.csv",header=T)
metadata$phenotype <- as.factor(as.character(metadata$phenotype))
metadata$block <- as.factor(as.character(metadata$block))
metadata$lineage <- as.factor(as.character(metadata$lineage))

# Import counts
tx2gene <- read.csv("tx2gene.csv",header=T)[,c(2,3)]
dirs.kallisto <- list.dirs("kallisto", recursive=FALSE)
files.list <- list()
for (i in 1:length(dirs.kallisto)){
  files.list[i] <- paste(unlist(dirs.kallisto[i]),"abundance.tsv",sep="/")
}

txi <- tximport(unlist(files.list), type = "kallisto",
                tx2gene = tx2gene,geneIdCol="gene_id",txIdCol="tx_id")

datExpr <- as.data.frame(txi[["counts"]])
dirnames.l <- length(unlist(strsplit(unlist(dirs.kallisto),"/")[1]))
dirnames <- sapply(strsplit(unlist(dirs.kallisto),"/"), "[", dirnames.l)
names(datExpr) <- dirnames
write.csv(datExpr,"datExpr.csv")


# Differential expression analysis (all blocks)
metadata_new <- metadata
metadata_new$BL <- factor(paste0(metadata_new$block,metadata_new$lineage))
datExpr_new <- datExpr[,names(datExpr)%in%metadata_new$individual]

dds <- DESeqDataSetFromMatrix(countData=round(datExpr_new),
                              colData=metadata_new,
                              design=~BL+phenotype)
DESeq <- DESeq(dds)
vsd <- varianceStabilizingTransformation(DESeq, blind=F)
dds.transform <- data.frame(assay(vsd))

### outlier check
checkpca <- PcaGrid(t(dds.transform))
plot(checkpca)

cpca.plot <- data.frame(sd=checkpca$sd,index=1:length(checkpca$sd),
                        flag=as.character(checkpca$flag),
                        label=names(checkpca$flag))
g <- ggplot(cpca.plot, aes(x=index, y=sd)) + 
  geom_point(size=2) + ggtitle("Robust PCA") +
  theme_prism() + xlab("Index") + ylab("Score distance") +
  scale_x_continuous(breaks=c(seq.int(0,40,5))) +
  scale_y_continuous(breaks=c(seq.int(2,8,1)),limits=c(2,8)) +
  geom_hline(yintercept=checkpca@cutoff.sd,linetype="dashed", color = "red") +
  geom_text(aes(label=ifelse(flag==F,as.character(label),'')),
            hjust=.5,vjust=-1,size=5,fontface="bold")
g
ggsave("figS1.png",g,width=7,height=5)


metadata_new <- metadata_new[-c(17),]
datExpr_new <- datExpr[,names(datExpr)%in%metadata_new$individual]

dds <- DESeqDataSetFromMatrix(countData=round(datExpr_new),
                              colData=metadata_new,
                              design=~BL+phenotype)
DESeq <- DESeq(dds)
vsd <- varianceStabilizingTransformation(DESeq, blind=F)
dds.transform <- data.frame(assay(vsd))
checkpca <- PcaGrid(t(dds.transform))
plot(checkpca)

PCA <- plotPCA(vsd, intgroup=c("phenotype","BL"))
ylab <- PCA[["labels"]][["y"]]
xlab <- PCA[["labels"]][["x"]]
PCA <- PCA[["data"]]
names(PCA)[6] <- "individual"
PCA <- PCA[match(metadata_new$individual, PCA$individual),]
PCA <- left_join(PCA,metadata_new[,c(1,3)],by="individual")

PCA.plot <- ggplot(PCA, aes(x=PC1, y=PC2)) +
  geom_point(size=3, aes(color=phenotype,shape=BL)) +
  xlab(xlab)+
  ylab(ylab) +
  theme(text = element_text(size=20)) +
  xlim(c(-100,100))
PCA.plot

results <- data.frame(results(DESeq, contrast=c("phenotype",
                                                "unresponsive",
                                                "responsive")))
results$abslog2FoldChange <- abs(results$log2FoldChange)
write.csv(results,"DESeq2_results.csv")
write.csv(metadata_new,"metadata_noOutliers.csv",row.names=F)
write.csv(dds.transform,"datExpr_vst.csv")


### WGCNA
##### blocks 1-3 | signed 15 | 2 modules correlated with phenotype
metadata <- read.csv("metadata_noOutliers.csv",header=T)
metadata$phenotype <- factor(as.character(metadata$phenotype))
metadata$block <- factor(as.character(metadata$block))
metadata$lineage <- factor(as.character(metadata$lineage))
metadata$BL <- factor(paste0(metadata$block,metadata$lineage))

DEG_results <- read.csv("DESeq2_results.csv",row.names=1)

datExpr_WGCNA <- read.csv("datExpr_vst.csv",row.names=1)
datExpr_WGCNA <- datExpr_WGCNA[row.names(datExpr_WGCNA)%in%
                     row.names(DEG_results[!is.na(DEG_results$padj),]),]
datExpr_WGCNA <- datExpr_WGCNA[,match(metadata$individual, names(datExpr_WGCNA))]
metadata <- metadata[metadata$individual%in%names(datExpr_WGCNA),]

datExpr_SVM <- round(datExpr[row.names(datExpr)%in%row.names(datExpr_WGCNA),
                       names(datExpr)%in%names(datExpr_WGCNA)])
datExpr_SVM <- datExpr_SVM[,match(metadata$individual, names(datExpr_SVM))]
write.csv(datExpr_SVM,"datExpr_SVM.csv")


powers = c(c(10:20))
sft = pickSoftThreshold(t(datExpr_WGCNA),
                        powerVector = powers, 
                        verbose = 5,
                        networkType = "signed")

png("figS3.png",width=7,height=5,units="in",res=300)
par(mfrow = c(1,2))
cex1 = 0.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels=powers, cex=cex1,col="red")
dev.off()


adj_matrix <- adjacency(t(datExpr_WGCNA), power = 15, type = "signed")
TOM = TOMsimilarity(adj_matrix)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
MEList = moduleEigengenes(t(datExpr_WGCNA), colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.25
merge = mergeCloseModules(t(datExpr_WGCNA), dynamicColors, 
                          cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

png("figS4.png",width=7,height=5,units="in",res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

unique(mergedColors)

modules <- data.frame(geneID=paste0("LOC",row.names(datExpr_WGCNA)),
                      module=mergedColors)


datTraits <- metadata[,c(1,2,3,6)]
datTraits <- datTraits[datTraits$individual%in%names(datExpr_WGCNA),]
datTraits$phenotype <- as.character(datTraits$phenotype)
datTraits$block <- as.character(datTraits$block)
datTraits$BL <- as.character(datTraits$BL)
row.names(datTraits) <- datTraits$individual
datTraits$individual <- NULL
datTraits[datTraits$BL=="1A","BL"]=0
datTraits[datTraits$BL=="1B","BL"]=1
datTraits[datTraits$BL=="2A","BL"]=2
datTraits[datTraits$BL=="2B","BL"]=3
datTraits[datTraits$BL=="3A","BL"]=4
datTraits[datTraits$BL=="3B","BL"]=5
datTraits[datTraits$block==1,"block"]=0
datTraits[datTraits$block==2,"block"]=1
datTraits[datTraits$block==3,"block"]=2
datTraits[datTraits$phenotype=="unresponsive","phenotype"]=0
datTraits[datTraits$phenotype=="responsive","phenotype"]=1
datTraits$BL <- factor(datTraits$BL)
datTraits$block <- factor(datTraits$block)
datTraits$phenotype <- factor(datTraits$phenotype)


nGenes = ncol(t(datExpr_WGCNA))
nSamples = nrow(t(datExpr_WGCNA))
MEs0 = moduleEigengenes(t(datExpr_WGCNA), mergedColors)$eigengenes
MEs = orderMEs(MEs0)
# MEs.NG = length(names(MEs))-1
MEs.NG = length(names(MEs))
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
weight = as.data.frame(datTraits$phenotype)
names(weight) = "weight"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(t(datExpr_WGCNA), MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(t(datExpr_WGCNA), weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")


textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor[1:MEs.NG,],
               xLabels = names(datTraits),
               yLabels = names(MEs)[1:MEs.NG],
               ySymbols = names(MEs)[1:MEs.NG],
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[1:MEs.NG,],
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               legendLabel = expression(italic("r")),
               main = paste("Module-trait relationships"))

write.csv(modules,"WGCNA_modules.csv",row.names=F)