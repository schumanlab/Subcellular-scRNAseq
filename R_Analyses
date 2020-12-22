# Subcellular-scRNAseq

#1--------------------------------------------------------Loading Analysis---------------------------------------------------------------------
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(scater)
library(future)
library(Seurat)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(MASS)
library("reshape")
library(ISLR)
library(dunn.test)
library(ggpubr)
library(ggrepel)

options(future.globals.maxSize = 1500 * 1024^2)

#Assigning links
setwd()
Count_table<-"All_Samples.txt"
Samples_Table <- "Samples_Table.txt"

#Loading files
countfile<-read.delim(Count_table, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
metadata <- read.delim(Samples_Table, check.names = FALSE, header = TRUE)

#2--------------------------------------------------------Formating Files---------------------------------------------------------------------
countfile <- countfile[rowSums(countfile[, c(3:962)])>10, ]
GeneIDs <- row.names(countfile)
GeneIDs <- as.data.frame(GeneIDs)
GeneIDs <- separate(GeneIDs,"GeneIDs",into = c("Name","ID"),sep = "_XLOC_")
GeneIDs <- GeneIDs$ID
GeneNames <-countfile$GeneName
GeneNames <-as.character(GeneNames)
Coordinates <-countfile$Coordinates
Coordinates <-as.character(Coordinates)
Coordinates2 <-as.data.frame(Coordinates)

countfile <- dplyr::select(countfile, -1,-2)
countfile <- as.matrix(countfile)
Coordinates2 <- separate(Coordinates2,"Coordinates",into = c("Chromosome","Coordinates"),sep = ":")
CHR <-Coordinates2$Chromosome

sce <- SingleCellExperiment(list(counts=countfile))
rowData(sce)$GeneIDs <- GeneIDs
rowData(sce)$GeneNames <- GeneNames
rowData(sce)$Coordinates <- Coordinates
rowData(sce)$CHR <- CHR
mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
metadata$pct_counts_Mt<-sce$pct_counts_Mt

n <- match(colnames(sce), metadata[["SampleName"]])
stopifnot(all(!is.na(n)))
metadata <- metadata[n,]

rownames(sce) <- uniquifyFeatureNames(rowData(sce)$GeneIDs, rowData(sce)$GeneNames)
rownames(countfile) <- rownames(sce)

countfile <- countfile[-mito, ]#remove mitochondrial genes

#3--------------------------------------------------------Preparing Spike-Ins---------------------------------------------------------------------
ERCC.index <- grep(pattern = "^ERCC-", x = rownames(countfile), value = FALSE)
Total.ERCC <- Matrix::colSums(countfile[ERCC.index, ])
percent.ERCC <- Matrix::colSums(countfile[ERCC.index, ])/Matrix::colSums(countfile)
ERCC_Molecules_Total_Input <-311691
detected.fraction.ERCC <- ((Total.ERCC)/(Total.ERCC+ERCC_Molecules_Total_Input))
detected.perc.ERCC <- (detected.fraction.ERCC)*100

metadata$Total.ERCC <-Total.ERCC
metadata$percent.ERCC <-percent.ERCC
metadata$detected.perc.ERCC <-detected.perc.ERCC
countfile <- countfile[-ERCC.index, ]

rRNA_baits <- "Rn18|Rn45|rRNA|Rn5-8"
rRNA.index <- grep(pattern = rRNA_baits, x = rownames(countfile), value = FALSE)
rRNA <- Matrix::colSums(countfile[rRNA.index, ])
metadata$rRNA <- rRNA
countfile <- countfile[-rRNA.index, ]

#4--------------------------------------------------------Calculate QC Metrics---------------------------------------------------------------------
metadata$total_counts<-sce$total_counts
metadata$log_total_counts<-log(sce$total_counts)
predicted.sample.total <- ((sce$total_counts)/(detected.fraction.ERCC))
metadata$predicted.sample.total<-predicted.sample.total
metadata$log.predicted.sample.total<-log(predicted.sample.total)

#5--------------------------------------------------------Transfering to Seurat and adding metada---------------------------------------------------------------------
Seurat_All <- CreateSeuratObject(countfile, project = "All")
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$SampleName, col.name = 'SampleName')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Index, col.name = 'Index')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Neuron, col.name = 'Neuron')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Comparment, col.name = 'Comparment')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$ComparmentAv, col.name = 'ComparmentAv')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Preparation, col.name = 'Preparation')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Plate, col.name = 'Plate')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Collection, col.name = 'Collection')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$SequencingRun, col.name = 'SequencingRun')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$pct_counts_Mt, col.name = 'pct_counts_Mt')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$rRNA, col.name = 'rRNA')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$percent.ERCC, col.name = 'percent.ERCC')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Total.ERCC, col.name = 'Total.ERCC')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$detected.perc.ERCC, col.name = 'detected.perc.ERCC')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$predicted.sample.total, col.name = 'predicted.sample.total')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$log.predicted.sample.total, col.name = 'log.predicted.sample.total')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$total_counts, col.name = 'total_counts')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$log_total_counts, col.name = 'log_total_counts')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$GraphOrder, col.name = 'GraphOrder')
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata$Source, col.name = 'Source')

#6--------------------------------------------------------QC Filter 1: Failed Samples---------------------------------------------------------------------
is.mix <- grepl("^E2_S1_I3", Seurat_All@meta.data$orig.ident) #removal of sample known to be contaminated by other somata
Seurat_All <- Seurat_All[,!is.mix]
is.Empty <- grepl("^Empty", Seurat_All@meta.data$Comparment)
Seurat_All <- Seurat_All[,!is.Empty]
metadata.all <- Seurat_All[[c("nCount_RNA","nFeature_RNA")]] 
metadata.all$FeaCount <- (metadata.all$nFeature_RNA)/(metadata.all$nCount_RNA)
Seurat_All <- AddMetaData(object = Seurat_All, metadata = metadata.all$FeaCount, col.name = 'FeatCount')

#7------------------------------------------Subseting Samples---------------------------------------------
QC_Som <- subset(Seurat_All, Comparment=="Somata")
QC_Den <- subset(Seurat_All, Comparment=="Dendrites")
QC_NC <- subset(Seurat_All, Comparment=="NegativeControls")

Somata_Metadata <- QC_Som@meta.data
Dendrites_Metadata <- QC_Den@meta.data
NC_Metadata <- QC_NC@meta.data
Som_Avg_Counts <- mean(Somata_Metadata$nCount_RNA)
Den_Avg_Counts <- mean(Dendrites_Metadata$nCount_RNA)
NC_Avg_Counts <- mean(NC_Metadata$nCount_RNA)
Som_Avg_Genes <- mean(Somata_Metadata$nFeature_RNA)
Den_Avg_Genes <- mean(Dendrites_Metadata$nFeature_RNA)
NC_Avg_Genes <- mean(NC_Metadata$nFeature_RNA)

is.Somata <- grepl("^Somata", Seurat_All@meta.data$Comparment)
QC_Den_NC <- Seurat_All[,!is.Somata]
is.NC <- grepl("^NegativeControls", Seurat_All@meta.data$Comparment)
QC_Som_Den <- Seurat_All[,!is.NC]

#8------------------------------------------QC Graphs---------------------------------------------
VlnPlot(Seurat_All, features = c("nCount_RNA","nFeature_RNA"), cols = c("#66A61E","#7570B3","gray"), group.by = "GraphOrder", 
        log = TRUE, pt.size = 2)
Fig1A<-VlnPlot(Seurat_All, features = c("nFeature_RNA"), cols = c("#66A61E","#7570B3","gray"), group.by = "GraphOrder", log = TRUE) + 
  NoLegend() + theme(axis.text.y = element_text(size = 16))
VlnPlot(Seurat_All, features = c("rRNA", "pct_counts_Mt"), cols = c("#66A61E","#7570B3","gray"), group.by = "GraphOrder", 
        pt.size = 2, log = TRUE) + NoLegend() + theme(axis.text.y = element_text(size = 16))

FigS1EL <- FeatureScatter(QC_Som, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = c("#66A61E"), pt.size = 2, group.by = "Comparment") + NoLegend() + geom_hline(yintercept = 700, size = 1, colour = "#FF3721") + ylim(0,15000) + xlim(0,500000)
FigS1EM <- FeatureScatter(QC_Den, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = c("#7570B3"), pt.size = 2, group.by = "Comparment") + NoLegend() + geom_hline(yintercept = 55, size = 1, colour = "#FF3721") + ylim(0,15000) + xlim(0,500000)
FigS1ER <- FeatureScatter(QC_NC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = c("gray"), pt.size = 2, group.by = "Comparment") + NoLegend() + ylim(0,15000) + xlim(0,500000)
plot_grid(FigS1EL,FigS1EM,FigS1ER, ncol = 3)

#9------------------------------------------QC Filter 2: Number of Features---------------------------------------------
Seurat_Fourth <- subset(Seurat_All, SequencingRun=="Fourth")
Seurat_Fourth <- subset(Seurat_Fourth, subset = nFeature_RNA > 700)
is.Fourth <- grepl("^Fourth", Seurat_All@meta.data$SequencingRun)
Seurat_First <- Seurat_All[,!is.Fourth]
is.NC <- grepl("^NegativeControls", Seurat_First@meta.data$Comparment)
Seurat_First <- Seurat_First[,!is.NC] #removal of negative control samples
is.I5 <- grepl("^I5", Seurat_First@meta.data$Index)
Seurat_First <- Seurat_First[,!is.I5] #Removal of samples of Index 5 since they tended to have to few reads

hist(Seurat_First$Total.ERCC, breaks = 50, col = "gray", cex.axis=1.3) + abline(v=10000, col="red")
Seurat_First <- subset(Seurat_First, subset = Total.ERCC >10000 & nFeature_RNA > 55)#removal of samples in which less than 10K ERCCs where detected and/or less than 55 RNA species
Seurat_All <- merge(x = Seurat_First, y = Seurat_Fourth)

Seurat_Den <- subset(Seurat_All, Comparment=="Dendrites")
Seurat_Som <- subset(Seurat_All, Comparment=="Somata")
Seurat_Som <- subset(Seurat_Som, subset = nFeature_RNA >700)#removal of somata with less 700 RNA species
Seurat_All <- merge(x = Seurat_Som, y = Seurat_Den)

#10--------------------------------------------------------Removal of contaminated/unhealthy cells---------------------------------------------------------------------
Seurat_All <- SCTransform(Seurat_All, variable.features.n = 3000)
Seurat_All <- RunPCA(Seurat_All)
Seurat_All <- RunUMAP(Seurat_All, dims = 1:30)
Seurat_All <- FindNeighbors(Seurat_All, dims = 1:2, verbose = FALSE, reduction = "umap")
Seurat_All <- FindClusters(Seurat_All, verbose = FALSE)

plan("multiprocess", workers = 4)
Seurat_All.markers <- FindAllMarkers(Seurat_All, min.pct = 0.5, logfc.threshold = 0.5, only.pos = TRUE, test.use = "LR")
top10 <- Seurat_All.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(Seurat_All, features = top10$gene, label = TRUE, angle = 45, hjust = 0, size = 4) + NoLegend()

new.cluster.ids <- c("Dendrites","Dendrites","GABAergic","GABAergic","GABAergic","Glutamatergic","GABAergic",
                     "Glutamatergic","Apoptotic","Glia-contaminated") #IMPORTANT: groups and order may need to be updated to obtain same results
names(new.cluster.ids) <- levels(Seurat_All)
Seurat_All <- RenameIdents(Seurat_All, new.cluster.ids)
Seurat_All <- AddMetaData(object = Seurat_All, metadata = Seurat_All@active.ident, col.name = 'QC_type')

FigS1F <- DimPlot(Seurat_All, reduction = "umap", pt.size = 2, cols = c("#7570B3","deepskyblue","coral3","gray65","gray65"), label= FALSE, label.size = 4) + NoLegend()
FigS1G <- FeaturePlot(Seurat_All, features = c("Gfap","S100b","Fabp7","Vim"), ncol = 2, pt.size = 1) + NoLegend()

#11--------------------------------------------------------Removing Apoptotic & Glia Contamintaed---------------------------------------------------------------------
is.apoptotic <- grepl("^Apoptotic", Seurat_All@active.ident)
Seurat_Cleaned <- Seurat_All[,!is.apoptotic]
is.Glial <- grepl("^Glia-contaminated", Seurat_Cleaned@active.ident)
Seurat_Cleaned <- Seurat_Cleaned[,!is.Glial]

#12--------------------------------------------------------Somata Clustering--------------------------------------------------------------------
Seurat_Som <- subset(Seurat_Cleaned, Comparment=="Somata")
Seurat_Som <- SCTransform(Seurat_Som, variable.features.n = 3000, return.only.var.genes = FALSE)
Seurat_Som <- RunPCA(Seurat_Som, verbose = FALSE)
Seurat_Som <- RunUMAP(Seurat_Som, dims = 1:30, verbose = FALSE)
Seurat_Som <- FindNeighbors(Seurat_Som, dims = 1:15, verbose = FALSE)
Seurat_Som <- FindClusters(Seurat_Som, verbose = FALSE)
DimPlot(Seurat_Som, reduction = "umap", label = TRUE, pt.size = 4, label.size = 7) + NoLegend()
table(Idents(Seurat_Som))

#13--------------------------------------------------------IDing Somata clusters---------------------------------------------------------------------
new.cluster.ids <- c("GABAergic_1","GABAergic_2","GABAergic_3","Glutamatergic_1","Glutamatergic_2")#IMPORTANT: groups and order may need to be updated to obtain same results
names(new.cluster.ids) <- levels(Seurat_Som)
Seurat_Som <- RenameIdents(Seurat_Som, new.cluster.ids)
Seurat_Som <- AddMetaData(object = Seurat_Som, metadata = Seurat_Som@active.ident, col.name = 'Specific_type')
my_levels <- c("Glutamatergic_1","Glutamatergic_2","GABAergic_1","GABAergic_2","GABAergic_3") 
Idents(Seurat_Som) <- factor(Idents(Seurat_Som), levels = my_levels)
Fig1D <- DimPlot(Seurat_Som, reduction = "umap", pt.size = 5, cols = c("coral2","coral4","#1B9E77","#E6AB02","#386CB0")) + NoLegend()

Seurat_Som.markers <- FindAllMarkers(Seurat_Som, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")
write.table(Seurat_Som.markers, file="Somata_Clusters_Differential_Expression.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)
top8 <- Seurat_Som.markers %>% group_by(cluster) %>% top_n(n = -8, wt = p_val_adj)
gene.list <- top8$gene
FigS2D <- DoHeatmap(Seurat_Som, features = gene.list, label = FALSE, group.bar = TRUE,
          group.colors = c("coral2","coral4","#1B9E77","#E6AB02","#386CB0"), disp.max = 10, disp.min = 1.5) + 
  NoLegend() + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis(option = "cividis")

Fig1E <- FeaturePlot(Seurat_Som, features = c("Slc6a1","Cort","Id2","Cck","Sv2b","Zranb2"), pt.size = 0.5, 
            ncol = 3, slot = "scale.data")
FigS3A <- FeaturePlot(Seurat_Som, features = c("Gad1","Gad2","Cnr1","Cxcl14","Tac3","Htr3a","Vip","Penk","Lamp5",
                                     "Reln","Rgs10","Crh","Cartpt","Pnoc","Sst","Pthlh","Pvalb","Tac1"), 
            pt.size = 0.1, ncol = 3, slot = "scale.data")

#14--------------------------------------------------------Dendrites Clustering--------------------------------------------------------
#15-----------------Transfering Soma Types to Dendrites-----------------
Seurat_Den <- subset(Seurat_Cleaned, Comparment=="Dendrites")
Som_Metadata <- Seurat_Som@meta.data
Den_Metada <- Seurat_Den@meta.data

Som_Den_Intersect <- merge(Den_Metada[,c("Neuron", "SampleName")],
                           Som_Metadata[, c("Neuron", "Specific_type")], all.x = TRUE)
m <- match((Seurat_Den@meta.data$SampleName), Som_Den_Intersect[["SampleName"]])
stopifnot(all(!is.na(m)))
Som_Den_Intersect <- Som_Den_Intersect[m,]
levels <- levels(Som_Den_Intersect$Specific_type)
levels[length(levels)+1]<-"Orphan"
Som_Den_Intersect$Specific_type <-factor(Som_Den_Intersect$Specific_type, levels = levels)
Som_Den_Intersect$Specific_type[is.na(Som_Den_Intersect$Specific_type)] <- "Orphan"
Seurat_Den <- AddMetaData(object = Seurat_Den, metadata = Som_Den_Intersect$Specific_type, col.name = 'Specific_type')

Idents(object = Seurat_Den) <- "Specific_type"
table(Idents(Seurat_Den))

#16--------------------------------------------------------Removing Missed Glial-Contamination & Soma-Contaminated Dendrites---------------------------------------------------------------------
VlnPlot(Seurat_Den,features = c("Dbi","Fabp7","Gfap","S100b","Slc1a2","Vim"), slot = "counts", log = FALSE) # removing dendrites with some remnant glial contaminants
Seurat_Den <- subset(Seurat_Den, slot = 'counts', 
                     subset = Vim <25 & S100b <5 & Fabp7 <5 & Dbi <200 & Slc1a3 < 20 & Gfap < 75 & Slc1a2 < 20)
VlnPlot(Seurat_Den,features = c("Meg3"), slot = "counts", log = FALSE)
Seurat_Den <- subset(Seurat_Den, slot = 'counts', subset = Meg3 < 90)
VlnPlot(Seurat_Den, features = c("nFeature_RNA"), cols = c("#1B9E77","#E6AB02","#386CB0","coral2","coral4","gray"), log = FALSE) + 
  NoLegend() + theme(axis.text.y = element_text(size = 16))
Seurat_Den <- subset(Seurat_Den, subset = nFeature_RNA > 80 & nFeature_RNA < 2000)
table(Idents(Seurat_Den))

#17--------------------------------------------------------Reclustering samples---------------------------------------------------------------------
Seurat_ReCluster <- merge(x = Seurat_Som, y = Seurat_Den)
Seurat_ReCluster <- SCTransform(Seurat_ReCluster, variable.features.n = 3000)
Seurat_ReCluster <- RunPCA(Seurat_ReCluster)
Seurat_ReCluster <- RunUMAP(Seurat_ReCluster, dims = 1:30)
Seurat_ReCluster <- FindNeighbors(Seurat_ReCluster, dims = 1:2, reduction = "umap", verbose = FALSE)
Seurat_ReCluster <- FindClusters(Seurat_ReCluster, verbose = FALSE)

Seurat_ReCluster.markers <- FindAllMarkers(Seurat_ReCluster, min.pct = 0.5, logfc.threshold = 0.5, only.pos = TRUE, test.use = "LR")
top10 <- Seurat_ReCluster.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(Seurat_ReCluster, features = top10$gene, label = TRUE, angle = 45, hjust = 0, size = 4) + NoLegend()

new.cluster.ids <- c("GABAergic","GABAergic","Dendrites","Glutamatergic","Dendrites","GABAergic","GABAergic","Glutamatergic")
names(new.cluster.ids) <- levels(Seurat_ReCluster)
Seurat_ReCluster <- RenameIdents(Seurat_ReCluster, new.cluster.ids)
Seurat_ReCluster <- AddMetaData(object = Seurat_ReCluster, metadata = Seurat_ReCluster@active.ident, col.name = 'Recluster_type')

p1 <- DimPlot(Seurat_ReCluster, reduction = "umap", pt.size = 2, group.by = "Comparment", cols = c("#7570B3","#66A61E"), 
              label = FALSE, label.size = 8) + NoLegend()
FigS2A <- DimPlot(Seurat_ReCluster, reduction = "umap", pt.size = 2, group.by = "Specific_type", 
              cols = c("#1B9E77","#E6AB02","#386CB0","coral2","coral4","gray"), shape.by = "Comparment",
              label= FALSE, label.size = 4) + NoLegend()
plot_grid(p1, FigS2A)

Fig1C <- DimPlot(Seurat_ReCluster, reduction = "umap", pt.size = 3, group.by = "Recluster_type", 
        cols = c("deepskyblue","#7570B3","coral3"), label= FALSE, label.size = 4)+ NoLegend()

Seurat_ReCluster$Type.Comparment <- paste(Seurat_ReCluster$Specific_type, Seurat_ReCluster$Comparment, sep = ".")
Idents(object = Seurat_ReCluster) <- "Type.Comparment"
table(Idents(Seurat_ReCluster))

#18--------------------------------------------------------Downsampling Somata---------------------------------------------------------------------
Seurat_ReCluster_Som <- subset(Seurat_ReCluster, Comparment=="Somata")
Seurat_ReCluster_Den <- subset(Seurat_ReCluster, Comparment=="Dendrites")
Seurat_ReCluster_Som.counts = as.matrix(x = GetAssayData(object = Seurat_ReCluster_Som, assay = "RNA", slot = "counts"))
Seurat_ReCluster_SomDown1Counts <- SampleUMI(data = Seurat_ReCluster_Som.counts, max.umi = 7000)
Seurat_ReCluster_SomDown2Counts <- SampleUMI(data = Seurat_ReCluster_Som.counts, max.umi = 1000)
Seurat_ReCluster_SomDown3Counts <- SampleUMI(data = Seurat_ReCluster_Som.counts, max.umi = 4000)
Seurat_ReCluster_SomDown4Counts <- SampleUMI(data = Seurat_ReCluster_Som.counts, max.umi = 500)
Seurat_ReCluster_SomDown1 <- CreateSeuratObject(Seurat_ReCluster_SomDown1Counts, project = "DownMean")
Seurat_ReCluster_SomDown2 <- CreateSeuratObject(Seurat_ReCluster_SomDown2Counts, project = "DownMedian")
Seurat_ReCluster_SomDown3 <- CreateSeuratObject(Seurat_ReCluster_SomDown3Counts, project = "DownRandom")
Seurat_ReCluster_SomDown4 <- CreateSeuratObject(Seurat_ReCluster_SomDown4Counts, project = "DownRandom")

new.cluster.ids <- c("Somata","Somata","Somata","Somata")
names(new.cluster.ids) <- levels(Seurat_ReCluster_SomDown1)
Seurat_ReCluster_SomDown1 <- RenameIdents(Seurat_ReCluster_SomDown1, new.cluster.ids)
Seurat_ReCluster_SomDown2 <- RenameIdents(Seurat_ReCluster_SomDown2, new.cluster.ids)
Seurat_ReCluster_SomDown3 <- RenameIdents(Seurat_ReCluster_SomDown3, new.cluster.ids)
Seurat_ReCluster_SomDown4 <- RenameIdents(Seurat_ReCluster_SomDown4, new.cluster.ids)
Seurat_ReCluster_SomDown1 <- AddMetaData(object = Seurat_ReCluster_SomDown1, metadata = Seurat_ReCluster_SomDown1@active.ident, col.name = 'Comparment')
Seurat_ReCluster_SomDown2 <- AddMetaData(object = Seurat_ReCluster_SomDown2, metadata = Seurat_ReCluster_SomDown2@active.ident, col.name = 'Comparment')
Seurat_ReCluster_SomDown3 <- AddMetaData(object = Seurat_ReCluster_SomDown3, metadata = Seurat_ReCluster_SomDown3@active.ident, col.name = 'Comparment')
Seurat_ReCluster_SomDown4 <- AddMetaData(object = Seurat_ReCluster_SomDown4, metadata = Seurat_ReCluster_SomDown3@active.ident, col.name = 'Comparment')

VlnPlot(Seurat_ReCluster_SomDown1, features = c("nCount_RNA","nFeature_RNA"), log = TRUE, pt.size = 2)
VlnPlot(Seurat_ReCluster_SomDown2, features = c("nCount_RNA","nFeature_RNA"), log = TRUE, pt.size = 2)
VlnPlot(Seurat_ReCluster_SomDown3, features = c("nCount_RNA","nFeature_RNA"), log = TRUE, pt.size = 2)
VlnPlot(Seurat_ReCluster_SomDown4, features = c("nCount_RNA","nFeature_RNA"), log = TRUE, pt.size = 2)
VlnPlot(Seurat_ReCluster_Den, features = c("nCount_RNA","nFeature_RNA"), log = TRUE, pt.size = 2)

#19a------------------------Cluster DownMean 1------------------------
Seurat_ReCluster_Down1 <- merge(x = Seurat_ReCluster_SomDown1, y = Seurat_ReCluster_Den)
Seurat_ReCluster_Down1 <- SCTransform(Seurat_ReCluster_Down1, variable.features.n = 3000)
Seurat_ReCluster_Down1 <- RunPCA(Seurat_ReCluster_Down1)
Seurat_ReCluster_Down1 <- RunUMAP(Seurat_ReCluster_Down1, dims = 1:30)
Seurat_ReCluster_Down1 <- FindNeighbors(Seurat_ReCluster_Down1, dims = 1:30, verbose = FALSE)
Seurat_ReCluster_Down1 <- FindClusters(Seurat_ReCluster_Down1, verbose = FALSE)
p1 <- DimPlot(Seurat_ReCluster_Down1, reduction = "umap", pt.size = 2, group.by = "Comparment", cols = c("#7570B3","#66A61E"), label = FALSE, label.size = 8) + NoLegend()
FigS2B_TL <- DimPlot(Seurat_ReCluster_Down1, reduction = "umap", pt.size = 2, label= TRUE, label.size = 4) + NoLegend()
plot_grid(p1, FigS2B_TL)

#19b------------------------Cluster DownMean 2------------------------
Seurat_ReCluster_Down2 <- merge(x = Seurat_ReCluster_SomDown2, y = Seurat_ReCluster_Den)
Seurat_ReCluster_Down2 <- SCTransform(Seurat_ReCluster_Down2, variable.features.n = 3000)
Seurat_ReCluster_Down2 <- RunPCA(Seurat_ReCluster_Down2)
Seurat_ReCluster_Down2 <- RunUMAP(Seurat_ReCluster_Down2, dims = 1:30)
Seurat_ReCluster_Down2 <- FindNeighbors(Seurat_ReCluster_Down2, dims = 1:30, verbose = FALSE)
Seurat_ReCluster_Down2 <- FindClusters(Seurat_ReCluster_Down2, verbose = FALSE)
p1 <- DimPlot(Seurat_ReCluster_Down2, reduction = "umap", pt.size = 2, group.by = "Comparment", cols = c("#7570B3","#66A61E"), label = FALSE, label.size = 8) + NoLegend()
FigS2B_TR <- DimPlot(Seurat_ReCluster_Down2, reduction = "umap", pt.size = 2, label= TRUE, label.size = 4) + NoLegend()
plot_grid(p1, FigS2B_TR)

#19c------------------------Cluster DownMean 3------------------------
Seurat_ReCluster_Down3 <- merge(x = Seurat_ReCluster_SomDown3, y = Seurat_ReCluster_Den)
Seurat_ReCluster_Down3 <- SCTransform(Seurat_ReCluster_Down3, variable.features.n = 3000)
Seurat_ReCluster_Down3 <- RunPCA(Seurat_ReCluster_Down3)
Seurat_ReCluster_Down3 <- RunUMAP(Seurat_ReCluster_Down3, dims = 1:30)
Seurat_ReCluster_Down3 <- FindNeighbors(Seurat_ReCluster_Down3, dims = 1:30, verbose = FALSE)
Seurat_ReCluster_Down3 <- FindClusters(Seurat_ReCluster_Down3, verbose = FALSE)
p1 <- DimPlot(Seurat_ReCluster_Down3, reduction = "umap", pt.size = 2, group.by = "Comparment", cols = c("#7570B3","#66A61E"), label = FALSE, label.size = 8) + NoLegend()
FigS2B_BL <- DimPlot(Seurat_ReCluster_Down3, reduction = "umap", pt.size = 2, label= TRUE, label.size = 4) + NoLegend()
plot_grid(p1, FigS2B_BL)

#19d------------------------Cluster DownMean 4------------------------
Seurat_ReCluster_Down4 <- merge(x = Seurat_ReCluster_SomDown4, y = Seurat_ReCluster_Den)
Seurat_ReCluster_Down4 <- SCTransform(Seurat_ReCluster_Down4, variable.features.n = 3000)
Seurat_ReCluster_Down4 <- RunPCA(Seurat_ReCluster_Down4)
Seurat_ReCluster_Down4 <- RunUMAP(Seurat_ReCluster_Down4, dims = 1:30)
Seurat_ReCluster_Down4 <- FindNeighbors(Seurat_ReCluster_Down4, dims = 1:30, verbose = FALSE)
Seurat_ReCluster_Down4 <- FindClusters(Seurat_ReCluster_Down4, verbose = FALSE)
p1 <- DimPlot(Seurat_ReCluster_Down4, reduction = "umap", pt.size = 2, group.by = "Comparment", cols = c("#7570B3","#66A61E"), label = FALSE, label.size = 8) + NoLegend()
FigS2B_BR <- DimPlot(Seurat_ReCluster_Down4, reduction = "umap", pt.size = 2, label= TRUE, label.size = 4) + NoLegend()
plot_grid(p1, FigS2B_BR)

#20------------------------------Expression-Detection Correlation-----------------------------
Som_DotPlot_SCT <- DotPlot(object = Seurat_Som, assay = "RNA", features =rownames(Seurat_Som), group.by = "Comparment")
Som_SCT <-Som_DotPlot_SCT$data
Som_SCT$FractionExp <- Som_SCT$pct.exp/100
Den_DotPlot_SCT <- DotPlot(object = Seurat_Den, assay = "RNA", features =rownames(Seurat_Den), group.by = "Comparment")
Den_SCT <-Den_DotPlot_SCT$data
Den_SCT$FractionExp <- Den_SCT$pct.exp/100
Som_Den_SCT <- merge(x = Som_SCT, y = Den_SCT, by = "features.plot")
rownames(Som_Den_SCT) <- Som_Den_SCT$features.plot

Som_Counts <- GetAssayData(object = Seurat_Som, assay = "RNA")
Den_Counts <- GetAssayData(object = Seurat_Den, assay = "RNA")

Som_Counts <- as.data.frame(Som_Counts)
Den_Counts <- as.data.frame(Den_Counts)
Cell_Averages <- as.data.frame(rowMeans(Som_Counts))
names(Cell_Averages)[1] <- "Som_Exp"
Cell_Averages$Den_Exp <- rowMeans(Den_Counts)
Cell_Ave_Plot <- melt(Cell_Averages)
Som_Counts <- as.matrix(Som_Counts)
Den_Counts <- as.matrix(Den_Counts)
Cell_Averages$Som_Pct <- (rowSums(Som_Counts !=0))/ncol(Som_Counts)
Cell_Averages$Den_Pct <- (rowSums(Den_Counts !=0))/ncol(Den_Counts)

p1 <- ggplot(Cell_Averages, aes(Som_Exp, Som_Pct )) + geom_point(size=1.2, colour = "gray25") + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), colour="red", alpha =1) + 
  scale_x_log10() + ylim(0,1)
p2 <-ggplot(Cell_Averages, aes(Som_Exp, Den_Pct)) + geom_point(size=1.2, colour = "gray25") + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), colour="red", alpha =1) + 
  scale_x_log10() + ylim(0,1)
p3<-ggplot(Cell_Averages, aes(Den_Exp, Den_Pct)) + geom_point(size=1.2, colour = "gray25") + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), colour="red", alpha =1) + 
  scale_x_log10() + ylim(0,1)
plot_grid(p1,p2,p3, ncol = 3)

#--------------------------------------------------------Integration of Middleton Dataset---------------------------------------------------------------------
#21-------------------------Loading & Formating Middleton Files-------------------------------
Middleton_Count_table<-"Middleton_2019_expression.tsv"
Middleton_Samples_Table <- "Middleton_2019_metadata.tsv"
Middleton_countfile<-read.delim(Middleton_Count_table, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
rownames(Middleton_countfile) <- Middleton_countfile$GeneName
Middleton_countfile <- dplyr::select(Middleton_countfile, -1,-34)
Middleton_metadata <- read.delim(Middleton_Samples_Table, check.names = FALSE, header = TRUE)

Middleton_sce <- CreateSeuratObject(Middleton_countfile, min.cells = 0, min.features = 0)
Middleton_sce <- AddMetaData(object = Middleton_sce, metadata = Middleton_metadata$Compartment, col.name = 'Comparment')
Middleton_sce <- AddMetaData(object = Middleton_sce, metadata = Middleton_metadata$Neuron, col.name = 'Neuron')
Middleton_sce <- AddMetaData(object = Middleton_sce, metadata = Middleton_metadata$Source, col.name = 'Source')
Middleton_sce <- SCTransform(Middleton_sce)
Idents(object = Middleton_sce) <- "Compartment"

#22-------------------------Integrating to Middleton Dataset-------------------------------
Middleton_sce <- AddMetaData(object = Middleton_sce, metadata = Middleton_sce$Source, col.name = 'IntegrationConstrat')
Seurat_ReCluster <- AddMetaData(object = Seurat_ReCluster, metadata = Seurat_ReCluster$Specific_type, col.name = 'IntegrationConstrat')

Middleton_list <-c(Seurat_ReCluster,Middleton_sce)
names(Middleton_list)<-c("Seurat_ReCluster","Middleton_sce")
Middleton_features <- SelectIntegrationFeatures(object.list = Middleton_list, nfeatures = 3000)
Middleton_list <- PrepSCTIntegration(object.list = Middleton_list, anchor.features = Middleton_features)
reference_dataset <- which(names(Middleton_list) == "Middleton_sce")
Middleton_anchors <- FindIntegrationAnchors(object.list = Middleton_list, normalization.method = "SCT", 
                                         anchor.features = Middleton_features, reference = reference_dataset, k.filter = 32)
Middleton_integrated <- IntegrateData(anchorset = Middleton_anchors, normalization.method = "SCT")
Middleton_integrated <- RunPCA(object = Middleton_integrated, verbose = FALSE)
Middleton_integrated <- RunUMAP(object = Middleton_integrated, dims = 1:30)

Idents(object = Middleton_integrated) <- "IntegrationConstrat"
Our_cells <-Cells(Seurat_ReCluster)

FigS2C_L <- DimPlot(Middleton_integrated, group.by = "IntegrationConstrat", cols = c("#1B9E77","#E6AB02","#386CB0","coral2","coral4","transparent","transparent","gray40"),
            shape.by = "Comparment", label = FALSE, pt.size = 2, label.size = 6) + NoLegend()
FigS2C_R <- DimPlot(Middleton_integrated, group.by = "IntegrationConstrat", cols = c("gray","gray","gray","gray","gray","#7570B3","#66A61E","gray"),
            shape.by = "Comparment", label = FALSE, pt.size = 2, label.size = 6) + NoLegend()
plot_grid(FigS2C_L,FigS2C_R)

#--------------------------------------------------------Integration to Tissue Datasets---------------------------------------------------------------------
#-------------------------Zeisel 2015-------------------------------
#23-------------------------Loading & Formating Zeisel Files-------------------------------
Zeisel_Count_table<-"Zeisel_2015_expression.txt"
Zeisel_Samples_Table <- "Zeisel_2015_metadata.txt"
Zeisel_countfile<-read.delim(Zeisel_Count_table, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
Zeisel_metadata <- read.delim(Zeisel_Samples_Table, check.names = FALSE, header = TRUE)
rownames(Zeisel_countfile) = make.names(Zeisel_countfile[,1], unique=TRUE)
names(Zeisel_countfile)[1] <- "Gene_Names"
Zeisel_countfile <- dplyr::select(Zeisel_countfile, -1)
names(Zeisel_metadata)[1] <- "Sample_Name"

#24-------------------------Create Seurat Object for Zeisel-------------------------------
Zeisel_sce <- CreateSeuratObject(Zeisel_countfile, project = "Hippocampus Tissue", min.cells = 0, min.features = 0)
Zeisel_sce <- AddMetaData(object = Zeisel_sce, metadata = Zeisel_metadata$tissue, col.name = 'Tissue')
Zeisel_sce <- subset(Zeisel_sce, Tissue=="ca1hippocampus")
Zeisel_sce$orig.ident<-Zeisel_sce$Tissue

#25-------------------------Cluster Zeisel dataset-------------------------------
Zeisel_sce <- SCTransform(Zeisel_sce, return.only.var.genes = FALSE)
Clustering_Zeisel_sce <- RunPCA(Zeisel_sce)
Clustering_Zeisel_sce <- RunUMAP(Clustering_Zeisel_sce, dims = 1:30)
Clustering_Zeisel_sce <- FindNeighbors(Clustering_Zeisel_sce, dims = 1:30, verbose = FALSE)
Clustering_Zeisel_sce <- FindClusters(Clustering_Zeisel_sce, verbose = FALSE)
DimPlot(Clustering_Zeisel_sce, reduction = "umap", pt.size = 1, label = TRUE, label.size = 8) + NoLegend()

FeaturePlot(Clustering_Zeisel_sce, features = c("Gad1","Tbr1","Spink8","Mbp","Aldoc","Aif1","Cldn5","Acta2","Zmynd10"), 
            pt.size = 0.5, ncol = 3) #Markers used in Zeisel et al. 2015

Clustering_Zeisel_sce.markers <- FindAllMarkers(Clustering_Zeisel_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- Clustering_Zeisel_sce.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Clustering_Zeisel_sce, features = top10$gene, label = TRUE, size=4) + NoLegend()

new.cluster.ids <- c("Glut.3","Glut.2","Oligodendrocytes","Glut.5","Glut.4","Astrocytes",
                     "GABA_Cck+","Endothelial","GABA_Sst+","Microglia","Glut.1","Ependymal")#IMPORTANT: groups and order may need to be updated to obtain same results
names(new.cluster.ids) <- levels(Clustering_Zeisel_sce)
Clustering_Zeisel_sce <- RenameIdents(Clustering_Zeisel_sce, new.cluster.ids)
DimPlot(Clustering_Zeisel_sce, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) + NoLegend()

#-------------------------Harris 2018-------------------------------
#26-------------------------Loading & Formating Harris Files-------------------------------
Harris_Count_table<-"Harris_2018_expression.tsv"
Harris_Samples_Table <- "Harris_2018_cell_metadata.tsv"
Harris_countfile<-read.delim(Harris_Count_table, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
Harris_metadata <- read.delim(Harris_Samples_Table, check.names = FALSE, header = TRUE)

rownames(Harris_countfile) = make.names(Harris_countfile[,1], unique=TRUE)
names(Harris_countfile)[1] <- "Gene_Names"
names(Harris_countfile)[3665] <- "Empty"
Harris_countfile <- dplyr::select(Harris_countfile, -1,-3665)
names(Harris_metadata)[1] <- "Sample_Name"
Harris_sce <- CreateSeuratObject(Harris_countfile, project = "Hippocampus Inhibitory", min.cells = 0, min.features = 0, meta.data = Harris_metadata)
Harris_sce <- AddMetaData(object = Harris_sce, metadata = Harris_sce$orig.ident, col.name = 'Tissue')

#27-------------------------Cluster Harris dataset-------------------------------
Harris_sce <- SCTransform(Harris_sce, return.only.var.genes = FALSE)
Clustering_Harris_sce <- RunPCA(Harris_sce)
Clustering_Harris_sce <- RunUMAP(Clustering_Harris_sce, dims = 1:14)
Clustering_Harris_sce <- FindNeighbors(Clustering_Harris_sce, dims = 1:14, verbose = FALSE)
Clustering_Harris_sce <- FindClusters(Clustering_Harris_sce, verbose = FALSE)
DimPlot(Clustering_Harris_sce, reduction = "umap", label = TRUE, pt.size = 1, label.size = 5) + NoLegend()

FeaturePlot(Clustering_Harris_sce, features = c("Calb1","Calb2","Cck","Chrm2","Cnr1",
                                                "Cxcl14","Erbb4","Gabrd","Grm1","Lamp5",
                                                "Lhx6","Ndnf","Nos1","Npy","Ntng1",
                                                "Pcp4","Penk","Pvalb","Reln","Satb1",
                                                "Slc17a8","Sst","Tac1","Tac2","Vip"), 
            pt.size = 0.1, ncol = 5, slot = "scale.data") #Markers used in Harris et al. 2018

Clustering_Harris_sce.markers <- FindAllMarkers(Clustering_Harris_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Clustering_Harris_sce.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- Clustering_Harris_sce.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(Clustering_Harris_sce, features = top10$gene, label = TRUE) + NoLegend()

new.cluster.ids <- c("HS1","Ivy-MGE-NGF","RLMb-Cck-IN1","IS3","IS1","OBi","HS2","NB-Cck-IN","CGE-NGF1","Axo-axonic",
                     "RRH/Cck-In","RLMb-Cck-IN2","Basket/Bis","OLM","Trilaminar","RLMb-Cck-IN3","Sst-Projection","CGE-NGF2","Glial")
names(new.cluster.ids) <- levels(Clustering_Harris_sce)
Clustering_Harris_sce <- RenameIdents(Clustering_Harris_sce, new.cluster.ids)
DimPlot(Clustering_Harris_sce, reduction = "umap", label = TRUE, pt.size = 1, label.size = 5) + NoLegend()

#28-------------------------Integrating to Single Cell Datasets to Somata-------------------------------
Zeisel_sce <- AddMetaData(object = Zeisel_sce, metadata = Clustering_Zeisel_sce@active.ident, col.name = 'Old_Cell_type')
Zeisel_sce <- AddMetaData(object = Zeisel_sce, metadata = Clustering_Zeisel_sce$orig.ident, col.name = 'New_Cell_type')
Zeisel_sce <- AddMetaData(object = Zeisel_sce, metadata = Clustering_Zeisel_sce$Tissue, col.name = 'Comparment')
Zeisel_sce$Comparment <- "Somata"

Harris_sce <- AddMetaData(object = Harris_sce, metadata = Clustering_Harris_sce@active.ident, col.name = 'Old_Cell_type')
Harris_sce <- AddMetaData(object = Harris_sce, metadata = Clustering_Harris_sce$orig.ident, col.name = 'New_Cell_type')
Harris_sce <- AddMetaData(object = Harris_sce, metadata = Clustering_Harris_sce$Tissue, col.name = 'Comparment')
Harris_sce$Comparment <- "Somata"

Seurat_Som <- AddMetaData(object = Seurat_Som, metadata = Seurat_Som$Specific_type, col.name = 'New_Cell_type')
Seurat_Som <- AddMetaData(object = Seurat_Som, metadata = Seurat_Som$Source, col.name = 'Old_Cell_type')
Seurat_Som <- AddMetaData(object = Seurat_Som, metadata = Seurat_Som$Comparment, col.name = 'Tissue')

Som_list <-c(Zeisel_sce, Harris_sce, Seurat_Som)
names(Som_list)<-c("Zeisel_sce","Harris_sce","Seurat_Som")
Som_features <- SelectIntegrationFeatures(object.list = Som_list, nfeatures = 3000)
Som_list <- PrepSCTIntegration(object.list = Som_list, anchor.features = Som_features)
reference_dataset <- which(names(Som_list) == "Zeisel_sce")
Som_anchors <- FindIntegrationAnchors(object.list = Som_list, normalization.method = "SCT", 
                                     anchor.features = Som_features, reference = reference_dataset)
Som_integrated <- IntegrateData(anchorset = Som_anchors, normalization.method = "SCT")
Som_integrated <- RunPCA(object = Som_integrated, verbose = FALSE)
Som_integrated <- RunUMAP(object = Som_integrated, dims = 1:30)

Idents(object = Som_integrated) <- "Old_Cell_type"
Zeisel_cells <-Cells(Zeisel_sce)
Harris_cells <-Cells(Harris_sce)
Som_Cells <- Cells(Seurat_Som)

FigS3B_BL<-DimPlot(Som_integrated, cells =Zeisel_cells, label = TRUE, pt.size = 0.75, label.size = 4) + NoLegend()
FigS3B_BR<-DimPlot(Som_integrated, cells =Harris_cells, label = TRUE, pt.size = 0.75, label.size = 4) + NoLegend()
plot_grid(FigS3B_BL,FigS3B_BR)

FigS3B_TR <- DimPlot(Som_integrated, group.by = "Tissue", label = FALSE, pt.size = 0.75, label.size = 4)

Fig1F_L<-DimPlot(Som_integrated, label = TRUE, pt.size = 0.75, label.size = 4) + NoLegend()
Fig1F_R<-DimPlot(Som_integrated, group.by = "New_Cell_type", cols = c("gray","#1B9E77","#E6AB02","#386CB0","coral2","coral4","gray","black"),
            pt.size = 0.75, label = FALSE, label.size = 4) + NoLegend()
plot_grid(Fig1F_L,Fig1F_R)

#--------------------------------------------------------Cell-Type Specific Effects---------------------------------------------------------------------
#29--------------------------------------------------------Dendrites: Removing Orphan Dendrites---------------------------------------------------------------------
Idents(object = Seurat_Den) <- "Specific_type"
new.cluster.ids <- c("GABAergic_1","GABAergic_2","GABAergic_3","Glutamatergic","Glutamatergic","Orphan")
names(new.cluster.ids) <- levels(Seurat_Den)
Seurat_Den <- RenameIdents(Seurat_Den, new.cluster.ids)
Seurat_Den <- AddMetaData(object = Seurat_Den, metadata = Seurat_Den@active.ident, col.name = 'SemiGeneral_type')
table(Idents(Seurat_Den))

Idents(object = Seurat_Den) <- "SemiGeneral_type"
Seurat_Den_All <- Seurat_Den
is.Orphan <- grepl("^Orphan", Seurat_Den@active.ident)
Seurat_Den <- Seurat_Den[,!is.Orphan]

#30--------------------------------------------------------Average Expression according to cell type---------------------------------------------------------------------
Idents(object = Seurat_Som) <- "Specific_type"
new.cluster.ids <- c("GABAergic_1","GABAergic_2","GABAergic_3","Glutamatergic","Glutamatergic")
names(new.cluster.ids) <- levels(Seurat_Som)
Seurat_Som <- RenameIdents(Seurat_Som, new.cluster.ids)
Seurat_Som <- AddMetaData(object = Seurat_Som, metadata = Seurat_Som@active.ident, col.name = 'SemiGeneral_type')

VlnPlot(Seurat_Som, features = c("nCount_RNA","nFeature_RNA"), cols = c("#1B9E77","#E6AB02","#386CB0","coral3"), log = TRUE, pt.size = 2)
VlnPlot(Seurat_Den, features = c("nCount_RNA","nFeature_RNA"), cols = c("#1B9E77","#E6AB02","#386CB0","coral3"), log = TRUE, pt.size = 2)

#31--------------------------------------------------------Somata: Calculating_Metrics---------------------------------------------------------------------
avg.Seurat_Som <- exp(log1p(AverageExpression(Seurat_Som, verbose = FALSE)$SCT))
avg.Seurat_Som$gene <- rownames(avg.Seurat_Som)
avg.Seurat_Som$mean <- rowMeans(subset(avg.Seurat_Som, select = c(Glutamatergic, GABAergic_1, GABAergic_2, GABAergic_3)))
avg.Seurat_Som$meanGABA <- rowMeans(subset(avg.Seurat_Som, select = c(GABAergic_1, GABAergic_2, GABAergic_3)))
avg.Seurat_Som$meanGABA2.3 <- rowMeans(subset(avg.Seurat_Som, select = c(GABAergic_2, GABAergic_3)))
avg.Seurat_Som$lfc_Glu <- log(exp(avg.Seurat_Som$Glutamatergic)/exp(avg.Seurat_Som$meanGABA))
avg.Seurat_Som$lfc_GA1 <- log(exp(avg.Seurat_Som$GABAergic_1)/exp(avg.Seurat_Som$meanGABA2.3))
avg.Seurat_Som$lfc_GA2.3 <- log(exp(avg.Seurat_Som$GABAergic_3)/exp(avg.Seurat_Som$GABAergic_2))
write.table(avg.Seurat_Som, file="Somata_Cell_Types_AverageExpression2.tsv", sep="\t", quote=FALSE, col.names=TRUE)

#32--------------------------------------------------------Dendrites: Calculating_Metrics---------------------------------------------------------------------
avg.Seurat_Den <- exp(log1p(AverageExpression(Seurat_Den, verbose = FALSE)$SCT))
avg.Seurat_Den$gene <- rownames(avg.Seurat_Den)
avg.Seurat_Den$mean <- rowMeans(subset(avg.Seurat_Den, select = c(Glutamatergic, GABAergic_1, GABAergic_2, GABAergic_3)))
avg.Seurat_Den$meanGABA <- rowMeans(subset(avg.Seurat_Den, select = c(GABAergic_1, GABAergic_2, GABAergic_3)))
avg.Seurat_Den$meanGABA2.3 <- rowMeans(subset(avg.Seurat_Den, select = c(GABAergic_2, GABAergic_3)))
avg.Seurat_Den$lfc_Glu <- log(exp(avg.Seurat_Den$Glutamatergic)/exp(avg.Seurat_Den$meanGABA))
avg.Seurat_Den$lfc_GA1 <- log(exp(avg.Seurat_Den$GABAergic_1)/exp(avg.Seurat_Den$meanGABA2.3))
avg.Seurat_Den$lfc_GA2.3 <- log(exp(avg.Seurat_Den$GABAergic_3)/exp(avg.Seurat_Den$GABAergic_2))
write.table(avg.Seurat_Den, file="Dendrite_Cell_Types_AverageExpression.tsv", sep="\t", quote=FALSE, col.names=TRUE)

#33--------------------------------------------------------Genes of Interest---------------------------------------------------------------------
GOI <- c("Adarb1","Atp5f1b","Arpp21","Babam2","Bdnf","Calm1","Cck","Cdk5r1","Cnr1","Cplx2","Cort","Cox6a1","Cox8a","Fth1",
         "Gabra1","Gad2","Gap43","Gria2","Hpcal4","Lgi2","Malat1","Map1a","Meg3","Nap1l5","Nnat","Nsg1","Nxph1",
         "Psmb2","Psmd6","Rpl41","Rps16","Satb1","Serf2","Slc32a1","Sst")

#34--------------------------------------------------------Somata Cell type: Rank Expression---------------------------------------------------------------------
avg.Glu_Som <- avg.Seurat_Som[order(avg.Seurat_Som$Glutamatergic),]
avg.Glu_Som$Glu_Rank <- seq.int(nrow(avg.Glu_Som))
Fig2B_Glu <-ggplot(avg.Glu_Som, aes(Glu_Rank, Glutamatergic)) + geom_point(size=2, colour = "coral3") +
  geom_point(data = subset(avg.Glu_Som, gene %in% GOI), size=2, shape = 21, fill= "coral3", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.Glu_Som, gene %in% GOI), aes(label = gene), size = 3, max.overlaps = Inf,
                  nudge_x = -500, nudge_y = 0.5)

avg.GABA1_Som <- avg.Seurat_Som[order(avg.Seurat_Som$GABAergic_1),]
avg.GABA1_Som$GABA1_Rank <- seq.int(nrow(avg.GABA1_Som))
Fig2B_GABA1 <-ggplot(avg.GABA1_Som, aes(GABA1_Rank, GABAergic_1)) + geom_point(size=2, colour = "#1B9E77") +
  geom_point(data = subset(avg.GABA1_Som, gene %in% GOI), size=2, shape = 21, fill= "#1B9E77", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.GABA1_Som, gene %in% GOI), aes(label = gene), size = 3, max.overlaps = Inf,
                  nudge_x = -500, nudge_y = 0.5)

avg.GABA2_Som <- avg.Seurat_Som[order(avg.Seurat_Som$GABAergic_2),]
avg.GABA2_Som$GABA2_Rank <- seq.int(nrow(avg.GABA2_Som))
Fig2B_GABA2 <-ggplot(avg.GABA2_Som, aes(GABA2_Rank, GABAergic_2)) + geom_point(size=2, colour = "#E6AB02") +
  geom_point(data = subset(avg.GABA2_Som, gene %in% GOI), size=2, shape = 21, fill= "#E6AB02", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.GABA2_Som, gene %in% GOI), aes(label = gene), size = 3, max.overlaps = Inf,
                  nudge_x = -500, nudge_y = 0.5)

avg.GABA3_Som <- avg.Seurat_Som[order(avg.Seurat_Som$GABAergic_3),]
avg.GABA3_Som$GABA3_Rank <- seq.int(nrow(avg.GABA3_Som))
Fig2B_GABA3 <-ggplot(avg.GABA3_Som, aes(GABA3_Rank, GABAergic_3)) + geom_point(size=2, colour = "#386CB0") +
  geom_point(data = subset(avg.GABA3_Som, gene %in% GOI), size=2, shape = 21, fill= "#386CB0", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.GABA3_Som, gene %in% GOI), aes(label = gene), size = 3, max.overlaps = Inf,
                  nudge_x = -500, nudge_y = 0.5)

#35--------------------------------------------------------Dendrites Cell type: Rank Expression---------------------------------------------------------------------
Den_Glu_GOI <- c("Adarb1","Arpp21","Atp5f1b","Bdnf","Calm1","Cox6a1","Cox8a","Fth1","Gad2","Gria2","Hpcal4",
                 "Malat1","Map1a","Meg3","Nap1l5","Nnat","Nsg1","Psmb2","Rpl41","Serf2")
avg.Glu_Den <- avg.Seurat_Den[order(avg.Seurat_Den$Glutamatergic),]
avg.Glu_Den <- subset(avg.Glu_Den, Glutamatergic > 1.9)
avg.Glu_Den$Glu_Rank <- seq.int(nrow(avg.Glu_Den))
Fig2A_Glu<-ggplot(avg.Glu_Den, aes(Glu_Rank, Glutamatergic)) + geom_point(size=2, colour = "coral3") +
  geom_point(data = subset(avg.Glu_Den, gene %in% Den_Glu_GOI), size=2, shape = 21, fill= "coral3", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.Glu_Den, gene %in% Den_Glu_GOI), aes(label = gene), size = 3, max.overlaps = Inf,
                  nudge_x = -500, nudge_y = 0.5)

Den_GABA1_GOI <- c("Atp5f1b","Babam2","Bdnf","Calm1","Cck","Cnr1","Cort","Cox6a1","Cox8a","Cplx2","Fth1","Gad2",
                   "Gria2","Hpcal4","Lgi2","Malat1","Map1a","Meg3","Nap1l5","Nnat","Nsg1","Psmb2","Rpl41","Satb1","Serf2","Sst")
avg.GABA1_Den <- avg.Seurat_Den[order(avg.Seurat_Den$GABAergic_1),]
avg.GABA1_Den <- subset(avg.GABA1_Den, GABAergic_1 > 1.9)
avg.GABA1_Den$GABA1_Rank <- seq.int(nrow(avg.GABA1_Den))
Fig2A_GABA1 <-ggplot(avg.GABA1_Den, aes(GABA1_Rank, GABAergic_1)) + geom_point(size=2, colour = "#1B9E77") +
  geom_point(data = subset(avg.GABA1_Den, gene %in% Den_GABA1_GOI), size=2, shape = 21, fill= "#1B9E77", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.GABA1_Den, gene %in% Den_GABA1_GOI), aes(label = gene), size = 3, max.overlaps = Inf,
                  nudge_x = -500, nudge_y = 0.5)

Den_GABA2_GOI <- c("Adarb1","Arpp21","Atp5f1b","Babam2","Bdnf","Calm1","Cck","Cnr1","Cort","Cox6a1","Cox8a","Cplx2",
                   "Fth1","Gad2","Gria2","Hpcal4","Lgi2","Malat1","Map1a","Meg3","Nap1l5","Nnat","Nxph1",
                   "Psmb2","Rpl41","Serf2","Sst")
avg.GABA2_Den <- avg.Seurat_Den[order(avg.Seurat_Den$GABAergic_2),]
avg.GABA2_Den <- subset(avg.GABA2_Den, GABAergic_2 > 1.9)
avg.GABA2_Den$GABA2_Rank <- seq.int(nrow(avg.GABA2_Den))
Fig2A_GABA2 <-ggplot(avg.GABA2_Den, aes(GABA2_Rank, GABAergic_2)) + geom_point(size=2, colour = "#E6AB02") +
  geom_point(data = subset(avg.GABA2_Den, gene %in% Den_GABA2_GOI), size=2, shape = 21, fill= "#E6AB02", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.GABA2_Den, gene %in% Den_GABA2_GOI), aes(label = gene), size = 3, max.overlaps = Inf,
                  nudge_x = -500, nudge_y = 0.5)

Den_GABA3_GOI <- c("Adarb1","Atp5f1b","Babam2","Bdnf","Calm1","Cck","Cnr1","Cort","Cox6a1","Cox8a","Cplx2","Fth1","Gad2",
                   "Gria2","Hpcal4","Lgi2","Malat1","Map1a","Meg3","Nap1l5","Nnat","Nsg1","Nxph1","Rpl41","Psmb2","Sst")
avg.GABA3_Den <- avg.Seurat_Den[order(avg.Seurat_Den$GABAergic_3),]
avg.GABA3_Den <- subset(avg.GABA3_Den, GABAergic_3 > 1.9)
avg.GABA3_Den$GABA3_Rank <- seq.int(nrow(avg.GABA3_Den))
Fig2A_GABA3 <-ggplot(avg.GABA3_Den, aes(GABA3_Rank, GABAergic_3)) + geom_point(size=2, colour = "#386CB0") +
  geom_point(data = subset(avg.GABA3_Den, gene %in% Den_GABA3_GOI), size=2, shape = 21, fill= "#386CB0", colour = "black", stroke = 1) +
  theme_classic() + theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
    axis.title.y = element_blank(),) + scale_y_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(avg.GABA3_Den, gene %in% Den_GABA3_GOI), aes(label = gene), size = 4, max.overlaps = Inf, 
                  nudge_x = -500, nudge_y = 0.5)

plot_grid(Fig2A_Glu, Fig2A_GABA1, Fig2A_GABA2, Fig2A_GABA3, ncol = 4)
plot_grid(Fig2B_Glu, Fig2B_GABA1, Fig2B_GABA2, Fig2B_GABA3, ncol = 4)

#--------------------------------------------------------Dendrites: Cell-type Differential Expression---------------------------------------------------------------------
#36--------------------------------------------------------Dendrites Clustering---------------------------------------------------------------------
Seurat_Den <- SCTransform(Seurat_Den, variable.features.n = 3000, return.only.var.genes = FALSE)
Seurat_Den <- RunPCA(Seurat_Den)
Seurat_Den <- RunUMAP(Seurat_Den, dims = 1:30)
Seurat_Den <- FindNeighbors(Seurat_Den, dims = 1:15, verbose = FALSE)
Seurat_Den <- FindClusters(Seurat_Den, verbose = FALSE)
p1 <- DimPlot(Seurat_Den, reduction = "umap", pt.size = 4, cols = c("darkorange3","bisque4","red")) + NoLegend()
p2 <- DimPlot(Seurat_Den, reduction = "umap", pt.size = 4, group.by = "Specific_type",cols = c("#1B9E77","#E6AB02","#386CB0","coral2","coral4","gray")) + NoLegend()
plot_grid(p1,p2)

Idents(object = Seurat_Den) <- "SemiGeneral_type"
avg.Seurat_Den <- exp(log1p(AverageExpression(Seurat_Den, verbose = FALSE)$SCT))
avg.Seurat_Den$gene <- rownames(avg.Seurat_Den)
avg.Seurat_Den$mean <- rowMeans(subset(avg.Seurat_Den, select = c(Glutamatergic, GABAergic_1, GABAergic_2, GABAergic_3)))
avg.Seurat_Den$meanGABA <- rowMeans(subset(avg.Seurat_Den, select = c(GABAergic_1, GABAergic_2, GABAergic_3)))
avg.Seurat_Den$meanGABA2.3 <- rowMeans(subset(avg.Seurat_Den, select = c(GABAergic_2, GABAergic_3)))

#37a--------------------------------------------------------Differential Expression: Glutamatergic---------------------------------------------------------------------
Idents(object = Seurat_Den) <- "SemiGeneral_type"
Glu_Den.markers <- FindMarkers(Seurat_Den, ident.1 = "Glutamatergic", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')

RealGluDen <- as.data.frame(Glu_Den.markers$p_val)
names(RealGluDen)[1] <- "Real"
RealGluDen<-melt(RealGluDen)

Seurat_Den$MockGlu <- sample(1:2, size=97, prob = c(11,86), replace=TRUE)
Idents(object = Seurat_Den) <- "MockGlu"
MockGluDen_p <- c()
for(i in 1:1000) {
set.seed(Sys.time())
Idents(object = Seurat_Den) <- "MockGlu"
Seurat_Den$MockGlu <- sample(1:2, size=97, prob = c(11,86), replace=TRUE)
MockGlu.markers <- FindMarkers(Seurat_Den, ident.1 = "1",ident.2 = "2", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')
Temp <- as.data.frame(MockGlu.markers$p_val)
names(Temp)[1] <- i
Temp<-melt(Temp)
MockGluDen_p <- rbind(MockGluDen_p, Temp)
}

GluDen_p <- rbind(MockGluDen_p,RealGluDen)
GluDen_p$colors <- "Mock"
GluDen_p[which(GluDen_p$variable == "Real"),"colors"] <- "Formal"
MeanGluDen <- MockGluDen_p
MeanGluDen$colors <- "Average"
MeanGluDen$variable <- "Average"
FinalGluDen <- rbind(GluDen_p,MeanGluDen)
MeanGluDen <- MeanGluDen[order(MeanGluDen$value),]
MeanGluDen$pct <- 1/nrow(MeanGluDen)
MeanGluDen$Cummulative <- cumsum(MeanGluDen$pct)

MeanGluDen <- MeanGluDen[order(MeanGluDen$value),]
RealGluDen$colors <- "Real"
RealGluDen$pct <- 1/nrow(RealGluDen)
RealGluDen$Cummulative <- cumsum(RealGluDen$pct)
Real_vs_Mean <- rbind(MeanGluDen,RealGluDen)

FigS4C_L <- ggplot(FinalGluDen, aes(x=value, group = variable, color = colors, size = colors)) + geom_density() + 
  scale_colour_manual(values = c("black","blue","gray70")) + scale_size_manual(values = c(1.5,1.5,0.5)) +
  geom_vline(xintercept = 0.02, color = "firebrick", size = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.0001,1)) + annotation_logticks()

v0 <- 0.05
f1 <- approxfun(MeanGluDen$value, MeanGluDen$Cummulative)
Glu_Intersect <- optimize(function(t0) abs(f1(t0) - v0), interval = range(MeanGluDen$value))
Glu_Intersect$minimum

FigS4D_L <- ggplot(Real_vs_Mean, aes(value, Cummulative, color = colors)) + geom_line(size = 2) + 
  scale_colour_manual(values = c("black","blue")) +
  geom_hline(yintercept = v0, color = "firebrick", size = 1, linetype = "dashed") +
  geom_vline(xintercept = Glu_Intersect$minimum, color = "firebrick", size = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.0001,1)) + annotation_logticks()

Glu_Den.markers$gene <- rownames(Glu_Den.markers)
Glu_Den.markers <- merge(Glu_Den.markers, avg.Seurat_Den, by = "gene", all.y = TRUE)
Glu_Den.markers$colors <- "Unchanged"
Glu_Den.markers[which(Glu_Den.markers$p_val < Glu_Intersect$minimum & Glu_Den.markers$avg_logFC > 0),"colors"] <- "Glu_Up"
Glu_Den.markers[which(Glu_Den.markers$p_val < Glu_Intersect$minimum & Glu_Den.markers$avg_logFC < 0),"colors"] <- "Glu_Down"
Glu_Den.markers <- dplyr::arrange(Glu_Den.markers, desc(colors))
rownames(Glu_Den.markers) <- Glu_Den.markers$gene
Glu_GOI <- c("Adarb1","Arpp21","Bdnf","Cox5b","Hpcal4","Gad2","Lgi2","Nnat","Schip1","Slc32a1","Uqcrq")
Glu_Den_plot <- ggplot(Glu_Den.markers, aes(meanGABA, Glutamatergic)) + 
  geom_point(aes(fill=as.factor(colors )), size=3, shape = 21, stroke = 0) +
  geom_point(data = subset(Glu_Den.markers, gene %in% Glu_GOI), size=2.5, shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = c("cadetblue3","coral3","gray82")) + theme_classic() + 
  theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
                          axis.title.y = element_blank(),) + 
  scale_y_log10() + scale_x_log10() + annotation_logticks() +
  geom_text_repel(data = subset(Glu_Den.markers, gene %in% Glu_GOI), aes(label = gene), size = 3.5, max.overlaps = Inf, 
                  nudge_x = 0.1, nudge_y = 0.1)
table(Glu_Den.markers$colors)
cor(Glu_Den.markers$meanGABA,Glu_Den.markers$Glutamatergic)

#37b--------------------------------------------------------Differential Expression: GABAergic_1 vs GABAergic 2 & 3---------------------------------------------------------------------
Idents(object = Seurat_Den) <- "SemiGeneral_type"
is.Glutamatergic <- grepl("^Glutamatergic", Seurat_Den@active.ident)
Seurat_Den_GABA <- Seurat_Den[,!is.Glutamatergic]
GABA1_Den.markers <- FindMarkers(Seurat_Den_GABA, ident.1 = "GABAergic_1", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')

RealGABA1Den <- as.data.frame(GABA1_Den.markers$p_val)
names(RealGABA1Den)[1] <- "Real"
RealGABA1Den<-melt(RealGABA1Den)
Seurat_Den_GABA$MockGABA1 <- sample(1:2, size=86, prob = c(38,48), replace=TRUE)
Idents(object = Seurat_Den_GABA) <- "MockGABA1"
MockGABA1Den_p <- c()
for(i in 1:1000) {
  set.seed(Sys.time())
  Idents(object = Seurat_Den_GABA) <- "MockGABA1"
  Seurat_Den_GABA$MockGABA1 <- sample(1:2, size=86, prob = c(38,48), replace=TRUE)
  MockGABA1.markers <- FindMarkers(Seurat_Den_GABA, ident.1 = "1",ident.2 = "2", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')
  Temp <- as.data.frame(MockGABA1.markers$p_val)
  names(Temp)[1] <- i
  Temp<-melt(Temp)
  MockGABA1Den_p <- rbind(MockGABA1Den_p, Temp)
}

GABA1Den_p <- rbind(MockGABA1Den_p,RealGABA1Den)
GABA1Den_p$colors <- "Mock"
GABA1Den_p[which(GABA1Den_p$variable == "Real"),"colors"] <- "Formal"
MeanGABA1Den <- MockGABA1Den_p
MeanGABA1Den$colors <- "Average"
MeanGABA1Den$variable <- "Average"
FinalGABA1Den <- rbind(GABA1Den_p,MeanGABA1Den)
MeanGABA1Den <- MeanGABA1Den[order(MeanGABA1Den$value),]
MeanGABA1Den$pct <- 1/nrow(MeanGABA1Den)
MeanGABA1Den$Cummulative <- cumsum(MeanGABA1Den$pct)

RealGABA1Den$colors <- "Real"
RealGABA1Den$pct <- 1/nrow(RealGABA1Den)
RealGABA1Den$Cummulative <- cumsum(RealGABA1Den$pct)
Real_vs_Mean_GABA1 <- rbind(MeanGABA1Den,RealGABA1Den)
MeanGABA1Den <- MeanGABA1Den[!duplicated(MeanGABA1Den$value),]

FigS4C_M <- ggplot(FinalGABA1Den, aes(x=value, group = variable, color = colors, size = colors)) + geom_density() + 
  scale_colour_manual(values = c("black","blue","gray70")) + scale_size_manual(values = c(1.5,1.5,0.5)) +
  geom_vline(xintercept = 0.02, color = "firebrick", size = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.0001,1)) + annotation_logticks()

v0 <- 0.05
f1 <- approxfun(MeanGABA1Den$value, MeanGABA1Den$Cummulative)
GABA1_Intersect <- optimize(function(t0) abs(f1(t0) - v0), interval = range(MeanGABA1Den$value))
GABA1_Intersect$minimum

FigS4D_M <- ggplot(Real_vs_Mean_GABA1, aes(value, Cummulative, color = colors)) + geom_line(size = 2) + 
  scale_colour_manual(values = c("black","blue")) +
  geom_hline(yintercept = v0, color = "firebrick", size = 1, linetype = "dashed") +
  geom_vline(xintercept = GABA1_Intersect$minimum, color = "firebrick", size = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.0001,1)) + annotation_logticks()

GABA1_Den.markers$gene <- rownames(GABA1_Den.markers)
GABA1_Den.markers <- merge(GABA1_Den.markers, avg.Seurat_Den, by = "gene", all.y = TRUE)
GABA1_Den.markers$colors <- "Unchanged"
GABA1_Den.markers[which(GABA1_Den.markers$p_val < GABA1_Intersect$minimum & GABA1_Den.markers$avg_logFC > 0),"colors"] <- "GABA1_Up"
GABA1_Den.markers[which(GABA1_Den.markers$p_val < GABA1_Intersect$minimum & GABA1_Den.markers$avg_logFC < 0),"colors"] <- "GABA1_Down"
write.table(GABA1_Den.markers, file="GABAergic1_Dendrite_Differential_Expression.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)
GABA1_Den.markers <- dplyr::arrange(GABA1_Den.markers, desc(colors))
rownames(GABA1_Den.markers) <- GABA1_Den.markers$gene
GABA1_GOI <- c("Atp2b2","Babam2","Cck","Ccni","Cnr1","Cort","Gpr88","LOC252890-Cox7c","Satb1","Sst","Tmsb10","Unc119b")
GABA1_Den_plot<-ggplot(GABA1_Den.markers, aes(meanGABA2.3, GABAergic_1)) + 
  geom_point(aes(fill=as.factor(colors )), size=3, shape = 21, stroke = 0) +
  geom_point(data = subset(GABA1_Den.markers, gene %in% GABA1_GOI), size=2.5, shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = c("plum","#1B9E77","gray82")) + theme_classic() + 
  theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(),) + 
  scale_y_log10() + scale_x_log10() + annotation_logticks() +
  geom_text_repel(data = subset(GABA1_Den.markers, gene %in% GABA1_GOI), aes(label = gene), size = 3.5, max.overlaps = Inf, 
                  nudge_x = 0.1, nudge_y = 0.1)

table(GABA1_Den.markers$colors)
cor(GABA1_Den.markers$meanGABA2.3,GABA1_Den.markers$GABAergic_1)

#37c--------------------------------------------------------Differential Expression: GABAergic_2 vs GABAergic 3---------------------------------------------------------------------
Idents(object = Seurat_Den_GABA) <- "SemiGeneral_type"
is.GABA1 <- grepl("^GABAergic_1", Seurat_Den_GABA@active.ident)
Seurat_Den_GABA23 <- Seurat_Den_GABA[,!is.GABA1]
GABA23_Den.markers <- FindMarkers(Seurat_Den_GABA23, ident.1 = "GABAergic_2",ident.2 = "GABAergic_3", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')

RealGABA23Den <- as.data.frame(GABA23_Den.markers$p_val)
names(RealGABA23Den)[1] <- "Real"
RealGABA23Den<-melt(RealGABA23Den)
Seurat_Den_GABA23$MockGABA23 <- sample(1:2, size=48, prob = c(30,18), replace=TRUE)
Idents(object = Seurat_Den_GABA23) <- "MockGABA23"
MockGABA23Den_p <- c()
for(i in 1:1000) {
  set.seed(Sys.time())
  Idents(object = Seurat_Den_GABA23) <- "MockGABA23"
  Seurat_Den_GABA23$MockGABA23 <- sample(1:2, size=48, prob = c(30,18), replace=TRUE)
  MockGABA23.markers <- FindMarkers(Seurat_Den_GABA23, ident.1 = "1",ident.2 = "2", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')
  Temp <- as.data.frame(MockGABA23.markers$p_val)
  names(Temp)[1] <- i
  Temp<-melt(Temp)
  MockGABA23Den_p <- rbind(MockGABA23Den_p, Temp)
}

GABA23Den_p <- rbind(MockGABA23Den_p,RealGABA23Den)
GABA23Den_p$colors <- "Mock"
GABA23Den_p[which(GABA23Den_p$variable == "Real"),"colors"] <- "Formal"
MeanGABA23Den <- MockGABA23Den_p
MeanGABA23Den$colors <- "Average"
MeanGABA23Den$variable <- "Average"
FinalGABA23Den <- rbind(GABA23Den_p,MeanGABA23Den)
MeanGABA23Den <- MeanGABA23Den[order(MeanGABA23Den$value),]
MeanGABA23Den$pct <- 1/nrow(MeanGABA23Den)
MeanGABA23Den$Cummulative <- cumsum(MeanGABA23Den$pct)
RealGABA23Den$colors <- "Real"
RealGABA23Den$pct <- 1/nrow(RealGABA23Den)
RealGABA23Den$Cummulative <- cumsum(RealGABA23Den$pct)
Real_vs_Mean_GABA23 <- rbind(MeanGABA23Den,RealGABA23Den)
MeanGABA23Den <- MeanGABA23Den[!duplicated(MeanGABA23Den$value),]

FigS4C_R <- ggplot(FinalGABA23Den, aes(x=value, group = variable, color = colors, size = colors)) + geom_density() + 
  scale_colour_manual(values = c("black","blue","gray70")) + scale_size_manual(values = c(1.5,1.5,0.5)) +
  geom_vline(xintercept = 0.02, color = "firebrick", size = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.0001,1)) + annotation_logticks()

v0 <- 0.05
f1 <- approxfun(MeanGABA23Den$value, MeanGABA23Den$Cummulative)
GABA23_Intersect <- optimize(function(t0) abs(f1(t0) - v0), interval = range(MeanGABA23Den$value))
GABA23_Intersect$minimum

FigS4D_R <- ggplot(Real_vs_Mean, aes(value, Cummulative, color = colors)) + geom_line(size = 2) + 
  scale_colour_manual(values = c("black","blue")) +
  geom_hline(yintercept = v0, color = "firebrick", size = 1, linetype = "dashed") +
  geom_vline(xintercept = GABA23_Intersect$minimum, color = "firebrick", size = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.0001,1)) + annotation_logticks()

GABA23_Den.markers$gene <- rownames(GABA23_Den.markers)
GABA23_Den.markers <- merge(GABA23_Den.markers, avg.Seurat_Den, by = "gene", all.y = TRUE)
GABA23_Den.markers$colors <- "Unchanged"
GABA23_Den.markers[which(GABA23_Den.markers$p_val < GABA23_Intersect$minimum & GABA23_Den.markers$avg_logFC > 0),"colors"] <- "GABA2_Up"
GABA23_Den.markers[which(GABA23_Den.markers$p_val < GABA23_Intersect$minimum & GABA23_Den.markers$avg_logFC < 0),"colors"] <- "GABA3_Up"
write.table(GABA23_Den.markers, file="GABAergic2&3_Dendrite_Differential_Expression.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)
GABA23_Den.markers <- dplyr::arrange(GABA23_Den.markers, desc(colors))
rownames(GABA23_Den.markers) <- GABA23_Den.markers$gene
GABA23_GOI <- c("Adarb2-094873","Cck","Cnr1","Cplx2","Cyc1","Fabp5","Gng3","Mapk9","Nxph1","Sox2","Sst","Tubb2b")
GABA23_Den_plot<-ggplot(GABA23_Den.markers, aes(GABAergic_2, GABAergic_3)) + 
  geom_point(aes(fill=as.factor(colors )), size=3, shape = 21, stroke = 0) +
  geom_point(data = subset(GABA23_Den.markers, gene %in% GABA23_GOI), size=2.5, shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#E6AB02","#386CB0","gray82")) + theme_classic() + 
  theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(),) + 
  scale_y_log10() + scale_x_log10() + annotation_logticks() + 
  geom_text_repel(data = subset(GABA23_Den.markers, gene %in% GABA23_GOI), aes(label = gene), size = 3.5, max.overlaps = Inf, 
                  nudge_x = 0.1, nudge_y = 0.1)

table(GABA23_Den.markers$colors)
cor(GABA23_Den.markers$GABAergic_2,GABA23_Den.markers$GABAergic_3)

#38--------------------------------------------------------Somata: Cell-type Differential Expression---------------------------------------------------------------------
#38a--------------------------------------------------------Differential Expression: Glutamatergic---------------------------------------------------------------------
Glu_Som.markers <- FindMarkers(Seurat_Som, ident.1 = "Glutamatergic", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')
Glu_Som.markers$gene <- rownames(Glu_Som.markers)
Glu_Som.markers <- merge(Glu_Som.markers, avg.Seurat_Som, by = "gene", all.y = TRUE)
Glu_Som.markers$colors <- "Unchanged"
Glu_Som.markers[which(Glu_Som.markers$p_val_adj < 0.05 & Glu_Som.markers$avg_logFC > 0),"colors"] <- "Glu_Up"
Glu_Som.markers[which(Glu_Som.markers$p_val_adj < 0.05 & Glu_Som.markers$avg_logFC < 0),"colors"] <- "Glu_Down"
write.table(Glu_Som.markers, file="Glutamatergic_Somata_Differential_Expression.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)
Glu_Som.markers <- dplyr::arrange(Glu_Som.markers, desc(colors))
rownames(Glu_Som.markers) <- Glu_Som.markers$gene
Glu_Som_plot <- ggplot(Glu_Som.markers, aes(meanGABA, Glutamatergic)) + 
  geom_point(aes(fill=as.factor(colors )), size=3, shape = 21, stroke = 0) +
  geom_point(data = subset(Glu_Som.markers, gene %in% Glu_GOI), size=2.5, shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = c("cadetblue3","coral3","gray82")) + theme_classic() + 
  theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(),) + 
  scale_y_log10() + scale_x_log10() + annotation_logticks() +
  geom_text_repel(data = subset(Glu_Som.markers, gene %in% Glu_GOI), aes(label = gene), size = 3.5, max.overlaps = Inf, 
                  nudge_x = 0.1, nudge_y = 0.1)

table(Glu_Som.markers$colors)
cor(Glu_Som.markers$meanGABA,Glu_Som.markers$Glutamatergic)

#38b--------------------------------------------------------Differential Expression: GABAergic_1 vs GABAergic 2 & 3---------------------------------------------------------------------
is.Glutamatergic <- grepl("^Glutamatergic", Seurat_Som@active.ident)
Seurat_Som_GABA <- Seurat_Som[,!is.Glutamatergic]
GABA1_Som.markers <- FindMarkers(Seurat_Som_GABA, ident.1 = "GABAergic_1", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')
GABA1_Som.markers$gene <- rownames(GABA1_Som.markers)
GABA1_Som.markers <- merge(GABA1_Som.markers, avg.Seurat_Som, by = "gene", all.y = TRUE)
GABA1_Som.markers$colors <- "Unchanged"
GABA1_Som.markers[which(GABA1_Som.markers$p_val_adj < 0.05 & GABA1_Som.markers$avg_logFC > 0),"colors"] <- "GABA1_Up"
GABA1_Som.markers[which(GABA1_Som.markers$p_val_adj < 0.05 & GABA1_Som.markers$avg_logFC < 0),"colors"] <- "GABA1_Down"
write.table(GABA1_Som.markers, file="GABAergic1_Somata_Differential_Expression.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)
GABA1_Som.markers <- dplyr::arrange(GABA1_Som.markers, desc(colors))
rownames(GABA1_Som.markers) <- GABA1_Som.markers$gene
GABA1_Som_plot<-ggplot(GABA1_Som.markers, aes(meanGABA2.3, GABAergic_1)) + 
  geom_point(aes(fill=as.factor(colors )), size=3, shape = 21, stroke = 0) +
  geom_point(data = subset(GABA1_Som.markers, gene %in% GABA1_GOI), size=2.5, shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = c("plum","#1B9E77","gray82")) + theme_classic() + 
  theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(),) + 
  scale_y_log10() + scale_x_log10() + annotation_logticks() +
  geom_text_repel(data = subset(GABA1_Som.markers, gene %in% GABA1_GOI), aes(label = gene), size = 3.5, max.overlaps = Inf, 
                  nudge_x = 0.1, nudge_y = 0.1)
table(GABA1_Som.markers$colors)
cor(GABA1_Som.markers$meanGABA2.3,GABA1_Som.markers$GABAergic_1)

#38c--------------------------------------------------------Differential Expression: GABAergic_2 vs GABAergic 3---------------------------------------------------------------------
GABA2.3_Som.markers <- FindMarkers(Seurat_Som_GABA, ident.1 = "GABAergic_2",ident.2 = "GABAergic_3", min.pct = 0.1, logfc.threshold = 0.5, test.use = 'LR')
GABA2.3_Som.markers$gene <- rownames(GABA2.3_Som.markers)
GABA2.3_Som.markers <- merge(GABA2.3_Som.markers, avg.Seurat_Som, by = "gene", all.y = TRUE)
GABA2.3_Som.markers$colors <- "Unchanged"
GABA2.3_Som.markers[which(GABA2.3_Som.markers$p_val_adj < 0.05 & GABA2.3_Som.markers$avg_logFC > 0),"colors"] <- "GABA2_Up"
GABA2.3_Som.markers[which(GABA2.3_Som.markers$p_val_adj < 0.05 & GABA2.3_Som.markers$avg_logFC < 0),"colors"] <- "GABA3_Up"
write.table(GABA2.3_Som.markers, file="GABAergic2&3_Somata_Differential_Expression.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)
GABA2.3_Som.markers <- dplyr::arrange(GABA2.3_Som.markers, desc(colors))
rownames(GABA2.3_Som.markers) <- GABA2.3_Som.markers$gene
GABA2.3_Som_plot<-ggplot(GABA2.3_Som.markers, aes(GABAergic_2, GABAergic_3)) + 
  geom_point(aes(fill=as.factor(colors )), size=3, shape = 21, stroke = 0) +
  geom_point(data = subset(GABA2.3_Som.markers, gene %in% GABA23_GOI), size=2.5, shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#E6AB02","#386CB0","gray82")) + theme_classic() + 
  theme(text = element_text(size=20),legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(),) + 
  scale_y_log10() + scale_x_log10() + annotation_logticks() +
  geom_text_repel(data = subset(GABA2.3_Som.markers, gene %in% GABA23_GOI), aes(label = gene), size = 3.5, max.overlaps = Inf, 
                  nudge_x = 0.1, nudge_y = 0.1)
table(GABA2.3_Som.markers$colors)
cor(GABA2.3_Som.markers$GABAergic_2,GABA2.3_Som.markers$GABAergic_3)

#39--------------------------------------------------------Differential Expression Plots---------------------------------------------------------------------
Fig3A <- plot_grid(Glu_Som_plot, Glu_Den_plot, GABA1_Som_plot, GABA1_Den_plot, GABA2.3_Som_plot, GABA23_Den_plot, ncol = 2)

#40--------------------------------------------------------Plotting Genes---------------------------------------------------------------------
my_levels <- c("Glutamatergic","GABAergic_1","GABAergic_2","GABAergic_3") 
Idents(Seurat_Den) <- factor(Idents(Seurat_Den), levels = my_levels)
Idents(Seurat_Som) <- factor(Idents(Seurat_Som), levels = my_levels)

GOI <- c("Adarb1","Adarb2-094873","Arpp21","Atp2b2","Babam2","Bdnf","Cck","Ccni","Cnr1","Cort","Cox5b",
         "Cplx2","Cyc1","Fabp5","Gad2","Gng3","Gpr88","Hpcal4","Lgi2","LOC252890-Cox7c","Mapk9","Nnat",
         "Nxph1","Satb1","Schip1","Slc32a1","Sox2","Sst","Tmsb10","Tubb2b","Unc119b","Uqcrq")
         
FigS5A_L<- DoHeatmap(Seurat_Som, features = GOI, label = FALSE, slot = "counts", group.bar = TRUE,
          group.colors = c("coral3","#1B9E77","#E6AB02","#386CB0"), disp.max = 30, disp.min = 6) + NoLegend() + 
  theme(axis.text.y = element_text(size = 8)) + scale_fill_viridis(option = "cividis")

FigS5A_R<- DoHeatmap(Seurat_Den, features = GOI, label = FALSE, slot = "counts", group.bar = TRUE, 
          group.colors = c("coral3","#1B9E77","#E6AB02","#386CB0"), disp.min = 1) + NoLegend() + 
  theme(axis.text.y = element_text(size = 8)) + scale_fill_viridis(option = "cividis")

DefaultAssay(Seurat_Den) = "SCT"
Fig3B_L <- VlnPlot(Seurat_Som,features = c("Hpcal4","Cort","Cnr1","Adarb1","Cyc1","Babam2","Bdnf","Satb1","Adarb2-094878"), log = TRUE, 
        slot = 'counts', ncol = 3, cols = c("coral3","#1B9E77","#E6AB02","#386CB0"), pt.size = 0)
Fig3B_R <- VlnPlot(Seurat_Den,features = c("Hpcal4","Cort","Cnr1","Adarb1","Cyc1","Babam2","Bdnf","Satb1","Adarb2-094878"), log = TRUE,
        slot = 'counts', ncol = 3, cols = c("coral3","#1B9E77","#E6AB02","#386CB0"), pt.size = 0)

#--------------------------------------------------------Somata vs Dendrites---------------------------------------------------------------------
#41------------------------------Info Transfer-----------------------------
#Transfer identity of remaining dendrites to Seurat_All object
Den_Metadata <- Seurat_Den@meta.data
Seurat_Cell <- Seurat_All
Cell_Metada <- Seurat_Cell@meta.data
Den_Cell_Intersect <- merge(Cell_Metada[,c("Neuron", "SampleName")],
                            Den_Metadata[, c("Neuron", "SemiGeneral_type")], all.x = TRUE)

m <- match((Seurat_Cell@meta.data$SampleName), Den_Cell_Intersect[["SampleName"]])
stopifnot(all(!is.na(m)))
Den_Cell_Intersect <- Den_Cell_Intersect[m,]

Seurat_Cell <- AddMetaData(object = Seurat_Cell, metadata = Den_Cell_Intersect$SemiGeneral_type, col.name = 'SemiGeneral_type')
Idents(object = Seurat_Cell) <- "SemiGeneral_type"

NA_Values <-is.na(Seurat_Cell$SemiGeneral_type)
Seurat_Cell <- Seurat_Cell[,!NA_Values]

Seurat_Cell$Cell.Comparment <- paste(Seurat_Cell$Neuron, Seurat_Cell$Comparment, sep = ".")
Seurat_Cell$Type.Comparment <- paste(Seurat_Cell$SemiGeneral_type, Seurat_Cell$Comparment, sep = ".")
Seurat_Cell$Type.Cell.Comparment <- paste(Seurat_Cell$SemiGeneral_type, Seurat_Cell$Neuron, Seurat_Cell$Comparment, sep = ".")
Idents(object = Seurat_Cell) <- "Type.Comparment"

Shared_Som <- subset(Seurat_Cell, Comparment=="Somata")
Shared_Den <- subset(Seurat_Cell, Comparment=="Dendrites")

#42------------------------------Dendritic vs Somatic Expression 1-----------------------------
Som_Counts <- GetAssayData(object = Shared_Som, assay = "RNA")
Den_Counts <- GetAssayData(object = Shared_Den, assay = "RNA")

Som_Counts <- as.data.frame(Som_Counts)
Den_Counts <- as.data.frame(Den_Counts)
Cell_Averages <- as.data.frame(rowMeans(Som_Counts))
names(Cell_Averages)[1] <- "Som_Exp"
Cell_Averages$Den_Exp <- rowMeans(Den_Counts)
Cell_Ave_Plot <- melt(Cell_Averages)
Som_Counts <- as.matrix(Som_Counts)
Den_Counts <- as.matrix(Den_Counts)
Cell_Averages$Som_Pct <- (rowSums(Som_Counts !=0))/95
Cell_Averages$Den_Pct <- (rowSums(Den_Counts !=0))/95
Cell_Averages <- subset(Cell_Averages, Den_Exp > 0)

FigS6A <-ggplot(Cell_Ave_Plot, aes(x= value, fill = variable)) + geom_density(alpha = 1) + 
  scale_fill_manual(values = c("#66A61E","#7570B3")) + scale_x_log10() + annotation_logticks()

Fig4A <- ggplot(Cell_Averages, aes(Som_Exp, Den_Pct)) + geom_point(size=5, colour = "gray35") + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), colour="firebrick", alpha =1) #+ xlim(0,100) + ylim(0,1)
Fig4A<-LabelPoints(plot = Fig4A, points = c("Atp1b1","Atp5f1b","Calm1","Ccni","Cox6a1","Cox8a","Fth1",
                                            "Gria2","Malat1","Map1a","Meg3","Nsg1","Psmb2","Psmc5",
                                            "Psmd6","Rian","Rpl41","Rps29","Rpl3", "Rps3"), repel = TRUE)
plot_grid(Fig4A)

#43------------------------------Logistic Regression Model-----------------------------
Fraction_logistic <- glm(Den_Pct ~ Som_Exp, data = Cell_Averages, family = "binomial")
summary(Fraction_logistic)
ll.null <- Fraction_logistic$null.deviance/-2
ll.proposed <- Fraction_logistic$deviance/-2
(ll.null-ll.proposed)/ll.null

#44------------------------------Dendritic vs Somatic Expression 2-----------------------------
Som_Counts2 <- GetAssayData(object = Shared_Som, assay = "SCT", slot = "counts")
Den_Counts2 <- GetAssayData(object = Shared_Den, assay = "SCT", slot = "counts")

Som_Counts2 <- as.data.frame(Som_Counts2)
Den_Counts2 <- as.data.frame(Den_Counts2)
Som_Counts2$gene <- rownames(Som_Counts2)
Den_Counts2$gene <- rownames(Den_Counts2)
Som_Counts2 <- semi_join(Som_Counts2, Den_Counts2, by = "gene")
Den_Counts2 <- semi_join(Den_Counts2, Som_Counts2, by = "gene")
rownames(Som_Counts2) <- Som_Counts2$gene
rownames(Den_Counts2) <- Den_Counts2$gene
Som_Counts2 <- subset(Som_Counts2, select=-c(gene))
Den_Counts2 <- subset(Den_Counts2, select=-c(gene))

Cell_Averages2 <- as.data.frame(rowMeans(Som_Counts2))
names(Cell_Averages2)[1] <- "Som_Exp"
Cell_Averages2$Den_Exp <- rowMeans(Den_Counts2)
Som_Counts2 <- as.matrix(Som_Counts2)
Den_Counts2 <- as.matrix(Den_Counts2)
Cell_Averages2$Som_Pct <- (rowSums(Som_Counts2 !=0))/94
Cell_Averages2$Den_Pct <- (rowSums(Den_Counts2 !=0))/94

Som_Counts2 <- as.matrix(Som_Counts2)
Den_Counts2 <- as.matrix(Den_Counts2)
Sub_Som_Den <- Den_Counts2/Som_Counts2
Sub_Som_Den[is.nan(Sub_Som_Den)] <- 0
Sub_Som_Den[is.infinite(Sub_Som_Den)] <- 0
Cell_Averages2$Den_Rel <- (rowMeans(Sub_Som_Den))/100
Cell_Averages2$gene <- rownames(Cell_Averages2)
Cell_Averages2 <- subset(Cell_Averages2, Den_Pct > 0.4)
Cell_Averages3 <- subset(Cell_Averages2, Som_Exp < 95)

Fig4B <- ggplot(Cell_Averages3, aes(Som_Exp, Den_Exp)) + geom_point(size=5, colour = "gray35") + 
  geom_smooth(method = "glm", colour="firebrick", alpha =1) #+ scale_x_log10() + scale_y_log10()
Fig4B<-LabelPoints(plot = Fig4B, points =c("Atp1b1","Atp5f1b","Calm1","Ccni","Cox6a1","Cox8a","Fth1",
                                           "Gria2","Malat1","Map1a","Meg3","Nsg1","Psmb2","Psmc5",
                                           "Psmd6","Rian","Rpl41","Rps29","Rpl3", "Rps3"), repel = TRUE)
plot_grid(Fig4B)

#45--------------------------------------------------------Somata vs Dendrites Differential Expression---------------------------------------------------------------------
Idents(object = Seurat_Cell) <- "Comparment"
Seurat_Cell <- SCTransform(Seurat_Cell, return.only.var.genes = FALSE)
Dendrites_vs_Somata.markers <- FindMarkers(Seurat_Cell, ident.1 = "Dendrites", ident.2 = "Somata", min.pct = 0.25, logfc.threshold = 0, 
                                           test.use = "poisson", latent.vars = 'Neuron')
Dendrites_vs_Somata.markers$gene <- rownames(Dendrites_vs_Somata.markers)
Dendrites_vs_Somata.markers$colors <- "Constant"
Dendrites_vs_Somata.markers[which(Dendrites_vs_Somata.markers$p_val_adj < 0.05 & Dendrites_vs_Somata.markers$avg_logFC < -0.5),"colors"] <- "Somata"
Dendrites_vs_Somata.markers[which(Dendrites_vs_Somata.markers$p_val_adj < 0.05 & Dendrites_vs_Somata.markers$avg_logFC > 0.5),"colors"] <- "Dendrites"
write.table(Dendrites_vs_Somata.markers, file="Dendrites_vs_Somata_Differential_Expression_NegBinomial.tsv", sep="\t", quote=FALSE, col.names=TRUE)
Dendrites_vs_Somata.markers <- dplyr::arrange(Dendrites_vs_Somata.markers, desc(colors))
rownames(Dendrites_vs_Somata.markers) <- Dendrites_vs_Somata.markers$gene
table(Dendrites_vs_Somata.markers$colors)

avg.Seurat_Cell <- (AverageExpression(Seurat_Cell, verbose = FALSE)$SCT)
avg.Seurat_Cell$gene <- rownames(avg.Seurat_Cell)
Avg_Expression_Comparment <- merge(avg.Seurat_Cell, Dendrites_vs_Somata.markers, by = "gene", all.x=TRUE)
Avg_Expression_Comparment$colors[is.na(Avg_Expression_Comparment$colors)] <- "Constant"
Avg_Expression_Comparment <- subset(Avg_Expression_Comparment, Somata > 0)
Avg_Expression_Comparment <- subset(Avg_Expression_Comparment, Dendrites > 0)
Avg_Expression_Comparment <- Avg_Expression_Comparment[order(Avg_Expression_Comparment$colors),]
rownames(Avg_Expression_Comparment) <- Avg_Expression_Comparment$gene 
table(Avg_Expression_Comparment$colors)
Avg_Expression_Comparment$Cell <- (Avg_Expression_Comparment$Somata + Avg_Expression_Comparment$Dendrites)

Fig4C<-ggplot(Avg_Expression_Comparment, aes(Somata, Dendrites)) + geom_point(aes(color=as.factor(colors )), size=5) + 
  scale_colour_manual(values = c("gray82","#7570B3","#66A61E")) + NoLegend() + xlim(-2.5,6) + ylim(-2.5,6) + scale_y_log10() + scale_x_log10()
Fig4C<-LabelPoints(plot = Fig4C, points = c(c("Atp1b1","Atp5f1b","Calm1","Ccni","Cox6a1","Cox8a","Fth1",
                                                                                  "Gria2","Malat1","Map1a","Meg3","Nsg1","Psmb2","Psmc5",
                                                                                  "Psmd6","Rian","Rpl41","Rps29","Rpl3", "Rps3")), repel = TRUE, ynudge = 0.2)
plot_grid(Fig4C)
cor(Avg_Expression_Comparment$Somata, Avg_Expression_Comparment$Dendrites)

colnames(Avg_Expression_Comparment) <- c("Gene","SomataNormExp","DendritesNormExp","p-value","Avg_logFC","pct.Den","pct.Som",
                                         "adj_p-value","Evaluation","CellNormExp")
Avg_Expression_Comparment <- Avg_Expression_Comparment[order(Avg_Expression_Comparment$`p-value`),]
write.table(Avg_Expression_Comparment, file="Dendrites_and_Somata_Average_Expression.tsv", sep="\t", quote=FALSE, 
            col.names=TRUE, row.names = FALSE)
#46--------------------------------------------------------Dendrites: FISH Validation Dendritic Presence---------------------------------------------------------------------
Comp_FISH_file <- "tableX_DendriteSignalFig6_14Dec2020.txt"
ComptFISH <- read.delim(Comp_FISH_file, check.names = FALSE, header = TRUE)
colnames(ComptFISH) <- c("Gene","Dend.Length","DistalProximal")
ComptFISH$Gene <- factor(ComptFISH$Gene, levels = c("Eef1a1","Cox8a","Psmb2","Psmd6"))

comparisons <- list(c("Eef1a1","Psmb2"),c("Eef1a1","Psmd6"),
                    c("Cox8a","Psmb2"), c("Cox8a","Psmd6"))

FigS6D<-ggplot(ComptFISH, aes(Gene, DistalProximal)) + geom_violin(scale = "width", colour = "black", fill="gray50", lwd = 1) + 
  stat_summary(fun.y=mean, geom="crossbar", size=0.5, color="gray20") + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 3) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  stat_compare_means(label.x = 3, comparisons = comparisons)

FigS6E<-ggplot(ComptFISH, aes(Gene, Dend.Length)) + geom_violin(scale = "width", colour = "black", fill="gray50", lwd = 1) + 
  stat_summary(fun.y=mean, geom="crossbar", size=0.5, color="gray20") + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 3) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  stat_compare_means(label.x = 3, comparisons = comparisons)

plot_grid(FigS6D, FigS6E)

#47------------------------------Extracting Counts from specific Genes-----------------------------
DefaultAssay(Seurat_Cell) = "RNA"
Shared_Den <- subset(Seurat_Cell, Comparment=="Dendrites")
Shared_Som <- subset(Seurat_Cell, Comparment=="Somata")
Den_data <- GetAssayData(object = Shared_Den, slot = "counts")
Den_data <- as.data.frame(Den_data)
Som_data <- GetAssayData(object = Shared_Som, slot = "counts")
Som_data <- as.data.frame(Som_data)
Den_data[Den_data == 0] <- 0.9
Som_data[Som_data == 0] <- 0.9

#------------------------------Atp5f1b-----------------------------
Atp5f1b_Dendrites<- as.numeric(Den_data[grep("^Atp5f1b$", rownames(Den_data)),])
Atp5f1b_Somata<-as.numeric(Som_data[grep("^Atp5f1b$", rownames(Som_data)),])
Atp5f1b <- data.frame(Somata = Atp5f1b_Somata, Dendrites = Atp5f1b_Dendrites)
Atp5f1b$Neuron <-Shared_Den@meta.data$Neuron
Atp5f1b <- tidyr::gather(Atp5f1b, Cell, value, -Neuron)
Atp5f1b_plot <- ggplot(Atp5f1b, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Fth1-----------------------------
Fth1_Dendrites<- as.numeric(Den_data[grep("^Fth1$", rownames(Den_data)),])
Fth1_Somata<-as.numeric(Som_data[grep("^Fth1$", rownames(Som_data)),])
Fth1 <- data.frame(Somata = Fth1_Somata, Dendrites = Fth1_Dendrites)
Fth1$Neuron <-Shared_Den@meta.data$Neuron
Fth1 <- tidyr::gather(Fth1, Cell, value, -Neuron)
Fth1_plot <- ggplot(Fth1, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Cox8a-----------------------------
Cox8a_Dendrites<- as.numeric(Den_data[grep("^Cox8a$", rownames(Den_data)),])
Cox8a_Somata<-as.numeric(Som_data[grep("^Cox8a$", rownames(Som_data)),])
Cox8a <- data.frame(Somata = Cox8a_Somata, Dendrites = Cox8a_Dendrites)
Cox8a$Neuron <-Shared_Den@meta.data$Neuron
Cox8a <- tidyr::gather(Cox8a, Cell, value, -Neuron)
Cox8a_plot <- ggplot(Cox8a, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Nsg1-----------------------------
Nsg1_Dendrites<- as.numeric(Den_data[grep("^Nsg1$", rownames(Den_data)),])
Nsg1_Somata<-as.numeric(Som_data[grep("^Nsg1$", rownames(Som_data)),])
Nsg1 <- data.frame(Somata = Nsg1_Somata, Dendrites = Nsg1_Dendrites)
Nsg1$Neuron <-Shared_Den@meta.data$Neuron
Nsg1 <- tidyr::gather(Nsg1, Cell, value, -Neuron)
Nsg1_plot <- ggplot(Nsg1, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Psmb2-----------------------------
Psmb2_Dendrites<- as.numeric(Den_data[grep("^Psmb2$", rownames(Den_data)),])
Psmb2_Somata<-as.numeric(Som_data[grep("^Psmb2$", rownames(Som_data)),])
Psmb2 <- data.frame(Somata = Psmb2_Somata, Dendrites = Psmb2_Dendrites)
Psmb2$Neuron <-Shared_Den@meta.data$Neuron
Psmb2 <- tidyr::gather(Psmb2, Cell, value, -Neuron)
Psmb2_plot <- ggplot(Psmb2, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

plot_grid(Atp5f1b_plot,Fth1_plot,Cox8a_plot,Nsg1_plot,Psmb2_plot, ncol = 3)

#------------------------------Additional Genes-----------------------------
#------------------------------Map1a-----------------------------
Map1a_Dendrites<- as.numeric(Den_data[grep("^Map1a$", rownames(Den_data)),])
Map1a_Somata<-as.numeric(Som_data[grep("^Map1a$", rownames(Som_data)),])
Map1a <- data.frame(Somata = Map1a_Somata, Dendrites = Map1a_Dendrites)
Map1a$Neuron <-Shared_Den@meta.data$Neuron
Map1a <- tidyr::gather(Map1a, Cell, value, -Neuron)
Map1a_plot <- ggplot(Map1a, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Rpl41-----------------------------
Rpl41_Dendrites<- as.numeric(Den_data[grep("^Rpl41$", rownames(Den_data)),])
Rpl41_Somata<-as.numeric(Som_data[grep("^Rpl41$", rownames(Som_data)),])
Rpl41 <- data.frame(Somata = Rpl41_Somata, Dendrites = Rpl41_Dendrites)
Rpl41$Neuron <-Shared_Den@meta.data$Neuron
Rpl41 <- tidyr::gather(Rpl41, Cell, value, -Neuron)
Rpl41_plot <- ggplot(Rpl41, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Cox6a1-----------------------------
Cox6a1_Dendrites<- as.numeric(Den_data[grep("^Cox6a1$", rownames(Den_data)),])
Cox6a1_Somata<-as.numeric(Som_data[grep("^Cox6a1$", rownames(Som_data)),])
Cox6a1 <- data.frame(Somata = Cox6a1_Somata, Dendrites = Cox6a1_Dendrites)
Cox6a1$Neuron <-Shared_Den@meta.data$Neuron
Cox6a1 <- tidyr::gather(Cox6a1, Cell, value, -Neuron)
Cox6a1_plot <- ggplot(Cox6a1, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Tmsb4x-----------------------------
Tmsb4x_Dendrites<- as.numeric(Den_data[grep("^Tmsb4x$", rownames(Den_data)),])
Tmsb4x_Somata<-as.numeric(Som_data[grep("^Tmsb4x$", rownames(Som_data)),])
Tmsb4x <- data.frame(Somata = Tmsb4x_Somata, Dendrites = Tmsb4x_Dendrites)
Tmsb4x$Neuron <-Shared_Den@meta.data$Neuron
Tmsb4x <- tidyr::gather(Tmsb4x, Cell, value, -Neuron)
Tmsb4x_plot <- ggplot(Tmsb4x, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Calm3-----------------------------
Calm3_Dendrites<- as.numeric(Den_data[grep("^Calm3$", rownames(Den_data)),])
Calm3_Somata<-as.numeric(Som_data[grep("^Calm3$", rownames(Som_data)),])
Calm3 <- data.frame(Somata = Calm3_Somata, Dendrites = Calm3_Dendrites)
Calm3$Neuron <-Shared_Den@meta.data$Neuron
Calm3 <- tidyr::gather(Calm3, Cell, value, -Neuron)
Calm3_plot <- ggplot(Calm3, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Snhg11-----------------------------
Snhg11_Dendrites<- as.numeric(Den_data[grep("^Snhg11$", rownames(Den_data)),])
Snhg11_Somata<-as.numeric(Som_data[grep("^Snhg11$", rownames(Som_data)),])
Snhg11 <- data.frame(Somata = Snhg11_Somata, Dendrites = Snhg11_Dendrites)
Snhg11$Neuron <-Shared_Den@meta.data$Neuron
Snhg11 <- tidyr::gather(Snhg11, Cell, value, -Neuron)
Snhg11_plot <- ggplot(Snhg11, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Meg3-----------------------------
Meg3_Dendrites<- as.numeric(Den_data[grep("^Meg3$", rownames(Den_data)),])
Meg3_Somata<-as.numeric(Som_data[grep("^Meg3$", rownames(Som_data)),])
Meg3 <- data.frame(Somata = Meg3_Somata, Dendrites = Meg3_Dendrites)
Meg3$Neuron <-Shared_Den@meta.data$Neuron
Meg3 <- tidyr::gather(Meg3, Cell, value, -Neuron)
Meg3_plot <- ggplot(Meg3, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Atp1b1-----------------------------
Atp1b1_Dendrites<- as.numeric(Den_data[grep("^Atp1b1$", rownames(Den_data)),])
Atp1b1_Somata<-as.numeric(Som_data[grep("^Atp1b1$", rownames(Som_data)),])
Atp1b1 <- data.frame(Somata = Atp1b1_Somata, Dendrites = Atp1b1_Dendrites)
Atp1b1$Neuron <-Shared_Den@meta.data$Neuron
Atp1b1 <- tidyr::gather(Atp1b1, Cell, value, -Neuron)
Atp1b1_plot <- ggplot(Atp1b1, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Gria2-----------------------------
Gria2_Dendrites<- as.numeric(Den_data[grep("^Gria2$", rownames(Den_data)),])
Gria2_Somata<-as.numeric(Som_data[grep("^Gria2$", rownames(Som_data)),])
Gria2 <- data.frame(Somata = Gria2_Somata, Dendrites = Gria2_Dendrites)
Gria2$Neuron <-Shared_Den@meta.data$Neuron
Gria2 <- tidyr::gather(Gria2, Cell, value, -Neuron)
Gria2_plot <- ggplot(Gria2, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Rian-----------------------------
Rian_Dendrites<- as.numeric(Den_data[grep("^Rian$", rownames(Den_data)),])
Rian_Somata<-as.numeric(Som_data[grep("^Rian$", rownames(Som_data)),])
Rian <- data.frame(Somata = Rian_Somata, Dendrites = Rian_Dendrites)
Rian$Neuron <-Shared_Den@meta.data$Neuron
Rian <- tidyr::gather(Rian, Cell, value, -Neuron)
Rian_plot <- ggplot(Rian, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------First Plots-----------------------------
plot_grid(Map1a_plot,Rpl41_plot,Cox6a1_plot,Tmsb4x_plot,Calm3_plot,
          Snhg11_plot,Meg3_plot,Atp1b1_plot,Gria2_plot,Rian_plot, ncol = 5)

#------------------------------Ribosomal-Proteasome-----------------------------
#------------------------------Rps29-----------------------------
Rps29_Dendrites<- as.numeric(Den_data[grep("^Rps29$", rownames(Den_data)),])
Rps29_Somata<-as.numeric(Som_data[grep("^Rps29$", rownames(Som_data)),])
Rps29 <- data.frame(Somata = Rps29_Somata, Dendrites = Rps29_Dendrites)
Rps29$Neuron <-Shared_Den@meta.data$Neuron
Rps29 <- tidyr::gather(Rps29, Cell, value, -Neuron)
Rps29_plot <- ggplot(Rps29, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Rps27-----------------------------
Rps27_Dendrites<- as.numeric(Den_data[grep("^Rps27$", rownames(Den_data)),])
Rps27_Somata<-as.numeric(Som_data[grep("^Rps27$", rownames(Som_data)),])
Rps27 <- data.frame(Somata = Rps27_Somata, Dendrites = Rps27_Dendrites)
Rps27$Neuron <-Shared_Den@meta.data$Neuron
Rps27 <- tidyr::gather(Rps27, Cell, value, -Neuron)
Rps27_plot <- ggplot(Rps27, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Rps16-----------------------------
Rps16_Dendrites<- as.numeric(Den_data[grep("^Rps16$", rownames(Den_data)),])
Rps16_Somata<-as.numeric(Som_data[grep("^Rps16$", rownames(Som_data)),])
Rps16 <- data.frame(Somata = Rps16_Somata, Dendrites = Rps16_Dendrites)
Rps16$Neuron <-Shared_Den@meta.data$Neuron
Rps16 <- tidyr::gather(Rps16, Cell, value, -Neuron)
Rps16_plot <- ggplot(Rps16, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Rps26-----------------------------
Rps26_Dendrites<- as.numeric(Den_data[grep("^Rps26$", rownames(Den_data)),])
Rps26_Somata<-as.numeric(Som_data[grep("^Rps26$", rownames(Som_data)),])
Rps26 <- data.frame(Somata = Rps26_Somata, Dendrites = Rps26_Dendrites)
Rps26$Neuron <-Shared_Den@meta.data$Neuron
Rps26 <- tidyr::gather(Rps26, Cell, value, -Neuron)
Rps26_plot <- ggplot(Rps26, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Rpl38-----------------------------
Rpl38_Dendrites<- as.numeric(Den_data[grep("^Rpl38$", rownames(Den_data)),])
Rpl38_Somata<-as.numeric(Som_data[grep("^Rpl38$", rownames(Som_data)),])
Rpl38 <- data.frame(Somata = Rpl38_Somata, Dendrites = Rpl38_Dendrites)
Rpl38$Neuron <-Shared_Den@meta.data$Neuron
Rpl38 <- tidyr::gather(Rpl38, Cell, value, -Neuron)
Rpl38_plot <- ggplot(Rpl38, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Psmc5-----------------------------
Psmc5_Dendrites<- as.numeric(Den_data[grep("^Psmc5$", rownames(Den_data)),])
Psmc5_Somata<-as.numeric(Som_data[grep("^Psmc5$", rownames(Som_data)),])
Psmc5 <- data.frame(Somata = Psmc5_Somata, Dendrites = Psmc5_Dendrites)
Psmc5$Neuron <-Shared_Den@meta.data$Neuron
Psmc5 <- tidyr::gather(Psmc5, Cell, value, -Neuron)
Psmc5_plot <- ggplot(Psmc5, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Psmc3-----------------------------
Psmc3_Dendrites<- as.numeric(Den_data[grep("^Psmc3$", rownames(Den_data)),])
Psmc3_Somata<-as.numeric(Som_data[grep("^Psmc3$", rownames(Som_data)),])
Psmc3 <- data.frame(Somata = Psmc3_Somata, Dendrites = Psmc3_Dendrites)
Psmc3$Neuron <-Shared_Den@meta.data$Neuron
Psmc3 <- tidyr::gather(Psmc3, Cell, value, -Neuron)
Psmc3_plot <- ggplot(Psmc3, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Psmd6-----------------------------
Psmd6_Dendrites<- as.numeric(Den_data[grep("^Psmd6$", rownames(Den_data)),])
Psmd6_Somata<-as.numeric(Som_data[grep("^Psmd6$", rownames(Som_data)),])
Psmd6 <- data.frame(Somata = Psmd6_Somata, Dendrites = Psmd6_Dendrites)
Psmd6$Neuron <-Shared_Den@meta.data$Neuron
Psmd6 <- tidyr::gather(Psmd6, Cell, value, -Neuron)
Psmd6_plot <- ggplot(Psmd6, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Psma6-----------------------------
Psma6_Dendrites<- as.numeric(Den_data[grep("^Psma6$", rownames(Den_data)),])
Psma6_Somata<-as.numeric(Som_data[grep("^Psma6$", rownames(Som_data)),])
Psma6 <- data.frame(Somata = Psma6_Somata, Dendrites = Psma6_Dendrites)
Psma6$Neuron <-Shared_Den@meta.data$Neuron
Psma6 <- tidyr::gather(Psma6, Cell, value, -Neuron)
Psma6_plot <- ggplot(Psma6, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Psmd3-----------------------------
Psmd3_Dendrites<- as.numeric(Den_data[grep("^Psmd3$", rownames(Den_data)),])
Psmd3_Somata<-as.numeric(Som_data[grep("^Psmd3$", rownames(Som_data)),])
Psmd3 <- data.frame(Somata = Psmd3_Somata, Dendrites = Psmd3_Dendrites)
Psmd3$Neuron <-Shared_Den@meta.data$Neuron
Psmd3 <- tidyr::gather(Psmd3, Cell, value, -Neuron)
Psmd3_plot <- ggplot(Psmd3, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Second Plots-----------------------------
plot_grid(Rps29_plot,Rps27_plot,Rps16_plot,Rps26_plot,Rpl38_plot,
          Psmc5_plot,Psmc3_plot,Psmd6_plot,Psma6_plot,Psmd3_plot, ncol = 5)

#------------------------------Others-----------------------------
#------------------------------Cplx1-----------------------------
Cplx1_Dendrites<- as.numeric(Den_data[grep("^Cplx1$", rownames(Den_data)),])
Cplx1_Somata<-as.numeric(Som_data[grep("^Cplx1$", rownames(Som_data)),])
Cplx1 <- data.frame(Somata = Cplx1_Somata, Dendrites = Cplx1_Dendrites)
Cplx1$Neuron <-Shared_Den@meta.data$Neuron
Cplx1 <- tidyr::gather(Cplx1, Cell, value, -Neuron)
Cplx1_plot <- ggplot(Cplx1, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Uchl1-----------------------------
Uchl1_Dendrites<- as.numeric(Den_data[grep("^Uchl1$", rownames(Den_data)),])
Uchl1_Somata<-as.numeric(Som_data[grep("^Uchl1$", rownames(Som_data)),])
Uchl1 <- data.frame(Somata = Uchl1_Somata, Dendrites = Uchl1_Dendrites)
Uchl1$Neuron <-Shared_Den@meta.data$Neuron
Uchl1 <- tidyr::gather(Uchl1, Cell, value, -Neuron)
Uchl1_plot <- ggplot(Uchl1, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Gpx4-----------------------------
Gpx4_Dendrites<- as.numeric(Den_data[grep("^Gpx4$", rownames(Den_data)),])
Gpx4_Somata<-as.numeric(Som_data[grep("^Gpx4$", rownames(Som_data)),])
Gpx4 <- data.frame(Somata = Gpx4_Somata, Dendrites = Gpx4_Dendrites)
Gpx4$Neuron <-Shared_Den@meta.data$Neuron
Gpx4 <- tidyr::gather(Gpx4, Cell, value, -Neuron)
Gpx4_plot <- ggplot(Gpx4, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))
#------------------------------Calm1-----------------------------
Calm1_Dendrites<- as.numeric(Den_data[grep("^Calm1$", rownames(Den_data)),])
Calm1_Somata<-as.numeric(Som_data[grep("^Calm1$", rownames(Som_data)),])
Calm1 <- data.frame(Somata = Calm1_Somata, Dendrites = Calm1_Dendrites)
Calm1$Neuron <-Shared_Den@meta.data$Neuron
Calm1 <- tidyr::gather(Calm1, Cell, value, -Neuron)
Calm1_plot <- ggplot(Calm1, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Map1b-----------------------------
Map1b_Dendrites<- as.numeric(Den_data[grep("^Map1b$", rownames(Den_data)),])
Map1b_Somata<-as.numeric(Som_data[grep("^Map1b$", rownames(Som_data)),])
Map1b <- data.frame(Somata = Map1b_Somata, Dendrites = Map1b_Dendrites)
Map1b$Neuron <-Shared_Den@meta.data$Neuron
Map1b <- tidyr::gather(Map1b, Cell, value, -Neuron)
Map1b_plot <- ggplot(Map1b, aes(Cell, value)) + geom_violin(scale = "width", aes(fill = Cell), colour = "black", lwd = 1) + geom_point(size = 0) +
  geom_line(aes(group = Neuron), colour = "gray20", size=0.5) + scale_x_discrete(limits = c('Somata', 'Dendrites')) + 
  scale_fill_manual(values = c("#7570B3","#66A61E")) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limits = c(0.9,1000))

#------------------------------Final Plots-----------------------------
plot_grid(Cplx1_plot,Uchl1_plot,Gpx4_plot,Calm1_plot,Map1b_plot,Map1b_plot, ncol = 5)
