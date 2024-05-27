library(Matrix)
library(readr)
library(edgeR)
library(tidyverse)
library(Seurat)
library(dplyr)
library(patchwork)
library(EnhancedVolcano)



#Reads in data from directory...
spatial.data.mast <- Read10X(data.dir = getwd())

#Reads in csv that labels each cell class
cell_class <- read.csv('D:/Cone/KOR RNA-seq/Mouse Torpor Neurons/cell_class.csv')

#rownames needs to be exact same as seurat for this to work
cell_class <- as.data.frame(cell_class[,-1], row.names = cell_class[,1])
colnames(cell_class) <- c('cell.class')


#create searat object for analysis... seems to take long time...
scrna.mast <- CreateSeuratObject(counts = spatial.data.mast, meta.data = cell_class)


# regular workflow
scrna.mast <- NormalizeData(object = scrna.mast)
scrna.mast <- FindVariableFeatures(object = scrna.mast)
scrna.mast <- ScaleData(object = scrna.mast)
scrna.mast <- RunPCA(object = scrna.mast)
scrna.mast <- FindNeighbors(object = scrna.mast)
scrna.mast <- FindClusters(object = scrna.mast)
scrna.mast <- RunTSNE(object = scrna.mast)

# View metadata data frame, stored in object@meta.data
view_alldata <- scrna.mast[[]]

#Creates subset of just Oprk1+ cells from raw data set... very fast...
# only 1022 cells have at least 1 molecule of oprk1
subset_kor <-subset(x = scrna.mast, subset = Oprk1 > 0)

# View metadata data frame, stored in object@meta.data
view_subsetkor <- subset_kor[[]]



# all positive oprk1 cells - also controls for subset_kor
poscells_oprk1 <- WhichCells(scrna.mast, expression = Oprk1 > 0)


# all positive oprm1 cells - also controls for subset_kor
poscells_oprm1 <- WhichCells(scrna.mast, expression = Oprm1 > 0)


# all positive oprk1 cells - also controls for subset_kor
poscells_oprd1 <- WhichCells(scrna.mast, expression = Oprd1 > 0)




# all oprk1- cells
negcells_oprk1 <- WhichCells(scrna.mast, expression = Oprk1 == 0)

# all oprm1- cells
negcells_oprm1 <- WhichCells(scrna.mast, expression = Oprm1 == 0)

# all oprm1- cells
negcells_oprd1 <- WhichCells(scrna.mast, expression = Oprd1 == 0)


negcells_opioids <- WhichCells(scrna.mast, expression = Oprd1 == 0 & Oprk1 == 0 & Oprm1 == 0)


# ALL OPRK1+ CELLS VS ALL OTHER CELLS
oprk1.markers <- FindMarkers(scrna.mast, 
                             ident.1 = poscells_oprk1, 
                             ident.2 = negcells_oprk1, 
                             test.use = 't', # t-test
                             logfc.threshold = 0)



oprk1.markers.wilcox <- FindMarkers(scrna.mast, 
                                    ident.1 = poscells_oprk1, 
                                    ident.2 = negcells_oprk1, 
                                    test.use = 'wilcox', # wilcoxon rank sum test
                                    logfc.threshold = 0)


oprk1.markers.wilcox$p_val_adj <- p.adjust(oprk1.markers.wilcox$p_val, method = 'bonferroni')

EnhancedVolcano(oprk1.markers, 
                x = 'avg_log2FC', 
                y = 'p_val_adj', 
                lab = rownames(oprk1.markers),
                title = 'Differential Gene Expression in Kappa Opioid Receptor (KOR) Cells',
                subtitle = "Volcano Plot",
                labSize = 8.0)


EnhancedVolcano(oprk1.markers.wilcox, 
                x = 'avg_log2FC', 
                y = 'p_val_adj', 
                lab = rownames(oprk1.markers.wilcox),
                title = 'Differential Gene Expression in Kappa Opioid Receptor (KOR) Cells',
                subtitle = "Volcano Plot",
                drawConnectors = TRUE,
                labSize = 8.0,
                selectLab = c(rownames(head(oprd1.markers)), 'Pdyn', 'Gal', 'Tac1'), FCcutoff = 0.5, col=c('black', 'purple', 'blue', 'red3'))




# ALL OPRM1+ CELLS VS ALL OTHER CELLS
oprm1.markers <- FindMarkers(scrna.mast, 
                             ident.1 = poscells_oprm1, 
                             ident.2 = negcells_oprm1, 
                             test.use = 't', 
                             logfc.threshold = 0)


oprm1.markers.wilcox <- FindMarkers(scrna.mast, 
                                    ident.1 = poscells_oprm1, 
                                    ident.2 = negcells_oprm1, 
                                    test.use = 'wilcox', 
                                    logfc.threshold = 0)

EnhancedVolcano(oprm1.markers, 
                x = 'avg_log2FC', 
                y = 'p_val_adj', 
                lab = rownames(oprm1.markers), 
                drawConnectors = TRUE,
                title = 'Differential Gene Expression in Mu Opioid Receptor (MOR) Cells',
                subtitle = "Volcano Plot",
                labSize = 8.0,
                selectLab = c('Nts','Gal', 'Penk', rownames(head(oprm1.markers, 25))))


EnhancedVolcano(oprm1.markers.wilcox, 
                x = 'avg_log2FC', 
                y = 'p_val_adj', 
                lab = rownames(oprm1.markers.wilcox), 
                drawConnectors = TRUE,
                title = 'Differential Gene Expression in Mu Opioid Receptor (MOR) Cells',
                subtitle = "Volcano Plot",
                labSize = 8.0,
                selectLab = c('Nts','Gal', 'Penk', rownames(head(oprm1.markers, 25))), FCcutoff = 0.5, col=c('black', 'purple', 'blue', 'red3'))



# ALL OPRD1+ CELLS VS ALL OTHER CELLS
oprd1.markers <- FindMarkers(scrna.mast, 
                             ident.1 = poscells_oprd1, 
                             ident.2 = negcells_oprd1, 
                             test.use = 't', 
                             logfc.threshold = 0)


oprd1.markers.wilcox <- FindMarkers(scrna.mast, 
                                    ident.1 = poscells_oprd1, 
                                    ident.2 = negcells_oprd1, 
                                    test.use = 'wilcox', 
                                    logfc.threshold = 0)


EnhancedVolcano(oprd1.markers, 
                x = 'avg_log2FC', 
                y = 'p_val_adj', 
                lab = rownames(oprd1.markers), 
                drawConnectors = TRUE,
                labSize = 8.0,
                title = 'Differential Gene Expression in Delta Opioid Receptor (DOR) Cells',
                subtitle = "Volcano Plot",
                selectLab = c(rownames(head(oprd1.markers)), 'Gal', 'Penk'))


EnhancedVolcano(oprd1.markers.wilcox, 
                x = 'avg_log2FC', 
                y = 'p_val_adj', 
                lab = rownames(oprd1.markers.wilcox), 
                drawConnectors = TRUE,
                labSize = 8.0,
                title = 'Differential Gene Expression in Delta Opioid Receptor (DOR) Cells',
                subtitle = "Volcano Plot",
                pCutoff = 10e-2, 
                selectLab = c(rownames(head(oprd1.markers)), 'Gal', 'Penk'), FCcutoff = 0.5, col=c('black', 'purple', 'blue', 'red3'))



oprk1.markers.wilcox$gene_symbol <- row.names(oprk1.markers.wilcox)


# saves markers to csv file
write.csv(oprd1.markers.wilcox, "oprd1.markers.wilcox.csv")


FeaturePlot(scrna.mast, features = "Pdyn", pt.size = 0.75, cols = c("lightgrey", "red")) & NoLegend()


FeaturePlot(scrna.mast_filtered, features = c("Oprk1", "Nts"), cols = c("blue", "red"), blend = TRUE, pt.size = .75)


FeaturePlot(
  object = scrna.mast,
  features = c(
    "Oprd1", "Oprk1", "Oprm1"
  ),
  ncol = 3, pt.size = 0.75, cols = c("lightgrey", "red")) & NoLegend()


DimPlot(scrna.mast, label = TRUE, pt.size = 0.5, group.by = "cell.class") + ggtitle("")



# Filter out the "unstable" group from your Seurat object
scrna.mast_filtered <- subset(x = scrna.mast, subset = (cell.class == "Unstable" | cell.class == 'Ambiguous'), invert = TRUE)


# Create a DimPlot excluding the "unstable" group
DimPlot(scrna.mast_filtered, label = TRUE, pt.size = 0.75, group.by = "cell.class", label.size = 6) + ggtitle("") + ggtitle("") +
  theme(text = element_text(size = 14, face = "bold")) 






# Use DimPlot to generate a plot and extract the colors
colors_for_plot <- DimPlot(scrna.mast_filtered, label = TRUE, pt.size = 1, group.by = "cell.class", return.plot = TRUE)$label_params$fill

# Get the cell class labels and counts
cell_class_counts <- table(scrna.mast_filtered$cell.class)

cell_class_counts <- as.data.frame(cell_class_counts)


VlnPlot(scrna.mast_filtered, features = c("Oprd1"), group.by = 'cell.class', pt.size = 2) +
  NoLegend() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 5))

VlnPlot(scrna.mast_filtered, features = c("Oprk1"), group.by = 'cell.class', pt.size = 1) +
  NoLegend() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 5))

VlnPlot(scrna.mast_filtered, features = c("Oprm1"), group.by = 'cell.class', pt.size = 1) +
  NoLegend() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 3))


# Creates cell visualization for inhibitory cells and opioid receptor cell markers
feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprk1", "Pdyn"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 

feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprk1", "Tac1"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 

feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprk1", "Gal"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 





feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprd1", "Penk"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 

feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprd1", "Gal"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 



feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprm1", "Penk"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 

feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprm1", "Nts"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 

feat3 <-FeaturePlot(inhibitory_cells, features = c("Oprm1", "Gal"), blend = TRUE, pt.size = 2, cols = c("blue", "red"))

feat3[[3]] + NoLegend() 


feat3[[4]]




inhibitory_cells <- subset(scrna.mast_filtered, subset = cell.class == "Inhibitory")


excitatory_cells <- subset(scrna.mast_filtered, subset = cell.class == "Excitatory")

inhibitory_oprd1 <- WhichCells(inhibitory_cells, expression = Oprd1 > 0)

inhibitory_cells_oprk1 <- WhichCells(inhibitory_cells, expression = Oprk1 > 0)

inhibitory_cells_oprm1 <- WhichCells(inhibitory_cells, expression = Oprm1 > 0)


inhibitory_oprd1_neg <- WhichCells(inhibitory_cells, expression = Oprd1 > 0)

inhibitory_cells_oprk1_neg <- WhichCells(inhibitory_cells, expression = Oprk1 > 0)

inhibitory_cells_oprm1_neg <- WhichCells(inhibitory_cells, expression = Oprm1 > 0)



inhibitory_cells_oprd1 <- subset(scrna.mast_filtered, subset = (cell.class == "Inhibitory" & Oprd1 >0))
inhibitory_cells_oprk1 <- subset(scrna.mast_filtered, subset = (cell.class == "Inhibitory" & Oprk1 >0))
inhibitory_cells_oprm1 <- subset(scrna.mast_filtered, subset = (cell.class == "Inhibitory" & Oprm1 >0))



VlnPlot(inhibitory_cells_oprk1, features = c('Pdyn'), group.by = 'cell.class', pt.size = 1) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )  +
  scale_y_continuous(limits = c(0, 5))

VlnPlot(inhibitory_cells_oprk1, features = c('Tac1'), group.by = 'cell.class', pt.size = 1) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )  +
  scale_y_continuous(limits = c(0, 5))

VlnPlot(inhibitory_cells_oprk1, features = c('Gal'), group.by = 'cell.class', pt.size = 1) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )  +
  scale_y_continuous(limits = c(0, 5))






VlnPlot(inhibitory_cells_oprd1, features = c('Penk'), group.by = 'cell.class', pt.size = 1)+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 5))




VlnPlot(inhibitory_cells_oprd1, features = c('Gal'), group.by = 'cell.class', pt.size = 1)+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 5))





VlnPlot(inhibitory_cells_oprm1, features = c('Penk'), group.by = 'cell.class', pt.size = 1)+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 5))



VlnPlot(inhibitory_cells_oprm1, features = c('Nts'), group.by = 'cell.class', pt.size = 1)+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 5))


VlnPlot(inhibitory_cells_oprm1, features = c('Gal'), group.by = 'cell.class', pt.size = 1)+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = c(0, 5))




pbmc <- FindVariableFeatures(scrna.mast_filtered, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(scrna.mast_filtered)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1 + plot2