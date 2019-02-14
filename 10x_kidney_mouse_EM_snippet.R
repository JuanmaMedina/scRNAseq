# Count matrix from kidney cells sequenced with 10x data 
# (high hightroughput and low coverage)

library("Matrix")

# Barcode sequences associated with each cell (col names)
cellbarcodes5 <- read.table("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P4_5/barcodes.tsv")
cellbarcodes6 <- read.table("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P4_6/barcodes.tsv")
cellbarcodes7 <- read.table("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P7_5/barcodes.tsv")

# Genes (row names)
genenames5 <- read.table("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P4_5/genes.tsv")
genenames6 <- read.table("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P4_6/genes.tsv")
genenames7 <- read.table("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P7_5/genes.tsv")

# Due to the large size and sparsity of 10X data (up to 90% of the expression matrix may be 0s), 
# it is typically stored as a sparse matrix, stored in a .mtx file by Cell Ranger as a column of
# row coordinates, a column of column corodinates and a column of expression values > 0. 
molecules5 <- Matrix::readMM("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P4_5/matrix.mtx")
molecules6 <- Matrix::readMM("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P4_6/matrix.mtx")
molecules7 <- Matrix::readMM("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet/Kidney-10X_P7_5/matrix.mtx")

# Having the barcode sequences associated with each cell is a problem since each batch of 10x
# data uses the same pool of barcodes, so if we need to combine data from multiple 10x batches, 
# the cellbarcodes will not be unique. Hence, attach the batch ID to each cell barcode:
rownames(molecules5) <- genenames5[,1]
rownames(molecules6) <- genenames6[,1]
rownames(molecules7) <- genenames7[,1]

colnames(molecules5) <- paste("10X_P4_5", cellbarcodes5[,1], sep="_")
colnames(molecules6) <- paste("10X_P4_6", cellbarcodes6[,1], sep="_")
colnames(molecules7) <- paste("10X_P7_5", cellbarcodes7[,1], sep="_")

# Metadata
meta <- read.delim("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet_metadata.csv",
                   sep=",", header=TRUE)

# Fix mice IDs (from hyphens to underscores)
meta[meta$channel == "10X_P4_5" | meta$channel == "10X_P4_6" | meta$channel == "10X_P7_5",]
mouseID_38M <- "3_8_M"
mouseID_39M <- "3_9_M"
mouseID_357F <- "3_57_F"

# Annotations
ann <- read.delim("Jupyter/scRNAseq_course/tabula_muris/droplet/droplet_annotations.csv",
                  sep=",", header=TRUE)

# Correct formatting difference between the cellID in the annotation and the cellbarcodes
ann[,1] <- paste(ann[,1], "-1", sep="")

ann_subset5 <- ann[match(colnames(molecules5), ann[,1]),]
ann_subset6 <- ann[match(colnames(molecules6), ann[,1]),]
ann_subset7 <- ann[match(colnames(molecules7), ann[,1]),]

celltype5 <- ann_subset5[,3]
celltype6 <- ann_subset6[,3]
celltype7 <- ann_subset7[,3]

# Cell-metadata DF
cell_anns5 <- data.frame(mouse = rep(mouseID_38M, times=ncol(molecules5)), type=celltype5)
cell_anns6 <- data.frame(mouse = rep(mouseID_39M, times=ncol(molecules6)), type=celltype6)
cell_anns7 <- data.frame(mouse = rep(mouseID_357F, times=ncol(molecules7)), type=celltype7)

rownames(cell_anns5) <- colnames(molecules5)
rownames(cell_anns6) <- colnames(molecules6)
rownames(cell_anns7) <- colnames(molecules7)

# Sanity checks
# Genes are the same and in the same order across all batches
identical(rownames(molecules5), rownames(molecules6))
identical(rownames(molecules5), rownames(molecules7))

# No repeated cell IDs
sum(colnames(molecules5) %in% colnames(molecules6))
sum(colnames(molecules5) %in% colnames(molecules7))

# Combine them
all_molecules <- cbind(molecules5, molecules6, molecules7)
all_cell_anns <- as.data.frame(rbind(cell_anns5, cell_anns6, cell_anns7))
all_cell_anns$batch <- rep(c("10X_P4_5", "10X_P4_6","10X_P7_5"), times = c(nrow(cell_anns5), nrow(cell_anns6), nrow(cell_anns7)))


# SingleCell object
library("SingleCellExperiment")
library("scater")

all_molecules <- as.matrix(all_molecules)
sceset <- SingleCellExperiment(
  assays = list(counts = as.matrix(all_molecules)), 
  colData=all_cell_anns)

# Save the data
# saveRDS(sceset, "kidney_droplet.rds")

