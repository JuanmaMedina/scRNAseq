# Count matrix from kidney cells sequenced with FACS + Smartseq2 
# (low hightroughput and high coverage)
dat <- read.delim("Jupyter/scRNAseq_course/tabula_muris/facs/FACS/Kidney-counts.csv", 
                 sep=",", header=TRUE)
dat[1:5,1:5]

# Genes as row names
rownames(dat) <- dat[,1]
dat <- dat[,-1]

# Check for spike-ins that Smartseq2 can contain
rownames(dat)[grep("^ERCC-", rownames(dat))]

# Metadata: well, plate and mouse IDs
cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\.")
Well <- unlist(lapply(cell_info, function(x){x[1]}))
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))

# Check distribution of mice and plates
summary(factor(Mouse))
summary(factor(Plate))

# Check for confounding factors: all mice with the same ID were analyzed in the same plate
# The experimental batch (PCR plate) is completely confounded with donor mice
table(Mouse, Plate)

# Cell-type annotations
ann <- read.table("Jupyter/scRNAseq_course/tabula_muris/facs/FACS_annotations.csv", 
                  sep=",", header=TRUE)

# Kidney cells annotations
ann <- ann[match(cellIDs, ann[,1]),]

celltype <- ann[,3]

# QC and visualization analyses of scRNA-seq gene expression data
library("scater")
library("SingleCellExperiment")

# All cell annotations in the same DF, with cells as col names
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
rownames(cell_anns) <- colnames(dat)

# SingleCell object
sceset <- SingleCellExperiment(
  assays = list(counts = as.matrix(dat)), 
  colData = cell_anns)

# Variable to track spike-ins
isSpike(sceset, "ERCC") <- grepl("ERCC-", rownames(sceset))

