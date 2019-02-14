# S4 class for storing data from single-cell experiments, including methods to store and retrieve
# spike-in information, dimensionality reduction coordinates and size factors for each cell, 
# along with the usual metadata for genes and libraries

# The SingleCellExperiment class is that it is capable of storing data in normal or sparse
# matrix format, as well as HDF5 format, which allows large non-sparse matrices to be stored 
# & accessed on disk in an efficient manner rather than loading the whole thing into RAM.


library("SingleCellExperiment")

# Generate matrix of sampled random counts (raw count data)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)

# Genes as rows and cells as columns
rownames(counts) <- paste("gene", 1:10, sep = "")
colnames(counts) <- paste("cell", 1:10, sep = "")

# S4 object
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  rowData = data.frame(gene_names = paste("gene_name", 1:10, sep = "")),
  colData = data.frame(cell_names = paste("cell_name", 1:10, sep = ""))
)

sce

# Access sce count matrix
dim(counts(sce))
counts(sce)
