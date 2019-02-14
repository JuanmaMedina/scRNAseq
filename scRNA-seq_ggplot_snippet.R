library(ggplot2)
library(tidyverse)
library(pheatmap)

set.seed(1)

# Generate matrix of sampled random counts
counts <- as.data.frame(matrix(rpois(100, lambda = 10), ncol=10, nrow=10))

Gene_ids <- paste("gene", 1:10, sep = "")
colnames(counts) <- paste("cell", 1:10, sep = "")

counts <- data.frame(Gene_ids, counts)
counts

# Scatter plot: no correlation between gene expression in cells 1 and 2
ggplot(data = counts, mapping = aes(x = cell1, y = cell2)) + geom_point()

# Boxplot of every cell and their average expression across all genes
counts <- gather(counts, colnames(counts)[2:11], key = 'Cell_ID', value='Counts')
counts

ggplot(counts, aes(x=Cell_ID, y=Counts)) + geom_boxplot()

set.seed(2)

# Counts data
test <- matrix(rnorm(200), 20, 10)

# Generate variation
test[1:10, seq(1, 10, 2)] <- test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] <- test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] <- test[15:20, seq(2, 10, 2)] + 4

colnames(test) <- paste("Cell", 1:10, sep = "")
rownames(test) <- paste("Gene", 1:20, sep = "")
test

pheatmap(test)
