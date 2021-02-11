#### Packages ####
library(gplots)
library(RColorBrewer)
library(dplyr)
library(glue)

#### Functions ####
remove_na_dist <- function(mat) {
  dist.mat <- as.matrix(dist(mat))
  
  distNAs <- sum(is.na(dist.mat))
  
  if (distNAs > 0) {
    giveNAs <- which(is.na(dist.mat),
                     arr.ind=TRUE)
    
    tab <- sort(table(c(giveNAs)),decreasing=TRUE)
    checkNA <- sapply(1:length(tab),function(i){
      sum(is.na(as.matrix(dist(mat[-as.numeric(names(tab[1:i])), ]))))})
    rmv <- names(tab)[1:min(which(checkNA == 0))]
    
    return(mat[-as.numeric(rmv), ])
    
  } else {
    return(mat)
  }
}

#### Load data ####
prot.table <- read.csv("~/mrc/project/rna-seq/data/Inflammation_LPSinWTonly_Mus musculus_macrophage_foldChange.txt",
                          sep="\t")

#### Z-scoring ####
# Convert to matrix and z-scoring by gene (rows)
z.matrix <- t(scale(t(data.matrix(prot.table[, 3:ncol(z.table)])),
                    center=TRUE, scale=TRUE))

# Set gene symbol as matrix row names
rownames(z.matrix) <- prot.table[, 2]

#### Filtering ####
# Remove rows and columns with all NA
z.matrix <- z.matrix[rowSums(is.na(z.matrix)) != ncol(z.matrix),
                     colSums(is.na(z.matrix)) < nrow(z.matrix)]

# Remove rows giving NA in distance matrix for clustering
z.filt.matrix <- remove_na_dist(z.matrix)
print(glue("Genes remaining: {nrow(z.filt.matrix)} out of {nrow(prot.table)}"))

filt.gene.ids <- prot.table[match(rownames(z.filt.matrix), prot.table$symbol), c("gene_id")]

#### Heatmap ####
# Set up colour palette
palette <- colorRampPalette(c("blue", "white", "red"))(n=299)

# Plot and save static heatmap
png("~/mrc/project/rna-seq/processed/heatmap.png",
    width=3000, height=9000, res=300)
heatmap.2(z.filt.matrix,
          main="Fold-changes of proteostasis genes in mouse BMDM studies",
          Rowv=TRUE, Colv=TRUE,
          srtCol=45,
          labRow=FALSE,
          notecol="black",
          density.info="none",
          trace="none",
          col=palette,
          na.color="black",
          key=TRUE, keysize=1, lhei=c(1, 15), key.xlab="z-score")

dev.off()

## Interactive heatmap
library(heatmaply)

# Make custom hovertext

# Make interactive heatmap
heatmaply(z.filt.matrix,
          k_row=5,
          dist_method="euclidean", hclust_method="complete",
          seriate="mean",
          colors=palette,
          labRow=rownames(z.filt.matrix),
          row_dend_left=TRUE,
          plot_method="plotly",
          file="~/mrc/project/rna-seq/processed/heatmap.html")

browseURL("~/mrc/project/rna-seq/processed/heatmap.html")

