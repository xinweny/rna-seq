#### Packages ####
library(gplots)
library(RColorBrewer)
library(dplyr)
library(glue)
library(dendextend)

#### Functions ####
remove_na_dist <- function(mat) {
  dist.mat <- as.matrix(dist(mat))
  
  distNAs <- sum(is.na(dist.mat))
  
  if (distNAs > 0) {
    giveNAs <- which(is.na(dist.mat),
                     arr.ind=TRUE)
    
    tab <- sort(table(c(giveNAs)), 
                decreasing=TRUE)
    
    checkNA <- sapply(1:length(tab), function(i) {
      sum(is.na(as.matrix(dist(mat[-as.numeric(names(tab[1:i])), ]))))})
    
    rmv <- names(tab)[1:min(which(checkNA == 0))]
    
    return(mat[-as.numeric(rmv), ])
    
  } else {
    return(mat)
  }
}

#### Parameters ####
alpha <- 0.05

delete.columns <- c('GSE122292', 'GSE110243', 'GSE123596', 'GSE98563', 'GSE82043')

#### Load data ####
fc.table <- read.csv("~/mrc/project/rna-seq/data/Inflammation_LPSinWTonly_Mus musculus_macrophage_foldChange.txt",
                          sep="\t")
padj.table <- read.csv("~/mrc/project/rna-seq/data/Inflammation_LPSinWTonly_Mus musculus_macrophage_padj.txt",
                       sep="\t")

# Remove custom columns
if (length(delete.columns) > 0) {
  fc.table <- fc.table[, -which(names(fc.table) %in% delete.columns)]
  padj.table <- padj.table[, -which(names(padj.table) %in% delete.columns)]
}

# Convert to matrix
fc.matrix <- data.matrix(fc.table[, 3:ncol(fc.table)])

# Set gene symbol as matrix row names
rownames(fc.matrix) <- fc.table[, 2]

#### Filtering ####
# Convert non-significant gene fold-changes with p < alpha to NA
fc.matrix[which(padj.table[, -c(1:2)] > alpha)] <- NA

# Remove rows and columns with all NA
fc.matrix <- fc.matrix[rowSums(is.na(fc.matrix)) != ncol(fc.matrix),
                       colSums(is.na(fc.matrix)) != nrow(fc.matrix)]

# Rank genes by mean fold-change across GSEs
# fc.matrix[is.na(fc.matrix)] <- 0
# fc.rank <- sort(rowMeans(fc.matrix), decreasing=TRUE)

# Extract genes based on ranking
# fc.matrix <- fc.matrix[which(rownames(fc.matrix) %in% names(head(fc.rank, 5))), ]

# Remove rows giving NA in distance matrix for clustering
fc.filt.matrix <- remove_na_dist(fc.matrix)
print(glue("Genes remaining: {nrow(fc.filt.matrix)} out of {nrow(fc.table)}"))

filt.gene.ids <- fc.table[match(rownames(fc.filt.matrix), fc.table$symbol), c("gene_id")]

#### Heatmap ####
# Set up colour palette
palette <- colorRampPalette(c("blue", "white", "red"))(n=299)

# Clustering on rows and columns
Rowv <- fc.filt.matrix %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=5) %>% set("branches_lwd", 4) %>%
  rotate_DendSer(ser_weight=dist(fc.filt.matrix))
Colv <- fc.filt.matrix %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=2) %>% set("branches_lwd", 4) %>%
  rotate_DendSer(ser_weight=dist(t(fc.filt.matrix)))

# Plot and save static heatmap
png(file="~/mrc/project/rna-seq/processed/heatmap.png", 
    width=3000, height=9000, res=300)
hm <- heatmap.2(fc.filt.matrix,
                scale="none",
                main=glue("Fold-changes of proteostasis genes (n={nrow(fc.filt.matrix)})
                in mouse BMDM studies (N={ncol(fc.filt.matrix)}),
                          p < {alpha}"),
                Rowv=Rowv, Colv=TRUE,
                srtCol=45,
                labRow=FALSE, # row.names(fc.filt.matrix)
                notecol="black",
                density.info="none",
                trace="none",
                col=palette,
                symbreaks=TRUE,
                na.color="black",
                key=TRUE, keysize=1, lhei=c(1, 15), key.xlab="log2(fold-change)")
hm
dev.off()

#### Interactive heatmap ####
library(heatmaply)

# Make interactive heatmap
heatmaply(fc.filt.matrix,
          dist_method="euclidean", hclust_method="complete",
          scale="none",
          Rowv=Rowv, Colv=TRUE,
          seriate="mean",
          colors=palette,
          limits=c(-15, 15),
          labRow=rownames(fc.filt.matrix),
          row_dend_left=TRUE,
          na.rm=FALSE, na.value="black",
          key.title="log2(foldChange)",
          plot_method="plotly",
          branches_lwd=0.1,
          main=glue("Fold-changes of proteostasis genes (n={nrow(fc.filt.matrix)}, p < {alpha}) in mouse BMDM studies (N={ncol(fc.filt.matrix)})"),
          width=2000, height=4000,
          file="~/mrc/project/rna-seq/processed/heatmap.html")

browseURL("~/mrc/project/rna-seq/processed/heatmap.html")

#### Extract clusters ####
row.clusts <- sort(cutree(fc.filt.matrix %>% dist %>% hclust, k=5))

table(row.clusts)

select.clusts <- c(4, 5)
select.genes <- names(row.clusts[which(row.clusts %in% select.clusts)])

write.table(fc.table[which(fc.table$symbol %in% select.genes), ],
            file=glue("~/mrc/project/rna-seq/processed/heatmap_upreg_genes.txt"),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
