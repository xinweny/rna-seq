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
cluster.cols <- FALSE
k.cluster.rows <- FALSE

delete.columns <- c()
# BMDM - c('GSE122292', 'GSE110243', 'GSE123596', 'GSE98563', 'GSE82043')
# microglia - c('GSE117646', 'GSE79898', 'GSE98142', 'GSE105155', 'GSE153419', 'GSE75246')
select.columns <- c()
gse <- 'PROJ1742'
# GSE97538 - c('GSE97538_WT', 'GSE97538_A20')
# GSE143241 - c('GSE143241_WT', 'GSE143241_TTP')
# GSE134443 - c('GSE134443_WT', 'GSE134443_TRPM5')

heatmap.filename <- glue("~/mrc/project/rna-seq/processed/heatmap_{gse}.png")

#### Load data ####
fc.table <- read.csv(glue("~/mrc/project/rna-seq/data/{gse}_foldChange.txt"),
                          sep="\t")
padj.table <- read.csv(glue("~/mrc/project/rna-seq/data/{gse}_padj.txt"),
                       sep="\t")

# Remove/Select custom columns
if (length(delete.columns) > 0) {
  fc.table <- fc.table[, -which(names(fc.table) %in% delete.columns)]
  padj.table <- padj.table[, -which(names(padj.table) %in% delete.columns)]
}

if (length(select.columns) > 0) {
  select.columns <- c(c('gene_id', 'symbol'), select.columns)
  fc.table <- fc.table[, which(names(fc.table) %in% select.columns)]
  padj.table <- padj.table[, which(names(padj.table) %in% select.columns)]
}

# Convert to matrix
fc.matrix <- data.matrix(fc.table[, 3:ncol(fc.table)])

# Set gene symbol as matrix row names
rownames(fc.matrix) <- fc.table[, 2]

#### Filtering ####
# Convert non-significant gene fold-changes with p < alpha to NA
fc.matrix[which(padj.table[, -c(1:2)] > alpha)] <- 0
fc.matrix[is.na(padj.table[, -c(1:2)])] <- 0

# Remove rows and columns with all NA or 0
fc.matrix <- fc.matrix[rowSums(is.na(fc.matrix)) != ncol(fc.matrix),
                       colSums(is.na(fc.matrix)) != nrow(fc.matrix)]

fc.matrix <- fc.matrix[rowSums(fc.matrix) != 0, ]

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
k.rows <- 6
k.cols <- 2

Rowv <- fc.filt.matrix %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=k.rows) %>% 
  set("labels_col", k=k.rows) %>%
  set("branches_lwd", 4) %>%
  rotate_DendSer(ser_weight=dist(fc.filt.matrix))

Colv <- fc.filt.matrix %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=k.cols) %>% 
  set("branches_lwd", 4) %>%
  rotate_DendSer(ser_weight=dist(t(fc.filt.matrix)))

# Plot and save static heatmap
if (length(select.columns) == 0) {
  png(file=heatmap.filename, 
      width=3000, height=9000, res=300)
  hm <- heatmap.2(fc.filt.matrix,
                  scale="none",
                  main=glue("Fold-changes of proteostasis genes (n={nrow(fc.filt.matrix)})
                  over PROJ1742 timepoints (N={ncol(fc.filt.matrix)}),
                  p < {alpha}"),
                  Rowv=if (k.cluster.rows) Rowv else TRUE,
                  Colv=cluster.cols,
                  dendrogram=if (cluster.cols) "both" else "row",
                  srtCol=45, offsetCol=0.5,
                  margins=c(10, 5),
                  labRow=FALSE, # row.names(fc.filt.matrix)
                  notecol="black",
                  density.info="none",
                  trace="none",
                  col=palette,
                  symbreaks=TRUE,
                  na.color="black",
                  key=TRUE, keysize=1, lhei=c(1, 15), key.xlab="log2(foldChange)")
  hm
  dev.off()
} else {
  png(file=heatmap.filename, 
      width=3000, height=9000, res=300)
  hm <- heatmap.2(fc.filt.matrix,
                  scale="none",
                  main=glue("Fold-changes of proteostasis genes (n={nrow(fc.filt.matrix)})
                  in {gse} LPS treatment in WT and KO,
                  p < {alpha}"),
                  Rowv=TRUE, Colv=FALSE,
                  dendrogram='row',
                  srtCol=45, offsetCol=0.5, cexCol=2,
                  margins=c(12, 5),
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
}

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
row.clusts <- cutree(fc.filt.matrix %>% dist %>% hclust, k=k.rows)

row.clusts.cols <- labels_colors(Rowv)[row.clusts %>% sort %>% names]
table(row.clusts.cols)

select.clusts <- c("#A352D1", "#008FB7", "#009232")
select.genes <- names(row.clusts.cols[which(row.clusts.cols %in% select.clusts)])

write.table(fc.table[which(fc.table$symbol %in% select.genes), ],
            file=glue("~/mrc/project/rna-seq/processed/heatmap_BMDM_upreg_genes.txt"),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
