#### Packages ####
library(DESeq2)
library(glue)
library(biomaRt)
library(dplyr)

setwd("~/mrc/project/rna-seq")

#### Functions ####
add_ensembl_symbol <- function (table) {
  genes <- row.names(table)
  
  if (grepl("ENSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "hsapiens_gene_ensembl"
    symbol <- "hgnc_symbol"
  } else if (grepl("ENSMUSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "mmusculus_gene_ensembl"
    symbol <- "mgi_symbol"
  }
  
  mart <- useDataset(ensemblDataset, useMart("ENSEMBL_MART_ENSEMBL", host="http://www.ensembl.org"))
  geneList <- getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id", symbol),
                    values=genes,
                    mart=mart)
  
  geneList <- distinct(geneList, ensembl_gene_id, .keep_all=TRUE)
  
  row.names(geneList) <- geneList[, 1]
  geneList[, 1] <- NULL
  
  table$geneSymbol <- geneList[, 1][match(rownames(table), rownames(geneList))]
  newTable <- table
  
  return(newTable)
}

format_condition <- function (colnames) {
  replace <- c("_[0-9]*$", "_rep[0-9]*$", "^[A-Z]{3}[0-9]+_", "^[0-9]+_")
  
  for (r in replace) {
    colnames <- gsub(r, "", colnames)
  }
  
  return(colnames)
}

#### Load data ####
proj <- "PROJ1672"
counts <- read.table(glue("processed/{proj}_rawCounts.txt"), header=TRUE, sep='\t',
                          row.names=1, check.names=FALSE)

# Design
group <- "B3KO"
control <- "B3KO_T0"
treatment <- "B3KO_T3"

alpha <- 1e-10
lfcThresh <- 0

# counts$Length <- NULL
counts <- cbind(gene=rownames(counts), counts)
counts$gene <- gsub("\\.[0-9_A-Z]+$", "", counts$gene)

# Aggregate genes with same name
counts <- aggregate(counts[, -1], by=list(gene=counts$gene), FUN=sum)
rownames(counts) <- counts$gene
counts$gene <- NULL

# Group counts
counts.grp <- dplyr::select(counts, matches(group))

#### DESeq2 ####
# Create DDS object
colData <- data.frame(row.names=colnames(counts.grp),
                      condition=format_condition(colnames(counts.grp)))
dds <- DESeqDataSetFromMatrix(countData=counts.grp,
                              colData=colData,
                              design=~condition)

dds$condition <- relevel(dds$condition, ref=glue(control))

# PCA plot
rld <- vst(dds, blind=TRUE)
plotPCA(rld)

# Remove low count data
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", treatment, control),
                alpha=alpha, lfcThreshold=lfcThresh)

res <- add_ensembl_symbol(res)

# Save DESeq results table to output
deGenes <- as.data.frame(res) %>% arrange(padj, desc(log2FoldChange)) # order by adjusted p-value and FC
write.table(deGenes,
            file=glue("processed/{proj}_DESeq_{treatment}_vs_{control}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#### Visualisation ####
nUp <- nrow(filter(as.data.frame(res), padj < alpha & log2FoldChange > 0))
nDown <- nrow(filter(as.data.frame(res), padj < alpha & log2FoldChange < 0))

png(glue("processed/{proj}_MAplot_{treatment}_vs_{control}.png"))
DESeq2::plotMA(res, 
               main=glue("{proj}: {treatment} vs. {control}
          n={nUp + nDown}, UP={nUp}, DOWN={nDown}
          alpha={alpha}"))
dev.off()

