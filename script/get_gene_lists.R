#### Packages ####
library(glue)
library(dplyr)

#### Config ####
setwd("~/mrc/project/rna-seq")

#### Load data ####
gene.set <- read.table(glue("processed/GSE154573/GSE154573_DESeq_KO_vs_WT.txt"), header=TRUE, sep='\t',
                       row.names=1, check.names=FALSE)
background.set <- read.table(glue("processed/GSE154573/GSE154573_fpkm_filtered.txt"), header=TRUE, sep='\t',
                             row.names=1, check.names=FALSE)

#### Filtering ####
alpha <- 0.05
upreg.set <- gene.set %>% filter(padj < alpha & log2FoldChange > 0)
downreg.set <- gene.set %>% filter(padj < alpha & log2FoldChange < 0)

#### Save to list ####
fileConn <- file("processed/GSE154573/GSE154573_upregset.txt")
writeLines(rownames(upreg.set), fileConn)
close(fileConn)

fileConn <- file("processed/GSE154573/GSE154573_downregset.txt")
writeLines(rownames(downreg.set), fileConn)
close(fileConn)

fileConn <- file("processed/GSE154573/GSE154573_backgroundset.txt")
writeLines(rownames(background.set), fileConn)
close(fileConn)


