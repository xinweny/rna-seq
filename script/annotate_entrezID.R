#### Packages ####
library(biomaRt)
library(dplyr)

#### Config ####
setwd("/Users/Pomato/mrc/project/rna-seq")

#### Load data ####
deseq.res <- read.table("./processed/GSE154573/GSE154573_DESeq_KO_vs_WT.txt",
                   sep="\t", header=T)
gene.ids <- rownames(deseq.res)

#### Query ####
mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensembl.to.entrez <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                           filters="ensembl_gene_id",
                           values=gene.ids, mart=mart) %>% 
  distinct(ensembl_gene_id, .keep_all=TRUE)

row.names(ensembl.to.entrez) <- ensembl.to.entrez[, 1]
ensembl.to.entrez[, 1] <- NULL
deseq.res$entrezID <- ensembl.to.entrez[, 1][match(rownames(deseq.res), rownames(ensembl.to.entrez))]

write.table(deseq.res,
            file=glue("./processed/GSE154573/GSE154573_DESeq_KO_vs_WT.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

