#### Packages ####
library(DESeq2)
library(glue)
library(biomaRt)

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
  
  row.names(geneList) <- geneList[, 1]
  geneList[, 1] <- NULL
  
  table$geneSymbol <- geneList[, 1][match(rownames(table), rownames(geneList))]
  newTable <- table
  
  return(newTable)
}

#### Load data ####
proj <- "GSE147893"
countsTable <- read.table(glue("processed/{proj}_rawCounts.txt"), header=TRUE, sep='\t',
                          row.names=1, check.names=FALSE)
counts <- countsTable
counts$Length <- NULL

#### DESeq2 ####
# Create DDS object
colData <- data.frame(row.names=colnames(counts),
                      condition=colnames(counts))
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=colData,
                              design= ~ 1)

# Add gene length info
mcols(dds)$basepairs <- countsTable[, 'Length']

# Calculate FPKM
dds.fpkm <- fpkm(dds)

# Add gene symbol info
# dds.fpkm <- add_ensembl_symbol(dds.fpkm)

# Save to output
write.table(dds.fpkm, file=glue("processed/{proj}_fpkm.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
