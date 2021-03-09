#### Packages ####
library(dplyr)
library(ggplot2)
library(ggrepel)
library(glue)

#### Parameters ####
alpha <- 0.05

gse <- "GSE143241"
ko <- "TTP"

select.columns <- c(glue('{gse}_WT'), glue('{gse}_{ko}'))
# GSE97538 - c('GSE97538_WT', 'GSE97538_A20')
# GSE143241 - c('GSE143241_WT', 'GSE143241_TTP')
# GSE134443 - c('GSE134443_WT', 'GSE134443_TRPM5')

#### Load data ####
fc.table <- read.csv("~/mrc/project/rna-seq/data/Inflammation_LPS_EISA_Mus musculus_foldChange.txt",
                     sep="\t")
padj.table <- read.csv("~/mrc/project/rna-seq/data/Inflammation_LPS_EISA_Mus musculus_padj.txt",
                       sep="\t")

select.columns <- c(c('gene_id', 'symbol'), select.columns)
fc.table <- fc.table[, which(names(fc.table) %in% select.columns)]
padj.table <- padj.table[, which(names(padj.table) %in% select.columns)]

# Convert to matrix
fc.matrix <- data.matrix(fc.table[, 3:ncol(fc.table)])

# Set gene symbol as matrix row names
rownames(fc.matrix) <- fc.table[, 2]

#### Filtering ####
# Convert non-significant gene fold-changes with p > alpha to 0
fc.matrix[which(padj.table[, -c(1:2)] > alpha)] <- 0
fc.matrix[is.na(padj.table[, -c(1:2)])] <- 0

# Remove rows with all 0
fc.matrix <- fc.matrix[rowSums(fc.matrix) != 0, ]

fc.filt.table <- data.frame(fc.matrix)
colnames(fc.filt.table) <- c("WT_log2FC", "KO_log2FC")

fc.filt.table$difference <- fc.filt.table$KO_log2FC - fc.filt.table$WT_log2FC
fc.filt.table$direction <- ifelse(fc.filt.table$difference > 0, "UP", "DOWN")
fc.filt.table$WT_padj <- subset(padj.table, padj.table$symbol %in% rownames(fc.filt.table))[, 3]
fc.filt.table$KO_padj <- subset(padj.table, padj.table$symbol %in% rownames(fc.filt.table))[, 4]
fc.filt.table$symbol <- rownames(fc.filt.table)

rownames(fc.filt.table) <- NULL

up.genes <- fc.filt.table[which(fc.filt.table$direction == 'UP'), ]
down.genes <- fc.filt.table[which(fc.filt.table$direction == 'DOWN'), ]

#### Visualisation ####
ggplot(fc.filt.table, aes(x=WT_log2FC, y=KO_log2FC)) + 
  geom_abline(intercept=0, 
              slope=1,
              linetype="dashed", size=0.5, color="gray") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(aes(fill=direction),
             colour="black", pch=21) +
  geom_text_repel(data=. %>% mutate(label=ifelse(abs(difference) > 2, symbol, "")),
                  aes(label=label),
                  size=3, segment.alpha=0.5,
                  max.overlaps=Inf,
                  box.padding=0.5) +
  labs(title=glue("{gse}: Comparing fold-changes in LPS treatment in {ko} vs. WT"),
       subtitle=glue("UP={nrow(up.genes)}, DOWN={nrow(down.genes)}"),
       x="log2FoldChange (WT)", 
       y=glue("log2FoldChange ({ko} KO)"),
       fill="Direction (KO vs. WT)") +
  theme(legend.title=element_text(size=8))

