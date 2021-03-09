#### Packages ####
import pandas as pd
import os
import numpy as np

#### Config ####
os.chdir("/Users/pomato/mrc/project/rna-seq/")
mode = 'geneSymbol' # entrezID

#### Load data ####
go_cc = pd.read_csv(f"./data/GSE154573_KO_vs_WT_{mode}_GO_CC_mgi_functional_classification.tsv", header=0, sep='\t')[['Genes', 'Process~name']]

deseq = pd.read_csv("./processed/GSE154573/GSE154573_DESeq_KO_vs_WT.txt", header=0, sep='\t')

deseq['GO_CC'] = 'NA'
go_cc['Process~name'] = go_cc['Process~name'].apply(lambda x: x.split('~')[1])

deseq['entrezID'] = deseq['entrezID'].astype('str')
deseq['entrezID'].replace('\.0', '', regex=True, inplace=True)
deseq['entrezID'].replace('nan', 'NA', inplace=True)

for i in deseq[mode]:
    if i != 'NA':
        go_terms = list(go_cc[go_cc['Genes'].str.contains(f"\\b{i};", regex=True)]['Process~name'])
        deseq.loc[deseq[mode] == i, 'GO_CC'] = ','.join(go_terms)

if mode == 'geneSymbol':
    deseq.drop(['entrezID'], axis=1, inplace=True)

deseq.to_csv("./processed/GSE154573/GSE154573_DESeq_KO_vs_WT_GOCC.txt", sep='\t', header=True, index=True)