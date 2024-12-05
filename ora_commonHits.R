library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)

# Script for pathways over-representation analyses
geneNames <- readRDS("common_to_all_hits.RDS")

gm <- readRDS("allGenes.RDS") # intersection of all genes measured in the three datasets
b <- bitr(geneNames$ensembl_gene_id, fromType="ENSEMBL", toType=("ENTREZID"), OrgDb=org.Hs.eg.db) # map ENSEMBL IDs to ENTREZ IDs
bu <- bitr(gm$ensembl_gene_id, fromType="ENSEMBL", toType=("ENTREZID"), OrgDb=org.Hs.eg.db) # map ENSEMBL IDs to ENTREZ IDs
k <- enrichKEGG(as.character(b$ENTREZID),  universe =as.character(bu$ENTREZID), organism = "hsa")
kdf <- k@result # KEGG pathways enrichment results

go <- enrichGO(as.character(geneNames$ensembl_gene_id),  universe =as.character(gm$ensembl_gene_id), ont = "BP", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
godf <- go@result # GO "BP" pathways enrichment results
saveRDS(go, "gobp.RDS")

go <- enrichGO(as.character(geneNames$ensembl_gene_id),  universe =as.character(gm$ensembl_gene_id), ont = "MF", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
godfmf <- go@result # GO "MF" pathways enrichment results
go <- enrichGO(as.character(geneNames$ensembl_gene_id),  universe =as.character(gm$ensembl_gene_id), ont = "CC", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
godfcc <- go@result # GO "CC" pathways enrichment results


saveRDS(godf, "godf.RDS")
saveRDS(godfmf, "godfmf.RDS")
saveRDS(godfcc, "godfcc.RDS")
saveRDS(kdf, "keggdf.RDS")


# create dotplots
library(clusterProfiler)

godf <- readRDS("gobp.RDS")

dotplot(godf)
cnetplot(godf)
barplot(godf)
treeplot(godf)
emapplot(godf)

# create an emap plot
library(enrichplot)
edo <- pairwise_termsim(godf)
e <- emapplot(edo, showCategory=5)
pdf("gobp_emapplot.pdf")
print(e)
dev.off()
