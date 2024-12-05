# Bayesian Network of 16319G>A with associated genes
library(bnlearn)
library(dplyr)

intGenes <- c("MEOX2", "NPY2R", "GLP1R", "PDX1")
qn_sel <- qn[which(qn$id %in% intGenes), ] # select gene of interest expression values
rownames(qn_sel) <- qn_sel$id
df <- as.data.frame(t(rbind(qn_sel[, -1], mt_genotype))) # bind gene of interest expression values to the mitochondrial genotype of interest
colnames(df)[5] <- "16319G>A"
# set a black list: no gene to snp relationship is possible
bl = data.frame(from = c("NPY2R", "MEOX2", "PDX1", "GLP1R"), to = c("16319G>A"))
# set a white list: relationship supported by the litterature
wl = data.frame(from = c("GLP1R"), to = c("PDX1"))
set.seed(345)
bn1 <- boot.strength(df, algorithm = "hc", R=1000, algorithm.args = list(blacklist=bl, whitelist=wl)) # we appply 1000 bootstraps for evaluating the strength of conditional dependences
avg <- averaged.network(bn1,threshold = 0.6) # we average across the results from the 1000 bootstraps and keep only those relationships present in at least 60% of bootstraps

pdf("BN.pdf", width=11, height = 9)
strength.plot(avg, bn1) # we visualize the Bayesian Network
dev.off()
write.csv(bn1, "BN.csv")
saveRDS(bn1, "BN.RDS")

