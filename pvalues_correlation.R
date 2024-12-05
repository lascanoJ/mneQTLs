library(ggpmisc)
library(ggpubr)
library(ggplot2)

# Only considering mtSNP-gene pairs shared by all three datasets
dataset1 <- geuvadis_eur[commonset_all, ]
dataset2 <- gtex[commonset_all, ]
dataset3 <- geuvadis_yru[commonset_all, ]

# bind p-values in one dataframe
data <- rbind(dataset1$pvalue, dataset2$pvalue)
data2 <- rbind(dataset1$pvalue, dataset3$pvalue)
data <- cbind(data, data2)
data <- as.data.frame(t(data))
data <- -log10(data)
data$group <- c(rep("GTEX", 9), rep("Geuvadis Yru", 9))

g1 <- ggplot(data, aes(x=V1, y=V2, group = group)) +
  geom_smooth(method="lm", aes(color = group), alpha=0.3) + 
  geom_point( aes(x=V1, y=V2, color=group)) +
  stat_cor() + xlab("Geuvadis Eur dataset") + ylab("Other dataset")



# Datasets pairwise: considering mtSNP-gene pairs shared by two datasets
dataset1 <- geuvadis_eur[commonset1, ]
dataset2 <- gtex[commonset1, ]

# bind p-values in one dataframe

data <- cbind(-log10(dataset1$pvalue), -log10(dataset2$pvalue))
data <- as.data.frame(data)
g2 <- ggplot(data, aes(x=V1, y=V2)) +
  geom_smooth(method="lm", alpha=0.3) + 
  geom_point( aes(x=V1, y=V2)) +
  stat_cor() + xlab("Geuvadis Eur dataset") + ylab("GTEx dataset")

dataset1 <- geuvadis_eur[commonset2, ]
dataset3 <- geuvadis_yru[commonset2, ]

# bind p-values in one dataframe

data <- cbind(-log10(dataset1$pvalue), -log10(dataset3$pvalue))
data <- as.data.frame(data)
g3 <- ggplot(data, aes(x=V1, y=V2)) +
  geom_smooth(method="lm", alpha=0.3) + 
  geom_point( aes(x=V1, y=V2)) +
  stat_cor() + xlab("Geuvadis Yru dataset") + ylab("GTEx dataset")



dataset2 <- gtex[commonset3, ]
dataset3 <- geuvadis_yru[commonset3, ]

# bind p-values in one dataframe

data <- cbind(-log10(dataset3$pvalue), -log10(dataset2$pvalue))
data <- as.data.frame(data)
g4 <- ggplot(data, aes(x=V1, y=V2)) +
  geom_smooth(method="lm", alpha=0.3) + 
  geom_point( aes(x=V1, y=V2)) +
  stat_cor() + xlab("Geuvadis Yru dataset") + ylab("Geuvadis Eur dataset")




pdf("pvaluesShared_all.pdf")
print(g1)
dev.off()

pdf("pvaluesShared_geuveurgtex.pdf")
print(g2)
dev.off()

pdf("pvaluesShared_geuvyrugtex.pdf")
print(g3)
dev.off()

pdf("pvaluesShared_geuvyrueur.pdf")
print(g4)
dev.off()


