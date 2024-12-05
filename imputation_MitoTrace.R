library(MitoTrace)

listlcl <- list.files("GTExbam/LCL", full.names = T, pattern = ".bam$")
bams <- listlcl
fasta_loc <- "GTExbam/LCL/Homo_sapiens.GRCh38.dna.chromosome.MT.fa" # read reference mitochondrial sequence

mae_res <- MitoTrace(bam_list = listlcl, fasta = fasta_loc, chr_name = "MT")

af <- calc_allele_frequency(mae_res)
colnames(af) <- unlist(lapply(colnames(af), function(x) { a <- strsplit(x, "\\.") [[1]][1]}))

af[af<0.9] <- 0

af[af>=0.9] <- 1

afmat <- rowSums(af)/ncol(af)

af.sel <- af[which(afmat>0.05),]
saveRDS(af, "afall.RDS")
saveRDS(af.sel, "afsel.RDS")
write.csv(af, "afall.csv")
write.table(af, "afall.txt", quote=F, sep = "\t")
