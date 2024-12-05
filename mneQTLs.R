library(MatrixEQTL)
library(matrixStats)
library(dplyr)

###### Mitonuclear eQTL calculations

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use
useModel = modelLINEAR; 

# Genotype file name
SNP_file_name = "geno.txt"
snps_location_file_name = "snploc.txt"

# Gene expression file name
expression_file_name = "expr.txt";
gene_location_file_name = "geneloc.txt";

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# We calculate all associations

pvOutputThreshold_tra = 1;

# Error covariance matrix
errorCovariance = numeric();

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      
snps$fileOmitCharacters = "NA"; 
snps$fileSkipRows = 1;          
snps$fileSkipColumns = 1;       
snps$fileSliceSize = 2000;      
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1;         
gene$fileSkipColumns = 1;       
gene$fileSliceSize = 2000;      
gene$LoadFile("expr.txt");

## Load covariates
covariates_file_name="covmat.txt"
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      
cvrt$fileOmitCharacters = "NA";
cvrt$fileSkipRows = 1;          
cvrt$fileSkipColumns = 1;    
if(length(covariates_file_name)>0) {
  cvrt$LoadFile("covmat.txt");
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
set.seed(123)

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)


eqtl <- me$all$eqtls
eqtl$pos <- as.numeric(sapply(strsplit(as.character(eqtl$snps), "A|T|C|G"), `[`, 1))

eqtl$hit <- paste(eqtl$pos, eqtl$gene)

p <- plot(me)
pdf("qqplotGtex.pdf")
plot(me)
dev.off()
saveRDS(eqtl, "allNominalAssociations.RDS")


saveRDS(eqtl, "allNominal_geuv_y.RDS")
