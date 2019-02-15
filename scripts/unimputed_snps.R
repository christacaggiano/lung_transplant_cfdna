rm(list=ls())
setwd("~/Documents/UCSF_year2/research/lung_transplant/")

# install.packages("data.table")
# BiocManager::install("biomaRt", version = "3.8")

require(data.table)
require(biomaRt)

snps = fread("~/Documents/UCSF_year2/research/lung_transplant/data/ddcfDNA.csv")

snps_transposed = t(snps)
colnames(snps_transposed) = snps_transposed[1, ] 
snps_transposed = snps_transposed[-1, ]          
nrow(snps_transposed)
snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
searchDatasets(mart = snp_mart, pattern = "hsapiens")

ids = row.names(snps_transposed)
snp_attributes = c("refsnp_id", "chr_name", "chrom_start", "allele")

snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", values=ids, mart=snp_mart)
listAttributes(snp_mart)
snps = data.frame(snps_transposed)
snps = setDT(snps, keep.rownames = "refsnp_id")[]

snps[snps==0] = "NA" 
snps[snps==1] = "zero" 
snps[snps==2] = "one"
snps[snps==3] = "two" 

snps[snps=="zero"] = 0 
snps[snps=="one"] = 1 
snps[snps=="two"] = 2 

snps = subset(snps, snps$X1402872==0)

ids = c(snps$refsnp_id)



snp_locations$"chrom_end" = snp_locations$chrom_start+1 
snp_locations$"chr_name_test" = paste0("chr", snp_locations$chr_name)
snp_locations$chr_name = snp_locations$chr_name_test
snp_locations$chr_name_test = NULL 
snp_locations$allele_test = gsub("/", ":", snp_locations$allele)
snp_locations$allele = snp_locations$allele_test
snp_locations$allele_test = NULL 
snp_locations = snp_locations[,c(1,2,3,5,4)]

m = merge(snp_locations, snps, by="refsnp_id")

ind_1402 = m[,c("chr_name","chrom_start","chrom_end", "allele", "X1402871", "X1402872" )]
ind_1403 = m[,c("chr_name","chrom_start","chrom_end", "allele", "X1403631", "X1403632" )]

write.table(ind_1402, "1402872_alleles.txt", row.names = F, col.names = F, sep="\t", quote = F)
write.table(ind_1403, "1403632_alleles.txt", row.names = F, col.names = F, sep="\t", quote = F)



