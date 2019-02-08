

rm(list=ls())
setwd("~/Documents/UCSF_year2/research/lung_transplant/data")

df1 = read.table("140287_allele_counts.txt", sep = "\t", header=F) # imputed
names(df1) = c("CHR", "START", "END", "ALLELE", "DONOR", "RECIP", "A", "B", "OTHER")


df2=read.table("../140287_allele_counts_fwd.txt") # non-imputed
names(df2) = c("CHR", "START", "END", "ALLELE", "DONOR", "RECIP", "A", "B", "OTHER")

m = merge(df2, df1, by=c("CHR", "START"))

df3 = read.table("../1402872_alleles.txt")
names(df3) = c("CHR", "START", "END", "ALLELE", "DONOR", "RECIP") # original 

m = merge(df3, df2, by=c("CHR", "START"))

df1 = m[,c(1, 2, 3, 4, 5, 6, 7, 8, 9)]
names(df1) = c("CHR", "START", "END", "ALLELE", "DONOR", "RECIP", "A", "B", "OTHER")

df1_hom = subset(df1, df1$DONOR ==2 & df1$RECIP==0)
sum(df1_hom$B)
sum(df1_hom$A)
sum(df1_hom$B)/(sum(df1_hom$B) + sum(df1_hom$A))
(sum(df1_hom$B) + sum(df1_hom$A) + sum(df1_hom$OTHER))/nrow(df1_hom)

df1_hom = subset(df1, df1$DONOR ==1 & df1$RECIP==0)
sum(df1_hom$B)
sum(df1_hom$A)
sum(df1_hom$B)/(sum(df1_hom$B) + sum(df1_hom$A))
(sum(df1_hom$B) + sum(df1_hom$A) + sum(df1_hom$OTHER))/nrow(df1_hom)

df1_hom = subset(df1, df1$DONOR ==1 & df1$RECIP==2)
sum(df1_hom$B)
sum(df1_hom$A)
sum(df1_hom$B)/(sum(df1_hom$B) + sum(df1_hom$A))
(sum(df1_hom$B) + sum(df1_hom$A) + sum(df1_hom$OTHER))/nrow(df1_hom)

df1_hom = subset(df1, df1$DONOR==0 & df1$RECIP==2)
sum(df1_hom$B)
sum(df1_hom$A)
sum(df1_hom$A)/(sum(df1_hom$B) + sum(df1_hom$A))
(sum(df1_hom$B) + sum(df1_hom$A) + sum(df1_hom$OTHER))/nrow(df1_hom)

df1_hom = subset(df1, df1$DONOR==0 & df1$RECIP==0)
sum(df1_hom$B)
sum(df1_hom$A)
sum(df1_hom$A)/(sum(df1_hom$B) + sum(df1_hom$A))
(sum(df1_hom$B) + sum(df1_hom$A))/nrow(df1_hom)

df1_hom = subset(df1, df1$DONOR==2 & df1$RECIP==2)
sum(df1_hom$B)
sum(df1_hom$A)
nrow(df1_hom)
sum(df1_hom$B)/(sum(df1_hom$B) + sum(df1_hom$A))
(sum(df1_hom$B) + sum(df1_hom$A))/nrow(df1_hom)

df1_hom = subset(df1, df1$DONOR==1 & df1$RECIP==1)
sum(df1_hom$B)
sum(df1_hom$A)

sum(df1_hom$B)/(sum(df1_hom$B) + sum(df1_hom$A))
(sum(df1_hom$B) + sum(df1_hom$A) + sum(df1_hom$OTHER))/nrow(df1_hom)

m_pos = paste(m$CHR, m$START)
df3_pos = paste(df3$CHR, df3$START)

length(intersect(m_pos, df3_pos))

m_pos[duplicated(m_pos)==TRUE]
df3_pos[duplicated(df3_pos)==TRUE]

vars = read.table("../140287_vars2.txt", header=F)

df4 = read.table("../140287_hg38.txt", header = F)
names(df4) = c("CHR", "START", "END", "ALLELE", "DONOR", "RECIP")

df3_pos =  paste(df3$CHR, df3$START)
df4_pos =  paste(df4$CHR, df4$START)


l = df3_pos%in%df4_pos
length(l[l==TRUE])

nrow(df4)
nrow(df3)

df1_pos =  paste(df1$CHR, df1$START)
df2_pos =  paste(df2$CHR, df2$START)

a = setdiff(df1_pos, df2_pos)
a
length(df2_pos[a==FALSE])+length(df2_pos[a==TRUE])

subset


hom4 = subset(df4, df4$RECIP==2)
hom3 = subset(df3, df3$RECIP==2)

df1_pos =  paste(hom4$CHR, hom4$START)
df2_pos =  paste(hom3$CHR, hom3$START)

l = df1_pos%in%df2_pos
length(l[l==TRUE])

rm(list=ls())
require(data.table)

snps = fread("~/Documents/UCSF_year2/research/lung_transplant/data/ddcfDNA.csv")
unimputed_snps = t(snps)
nrow(unimputed_snps)
colnames(unimputed_snps) = unimputed_snps[1, ] 
unimputed_snps = unimputed_snps[-1, ]

unimputed = data.frame(unimputed_snps)
unimputed = setDT(unimputed, keep.rownames = "refsnp_id")[]

unimputed[unimputed==0] = "NA" 
unimputed[unimputed==1] = "zero" 
unimputed[unimputed==2] = "one"
unimputed[unimputed==3] = "two" 

unimputed[unimputed=="zero"] = 0 
unimputed[unimputed=="one"] = 1 
unimputed[unimputed=="two"] = 2 

unimputed = subset(unimputed, unimputed$X1402872==1)

unimputed_ids = c(unimputed$refsnp_id)



imputed = fread("~/Desktop/zaitlen_lab_desktop/140287_vars.csv", sep = ",")
nrow(imputed)
imputed_ids = strsplit(imputed$SNP, ":")
imputed = sapply(strsplit(imputed$SNP, ":"), `[`, 1)

length(imputed_ids)
length(unimputed_ids)

length(union(imputed_ids, unimputed_ids))
length(intersect(imputed_ids, unimputed_ids))
length(setdiff(imputed_ids, unimputed_ids))
length(setequal(imputed_ids, unimputed_ids))

i = unimputed_ids%in%imputed_ids
length(i[i==TRUE])

m = merge(imputed, )




