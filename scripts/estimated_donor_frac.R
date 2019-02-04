rm(list=ls())
setwd("~/Desktop/zaitlen_lab_desktop/")

df1 = read.table("140363_allele_counts.txt", sep = "\t", header=F)
names(df1) = c("CHR", "START", "END", "ALLELE", "DONOR", "RECIP", "A", "B", "OTHER")

sum(df1$B)

df1$REF_FRAC = df1$A/(df1$A + df1$B + df1$OTHER)
df1$total_count = df1$A + df1$B + df1$OTHER

sum(df1$total_count)/nrow(df1)

df1 = subset(df1, df1$REF_FRAC!="NaN")
df1 = subset(df1, df1$OTHER!=1)
nrow(df1)

df1$DONOR = round(df1$DONOR)
df1$RECIP = round(df1$RECIP)


df1_donor_homB = subset(df1, df1$DONOR==2)
df1_donor_homA = subset(df1, df1$DONOR==0)
df1_donor_het = subset(df1, df1$DONOR==1)
nrow(df1_donor_het)

hetA = sum(df1_donor_het$A)
hetB = sum(df1_donor_het$B)
homA_correct = sum(df1_donor_homA$A)
homB_correct = sum(df1_donor_homB$B)
homA_incorrect = sum(df1_donor_homA$B)
homB_incorrect = sum(df1_donor_homB$A)

hom_correct = homA_correct + homB_correct
hom_incorrect = homA_incorrect + homB_incorrect

(2*(hetB) + hom_correct)/(hetA+hetB+hom_correct+hom_incorrect)

(sum(df1$OTHER))/(hetA+hetB+hom_correct+hom_incorrect)


