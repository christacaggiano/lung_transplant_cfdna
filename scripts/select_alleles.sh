#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r y                            
#$ -j y
#$ -l mem_free=10G                 
#$ -l arch=linux-x64               
#$ -l netapp=5G,scratch=5G         
#$ -l h_rt=336:00:00

conda activate py36
python select_alleles.py 140287_allele_counts.txt processed/1/unsortedButMerged_ForBismark_file/1.sorted.bam 140287_hg38.txt
python select_alleles.py 140363_allele_counts.txt processed/2/unsortedButMerged_ForBismark_file/2.sorted.bam 140363_hg38.txt