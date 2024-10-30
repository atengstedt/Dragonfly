# Calculate LD decay for the NEZ individuals (largest group)

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State6
done | sort | uniq -c
##########################

#Requires the following packages
#vcftools
#r-base

#-----------Subset longest scaffold VCF file to include only NEZ individuals and apply MAF filter
#requires phased data
input="/faststorage/project/Coregonus/Aja/Dragonfly/shapeit"
output="/faststorage/project/Coregonus/Aja/Dragonfly/LD-decay"

sbatch -A Coregonus -t 12:00:00 --job-name maf --wrap\
 "vcftools --gzvcf ${input}/Contig9961_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.shapeit.vcf.gz --keep /faststorage/project/Coregonus/Aja/Dragonfly/LD-decay/NEZ2023.txt --max-missing 1 --maf 0.15 --recode --out ${output}/Contig9961_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.shapeit.NEZ2023.maf0.15"

# Calculate R-squared
folder="/faststorage/project/Coregonus/Aja/Dragonfly/LD-decay/"

sbatch -A Coregonus -t 12:00:00 --job-name calc_r2 --wrap\
 "vcftools --vcf ${folder}/Contig9961_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.shapeit.NEZ2023.maf0.15.recode.vcf --hap-r2 --out ${folder}/Contig9961.NEZ2023.maf0.15"
 
#================Plot decay of LD using R script from https://jujumaan.com/2017/07/15/linkage-disequilibrium-decay-plot/
folder="/faststorage/project/Coregonus/Aja/Dragonfly/LD-decay/"

sbatch -A Coregonus -t 12:00:00 --mem 164G --job-name rscript --wrap\
 "Rscript ${folder}/plot.R ${folder}/Contig9961.NEZ2023.maf0.15.hap.ld ${folder}/Contig9961.NEZ2023.maf0.15.tiff"

#============Remove sites with missing data and thin VCF so no sites are within <0.28 Mb from each other
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/"

sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name vcftools --wrap\
 "vcftools --vcf ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.vcf --max-missing 1 --thin 280000 --recode --recode-INFO-all --out ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned"
 
#5,835 SNPs (1 SNP per 280.000 bp)
