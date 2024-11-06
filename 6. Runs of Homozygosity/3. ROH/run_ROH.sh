##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#bcftools
#vcftools
#plink
#r-base
#r-detectruns

#===========Filter to retain only the 100 longest scaffolds
list="/faststorage/project/Coregonus/Aja/Dragonfly/scaffold_names.100longest_noSexChr.bed"
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf.gz"
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.vcf"

sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name bcftools --wrap\
 "bcftools view -R $list $input -O v -o $output"
 
#2,750,262 SNPs out of 7,724,101
#2,768,393 SNPs out of 7,790,422 (Dragonfly2 with reference individual)

#===========Discard low mappability regions
sbatch -A Coregonus -t 12:00:00 --wrap\
 "vcftools --vcf /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.vcf --exclude-bed /faststorage/project/Coregonus/Aja/Dragonfly/GenMap/dragonfly.contig.genmap.merged.bed --recode --recode-INFO-all --out /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask"

#2,424,500 SNPs
#2,437,876 SNPs (Dragonfly2 with reference individual)

#===========Discard repetitive regions
sbatch -A Coregonus -t 12:00:00 --wrap\
 "vcftools --vcf /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.recode.vcf --exclude-bed /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/AesVir_genomic.repeats.merged.bed --recode --recode-INFO-all --out /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask"
 
#1,628,513 SNPs
#1,636,441 SNPs (Dragonfly2 with reference individual)

#==========Remove sites with missing data
sbatch -A Coregonus -t 12:00:00 --wrap\
 "vcftools --vcf /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.recode.vcf --max-missing 1 --recode --recode-INFO-all --out /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1"
 
#1,121,350 SNPs (length of 100 scaffolds = 566 Mb, so approx 1 SNP per 500 bp)
#1,117,911 SNPs (Dragonfly2 with reference individual)

#--------------convert VCF to ped/map
file="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1.recode.vcf"
ID_file="/faststorage/project/Coregonus/Aja/Dragonfly/plink/FID2.txt"
out="/faststorage/project/Coregonus/Aja/Dragonfly/ROH/input_files/"

mkdir $out

sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name conv_pedmap --wrap\
 "plink --vcf ${file} --const-fid --allow-extra-chr --update-ids $ID_file --het --recode --out $out/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1"

#remove prefix from scaffold names
sed -i 's/Contig//g' "/faststorage/project/Coregonus/Aja/Dragonfly/ROH/input_files/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1.map"

#Move into ROH folder

#submit job
folder="/faststorage/project/Coregonus/Aja/Dragonfly/ROH/"
cd ${folder}

for roh in run1D run1E run2D run2E #run run0 run1A run1B run1C
do
mkdir $roh
sbatch -A Coregonus -t 12:00:00 --mem 16G --job-name ${roh} --wrap\
 "Rscript ${roh}.R"
done
