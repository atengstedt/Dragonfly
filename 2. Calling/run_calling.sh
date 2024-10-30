##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#bcftools
#r-base
#vcftools
#samtools

#===========Call variants for each scaffold
mkdir "/faststorage/project/Coregonus/Aja/Dragonfly/VCF"

ref="/faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"
from="/faststorage/project/Coregonus/Aja/Dragonfly/mem"
to="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/raw_scaffolds2"

mkdir $to

cat scaffold_names_ref.txt | while read line 
do
sbatch -A Coregonus -t 4:00:00 --job-name mpileup_call --wrap\
 "bcftools mpileup -Ou -r $line -q 30 -Q 25 --annotate FORMAT/DP,FORMAT/AD -f $ref ${from}/*_purged.bam | bcftools call -f GQ -v -m -O v -o ${to}/${line}.vcf"
done

#===========Merge the 1499 scaffolds >100 kb in length
in="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/raw_scaffolds2"

cd ${in} 

out="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_raw.vcf"
name_file="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/raw_scaffolds/vcflist.txt"

sbatch -A Coregonus -t 12:00:00 --mem 8G -c 8 --begin now+1hour --wrap "bcftools concat --threads 8 -f ${name_file} -O v -o ${out}"

#-----------Line count (how many variants remain?)
file="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_raw.vcf"

sbatch -A Coregonus -t 12:00:00 --wrap "bcftools view -H ${file} | wc -l"

#14,227,015 SNPs 
#14,241,806 SNPs (Dragonfly2 with reference individual)

#===========Inspect SNP depth distribution
awk '($0!~/^#/)&&($8~/^DP/){split($8,a,"[=;]");print a[2]}'\
 Dragonfly2_raw.vcf > Dragonfly2_raw.snp.depth

#-----------Plot (R)
a<-scan("/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_raw.snp.depth")
hist(a[a<4000],200,main=NULL,xlab="SNP depth")
abline(v=600,col="red")  #choose lower threshold
abline(v=1900,col="red")  #choose upper threshold
savePlot("/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_SNP_depth.tiff","png")

#===========Filter (indels, monomorphic+multi-allelic sites, depth, MapQ)
#use thresholds determined from plot
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_raw.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered"

sbatch -A Coregonus -t 12:00:00 --job-name vcf.sh --wrap\
 "vcfutils.pl varFilter -Q 30 -d 600 -D 1900 ${input} | vcftools --vcf - --remove-indels --max-alleles 2 --maf 0.000001 --recode --recode-INFO-all --out ${output}"
 
#8,943,371 SNPs
#8,995,078 SNPs (Dragonfly 2 with reference individual)

#===========Reheader (names are originally shown as e.g. mem/P192139_purged.bam and I want it to just be P192139)
file="/faststorage/project/Coregonus/Aja/Dragonfly/individuals2.txt"
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.recode.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.reheader.vcf"

sbatch -A Coregonus -t 12:00:00 --wrap\
 "bcftools reheader -s $file $input > $output"
 
#mv reheadered file to original path
mv ${output} ${input}

#--------Filter for minGQ and minDP to exclude low-quality SNPs and SNPs with very high or very low read depth
#--------Removes no sites, just changes offending genotypes to missing (./.)
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF"

sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name vcftools --wrap\
 "vcftools --vcf ${folder}/Dragonfly2_filtered.recode.vcf --minGQ 20 --minDP 6 --maxDP 45 --recode --recode-INFO-all --out ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45"
 
#===========Filter for HWE using Lins script
#https://github.com/shenglin-liu/VCF_HWF
sbatch -A Coregonus -t 12:00:00 --mem 100G --wrap\
 "Rscript /faststorage/home/anbt/software/VCF_HWE/VCF_HWF.r /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.recode.vcf /faststorage/project/Coregonus/Aja/Dragonfly/HWE/popmap2_region.txt 6 /faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf"

#------------Extract header from original VCF file
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.recode.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.header.vcf"

sbatch -A Coregonus -t 12:00:00 --wrap\
 "bcftools view --header-only $input > $output"

#------------Append new header to HWE filtered VCF
cat "/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.header.vcf" "/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf" > "/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.new-header.vcf"

#mv reheadered file to original path
mv "/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.new-header.vcf" "/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf"

#8,688,180 SNPs
#8,769,291 SNPs (Dragonfly 2 with reference individual)

###Run preliminary PCA (run_PCA_prelim.sh) to check population structure (presence of sex chromosomes)

#zip and index VCF
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/"
bgzip ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf
tabix ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf.gz

#--------Exclude scaffolds located on sex chromosome
list="/faststorage/project/Coregonus/Aja/Dragonfly/scaffold_names.100kb_noSexChr.bed"
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf.gz"
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.vcf"

sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name bcftools --wrap\
 "bcftools view -R $list $input -O v -o $output"
 
#7,724,101 SNPs
#7,790,422 SNPs (Dragonfly 2 with reference individual)

#=========Annotate SNP IDs in VCF file
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/"

sbatch -A Coregonus -t 12:00:00 --mem 8G --wrap\
 "bcftools annotate ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.vcf --set-id '%CHROM\_%POS' -O v -o ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.ann.vcf"

#move annotated file to original path
#mv {folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.ann.vcf ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.vcf
