# Read-based phasing using Whatshap combined with population-level phasing with shapeit
# Preparation for analysis of LD decay

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#vcftools
#samtools
#bcftools
#whatshap
#shapeit4

#--------------Extract each individual to its own vcf file
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/individuals"

mkdir ${output}

cat individuals2.txt | while read line
do
sbatch -A Coregonus -t 12:00:00 --job-name indv --wrap\
 "vcftools --vcf ${input} --indv ${line} --recode --out ${output}/${line}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr"
done

#===============Phasing with Whatshap
ref="/faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"
vcf="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/individuals/"
bam="/faststorage/project/Coregonus/Aja/Dragonfly/mem/"
folder="/faststorage/project/Coregonus/Aja/Dragonfly/whatshap"

mkdir ${folder}

cat individuals2.txt | while read line
do
sbatch -A Coregonus -t 12:00:00 --mem 12G --job-name whatshap --wrap\
 "whatshap phase -o ${folder}/${line}.whatshap.vcf --ignore-read-groups --sample=${line} --reference=${ref} ${vcf}/${line}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.recode.vcf ${bam}/${line}_purged.bam"
done

#--------------Zip individual vcf files
folder="/faststorage/project/Coregonus/Aja/Dragonfly/whatshap"

cat individuals2.txt | while read line
do
sbatch -A Coregonus -t 12:00:00 --job-name bgzip --wrap\
 "bgzip ${folder}/${line}.whatshap.vcf"
done

#--------------Index individual vcf files
folder="/faststorage/project/Coregonus/Aja/Dragonfly/whatshap"

cat individuals2.txt | while read line
do
sbatch -A Coregonus -t 12:00:00 --job-name tabix --wrap\
 "tabix ${folder}/${line}.whatshap.vcf.gz"
done

#--------------Merge individual VCF files into one
folder="/faststorage/project/Coregonus/Aja/Dragonfly/whatshap"

sbatch -A Coregonus -t 12:00:00 --mem 8G -c 8 --job-name merge --wrap "bcftools merge ${folder}/*.whatshap.vcf.gz --threads 8 -O z -o ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.whatshap.vcf.gz"

#--------------Index merged vcf file
file="/faststorage/project/Coregonus/Aja/Dragonfly/whatshap/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.whatshap.vcf.gz"

sbatch -A Coregonus -t 12:00:00 --job-name tabix --wrap "tabix ${file}"

#==============Phasing with Shapeit4 (separate scaffolds/chromosomes) 
input="/faststorage/project/Coregonus/Aja/Dragonfly/whatshap/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.whatshap.vcf.gz"
output="/faststorage/project/Coregonus/Aja/Dragonfly/shapeit"

mkdir ${output}
mkdir ${output}/log

cat scaffold_names.100kb_noSexChr.txt | while read scaffold
sbatch -A Coregonus -t 12:00:00 --job-name shapeit --mem 4G -c 4 --wrap\
 "shapeit4 --input ${input} --region $scaffold --use-PS 0.0001 --sequencing --thread 4 --log ${output}/log/${line}.shapeit.log --output ${output}/${scaffold}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.shapeit.vcf.gz"
done

# flag --use-PS 0.0001 informs phasing using PS field from read based phasing, with an expected error rate of 0.0001 (default/recommended in guidelines for SHAPEIT4)
