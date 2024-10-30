##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#bwa
#samtools
#r-base

#===========Index reference genome
ref="/faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"

sbatch -A Coregonus -t 12:00:00 --mem 8G --wrap\
 "bwa index $ref"

#===========Map reads for each individual
ref="/faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"
from="/faststorage/project/Coregonus/Backup/Dragonfly/"
to="/faststorage/project/Coregonus/Aja/Dragonfly/mem"

mkdir $to

cat individuals.txt | while read individual
do
sbatch -A Coregonus -t 12:00:00 --mem 16G -c 16 --wrap\
 "bwa mem -M -t 16 $ref ${from}/${individual}x1.fq.gz ${from}/${individual}x2.fq.gz > ${to}/${individual}.sam"
done

#===========Convert SAM to BAM, clean up and sort
folder="/faststorage/project/Coregonus/Aja/Dragonfly/mem"

mkdir ${folder}/tmp1
mkdir ${folder}/tmp2

cat individuals.txt | while read individual
do
sbatch -A Coregonus -t 24:00:00 --job-name samtools --mem 16G --wrap\
 "samtools sort -n -T ${folder}/tmp1/${individual} ${folder}/${individual}.sam | \
 samtools fixmate - - -m | \
 samtools sort - -T ${folder}/tmp2/${individual} \
 -O bam -o ${folder}/${individual}.bam"
done

#===========Remove PCR duplicates
folder="/faststorage/project/Coregonus/Aja/Dragonfly/mem"

cat individuals.txt | while read individual
do
sbatch -A Coregonus -t 12:00:00 --job-name markdup --mem 8G -c 4 --wrap\
 "samtools markdup -@ 4 -r ${folder}/${individual}.bam ${folder}/${individual}_purged.bam"
done

#===========Index BAM files
folder="/faststorage/project/Coregonus/Aja/Dragonfly/mem"

cat individuals.txt | while read individual
do
sbatch -A Coregonus -t 12:00:00 --job-name index --mem 8G --wrap\
 "samtools index ${folder}/${individual}_purged.bam"
done

#===========Stats
from="/faststorage/project/Coregonus/Aja/Dragonfly/mem"
to="/faststorage/project/Coregonus/Aja/Dragonfly/STATS_mem"

mkdir $to

for bam in `ls $from | grep "_purged.bam$"`
do
sbatch -A Coregonus -t 12:00:00 --wrap "samtools flagstat ${from}/${bam} > ${to}/${bam/.bam/.flagstat}"
sbatch -A Coregonus -t 12:00:00 --wrap "samtools stats -c 1,1000,1 ${from}/${bam} > ${to}/${bam/.bam/.stats}"
done

#-----------Summary of flagstat
echo -e "Indi\tTotal\tPrimary\tSec\tSup\tDup\tPDup\tMapped\tPMapped\tPaired\tR1\tR2\tPMapped\tIMMapped\tSingle\tMdiff\tMapQ>=5" > STATS_mem/zzflagstat
for file in `ls STATS_mem | grep "\\.flagstat$"`
do
echo -e -n ${file/.flagstat/}"\t" >> STATS_mem/zzflagstat
awk 'NR==1{printf($1);next}{printf("\t%d",$1)}END{printf("\n")}' STATS_mem/${file}\
 >> STATS_mem/zzflagstat
done

#-----------Summary (R)
# Run "STATS.r"


