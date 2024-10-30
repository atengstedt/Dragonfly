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
#r-base
#samtools
#pixy

#============VCF with monomorphic sites
#script uses the same pooled depth thresholds determined in run_calling.sh
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2"
script="/faststorage/project/Coregonus/Aja/Dragonfly/VCF_filter_2.sh"

mkdir $folder

while read line
do
sbatch -A Coregonus -t 100:00:00 $script $line
done < scaffold_names.100kb_noSexChr.txt

#---------Filter for minGQ and minDP (doesn't remove sites, just changes offending genotypes to ./.)
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/"

cat scaffold_names.100kb_noSexChr.txt | while read line
do
sbatch -A Coregonus -t 24:00:00 --wrap\
 "vcftools --vcf ${folder}/${line}.vcf --minGQ 20 --minDP 6 --maxDP 45 --recode --recode-INFO-all --out ${folder}/${line}.minDP6.maxDP45.minGQ20"
done

#===========Filter for HWE
#https://github.com/shenglin-liu/VCF_HWF
cat scaffold_names.100kb_noSexChr.txt | while read line
do
sbatch -A Coregonus -t 12:00:00 --mem 10G --wrap\
 "Rscript /faststorage/home/anbt/software/VCF_HWE/VCF_HWF.r /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/${line}.minDP6.maxDP45.minGQ20.recode.vcf /faststorage/project/Coregonus/Aja/Dragonfly/HWE/popmap2_region.txt 6 /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/${line}.minDP6.maxDP45.minGQ20.HWE6.vcf"
done

#------------Extract header from one of the original VCF files
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/Contig10003.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/Contig10003.header.vcf"

sbatch -A Coregonus -t 12:00:00 --wrap\
 "bcftools view --header-only $input > $output"

#------------Append new header to HWE filtered VCFs
cat scaffold_names.100kb_noSexChr.txt | while read line
do
sbatch -A Coregonus -t 12:00:00 --mem 8G --wrap\
 "cat /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/Contig10003.header.vcf /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/${line}.minDP6.maxDP45.minGQ20.HWE6.vcf > /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/${line}.minDP6.maxDP45.minGQ20.HWE6.new-header.vcf"
done

#mv reheader files to original paths
cat scaffold_names.100kb_noSexChr.txt | while read line
do
mv /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/${line}.minDP6.maxDP45.minGQ20.HWE6.new-header.vcf /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/${line}.minDP6.maxDP45.minGQ20.HWE6.vcf
done

#===========Merge scaffold VCF files
in="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2"

cd ${in} 

out="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/AllSites2_filtered.minDP6.maxDP45.minGQ20.HWE6.vcf.gz"
name_file="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/vcflist.txt"

sbatch -A Coregonus -t 24:00:00 --mem 8G -c 32 --wrap "bcftools concat --threads 32 -f ${name_file} -O z -o ${out}"

#---------Get overview of missing data per individual
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/"

sbatch -A Coregonus -t 12:00:00 --wrap\
 "vcftools --vcf ${folder}/AllSites2_filtered.minDP6.maxDP45.minGQ20.HWE6.vcf --missing-indv --out ${folder}/AllSites2_filtered.minDP6.maxDP45.minGQ20.HWE6"
 
#1,239,636,584 sites in total

#-----------Calculate HO per scaffold.
script="/faststorage/project/Coregonus/Aja/Dragonfly/VCF_HO_2.r"
folder="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2"

cat scaffold_names.100kb_noSexChr.txt | while read line 
do
sbatch -A Coregonus -t 12:00:00 --job-name calc_HO --mem 64G --wrap "Rscript ${script} ${folder}/${line}.minDP6.maxDP45.minGQ20.HWE6.vcf"
done

cat calc_HO*.out | grep "^/faststorage" > VCF2_HO_2.txt

#=============Plot HO from VCF (R)
par(mar=c(10,4,2,2))

a<-read.table("VCF2_HO_2.txt",sep="\t",stringsAsFactors=F)

popmap<-read.table("individuals2.fixed.txt",stringsAsFactors=F)

#popmap[,1]<-sub("^mem/(.*)_sorted.bam$","\\1",popmap[,1])

HO<-as.matrix(a[,-c(1:2)])

colnames(HO)<-popmap[,1]

rownames(HO)<-a[,1]

HO<-colSums(HO*a[,2])/sum(a[,2])

write.table(HO, "VCF2_HO_2-values.txt")

barplot(HO,las=3)

savePlot("VCF2_HO_2.tiff","png") 

#--------------index merged file
sbatch -A Coregonus -t 12:00:00 --wrap\
 "tabix /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/AllSites_filtered.minDP6.maxDP45.minGQ20.HWE6.vcf.gz"

#==============Nucleotide diversity calculated with pixy
output="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/pixy10kb"

#mkdir $output

sbatch -A Coregonus -t 96:00:00 --mem 24G --job-name pixy --wrap\
 "pixy --stats pi --vcf /faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/AllSites2_filtered.minDP6.maxDP45.minGQ20.HWE6.vcf.gz \
 --populations ${output}/populations2.txt \
 --window_size 10000 --output_folder ${output} --output_prefix AllSites2_filtered.minDP6.maxDP45.minGQ20.HWE6.pixy"

#calculate pi across the entire genome by summing column "count_diffs" and dividing by the sum of "count_comparisons"
#as described on the pixy webpage https://pixy.readthedocs.io/en/latest/output.html

awk 'NR>1 { sum7[$1]+=$7; sum8[$1]+=$8 } END { for (pop in sum8) { if (sum8[pop] != 0 && sum7[pop] != 0) print pop, sum7[pop]/sum8[pop] } }' "/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2/pixy10kb/AllSites2_filtered.minDP6.maxDP45.minGQ20.HWE6.pixy_pi.txt"


