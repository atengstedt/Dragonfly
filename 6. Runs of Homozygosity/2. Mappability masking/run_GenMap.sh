# Make mappability mask of genome

##########################
# Check job info in batch
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
echo $id
done | sort | uniq -c
##########################

#Requires the following packages
#genmap
#bedtools

#index reference genome
folder="/faststorage/project/Coregonus/Aja/Dragonfly/GenMap/"
mkdir ${folder}

sbatch -A Coregonus -t 12:00:00 --mem 50G --job-name genmap_index --wrap\
 "genmap index -F /faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa -I ${folder}/index"
 
#compute mappability
folder="/faststorage/project/Coregonus/Aja/Dragonfly/GenMap/"

sbatch -A Coregonus -t 12:00:00 --mem 100G --job-name genmap_map --wrap\
 "genmap map -K 100 -E 2 -I ${folder}/index -O ${folder} -bg"

#convert bedgraph to bed file (lists all regions with mappability < 1, which we want to remove)
awk -F$'\t' 'BEGIN { OFS="\t" } { if ($4 + 0.0 < 1) print $1, $2, $3}' "/faststorage/project/Coregonus/Aja/Dragonfly/GenMap/dragonfly.contig.genmap.bedgraph" > "/faststorage/project/Coregonus/Aja/Dragonfly/GenMap/dragonfly.contig.genmap.bed"

#merge book-ended features in bed file (increases speed of subsequent filtering process)
bedtools merge -i "/faststorage/project/Coregonus/Aja/Dragonfly/GenMap/dragonfly.contig.genmap.bed" > "/faststorage/project/Coregonus/Aja/Dragonfly/GenMap/dragonfly.contig.genmap.merged.bed"
