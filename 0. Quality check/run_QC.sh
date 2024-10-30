##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#samtools
#fastqc
#multiqc
#r-base

#===========Create a list of individual names
from="/faststorage/project/Coregonus/Dragonfly_genomes/"

ls -d ${from}/*/ | sed 's#.*/\([^/]*\)/$#\1#' > "/faststorage/project/Coregonus/Aja/Dragonfly/individuals.txt"

#OBS: Merged and compressed reads (produced as described below) have been moved to Backup folder, where they can be accessed directly for mapping

#===========Merge lane files for each individual read (forward/reverse)
#from="/faststorage/project/Coregonus/Dragonfly_genomes/"
#to="/faststorage/project/Coregonus/Aja/Dragonfly/reads"

#mkdir $to

#cat individuals.txt | while read line
#do
#sbatch -A Coregonus -t 12:00:00 --wrap "cat ${from}/${line}/*_1.fq > ${to}/${line}x1.fq"
#sbatch -A Coregonus -t 12:00:00 --wrap "cat ${from}/${line}/*_2.fq > ${to}/${line}x2.fq"
#done

#===========Compress read files
#folder="/faststorage/project/Coregonus/Aja/Dragonfly/reads/"

#cat individuals.txt | while read line
#do
#sbatch -A Coregonus -t 12:00:00 --wrap "bgzip ${folder}/${line}x1.fq"
#sbatch -A Coregonus -t 12:00:00 --wrap "bgzip ${folder}/${line}x2.fq"
#done

#OBS: Merged and compressed reads have been moved to Backup folder, where they can be accessed directly for mapping

#===========FastQC
from="/faststorage/project/Coregonus/Backup/Dragonfly/" # Backup folder
to="/faststorage/project/Coregonus/Aja/Dragonfly/QC"

mkdir $to

for file in ${from}/*.fq
do
sbatch -A Coregonus -t 12:00:00 --wrap\
 "fastqc -o $to $file"
done

#==========MultiQC (aggrevate output of fastQC into one report)
cd "/faststorage/project/Coregonus/Aja/Dragonfly/"

multiqc .

#-----------Info from FastQC result (run in R)
setwd("/faststorage/project/Coregonus/Aja/Dragonfly/QC/")

files<-grep("\\_fastqc.html$",dir(),value=T)

n<-length(files)

l.read<-integer(n)
n.read<-l.read
encoding<-character(n)

for(i in 1:n)
{
	html<-scan(files[i],what="",quiet=T,sep="\n")
	html<-grep("<td>Sequence length</td>",html,value=T)
	a<-unlist(strsplit(html,"<[[:alpha:]/]{1,}>"))
	l.read[i]<-a[which(a=="Sequence length")+2]
	n.read[i]<-a[which(a=="Total Sequences")+2]
	encoding[i]<-a[which(a=="Encoding")+2]
}
info<-cbind(n.read,l.read,encoding)
rownames(info)<-files
write.table(info,"info",sep="\t",quote=F)


#===========Merge lane files for individual used for reference genome
sbatch -A Coregonus -t 12:00:00 --wrap "cat Coregonus/Aja/Dragonfly_harddisk/DRAfijwR/CleanData/*_1.fq.gz > Coregonus/Backup/Dragonfly/L210930x1-FULL.fq.gz"
sbatch -A Coregonus -t 12:00:00 --wrap "cat Coregonus/Aja/Dragonfly_harddisk/DRAfijwR/CleanData/*_2.fq.gz > Coregonus/Backup/Dragonfly/L210930-FULLx2.fq.gz"

#===========Subsample reads to get same coverage as other individuals (average number of reads across the other individuals)
sbatch -A Coregonus -t 12:00:00 --mem 200G --wrap "seqtk sample -s100 /faststorage/project/Coregonus/Backup/Dragonfly/L210930-FULLx1.fq.gz 136645215 > /faststorage/project/Coregonus/Backup/Dragonfly/NHMA-ENT-221278x1.fq"
sbatch -A Coregonus -t 12:00:00 --mem 200G --wrap "seqtk sample -s100 /faststorage/project/Coregonus/Backup/Dragonfly/L210930-FULLx2.fq.gz 136645215 > /faststorage/project/Coregonus/Backup/Dragonfly/NHMA-ENT-221278x2.fq"

#-----------Zip subsets
sbatch -A Coregonus -t 12:00:00 --wrap "bgzip /faststorage/project/Coregonus/Backup/Dragonfly/NHMA-ENT-221278x1.fq"
sbatch -A Coregonus -t 12:00:00 --wrap "bgzip /faststorage/project/Coregonus/Backup/Dragonfly/NHMA-ENT-221278x2.fq"

#===========FastQC
from="/faststorage/project/Coregonus/Backup/Dragonfly/"
to="/faststorage/project/Coregonus/Aja/Dragonfly/QC"

sbatch -A Coregonus -t 12:00:00 --wrap\
 "fastqc -o $to /faststorage/project/Coregonus/Backup/Dragonfly/NHMA-ENT-221278x1.fq.gz"
sbatch -A Coregonus -t 12:00:00 --wrap\
 "fastqc -o $to /faststorage/project/Coregonus/Backup/Dragonfly/NHMA-ENT-221278x2.fq.gz"
