##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#bcftools
#seqtk
#psmc
#r-base

#===========Sequences
ref="/faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"
from="/faststorage/project/Coregonus/Aja/Dragonfly/mem"
to="/faststorage/project/Coregonus/Aja/Dragonfly/SEQ"
fdep="/faststorage/project/Coregonus/Aja/Dragonfly/STATS_mem/zstats"

mkdir ${to}

samtools faidx $ref

cat individuals2.txt | while read individual
do
dep=`grep -P "^${individual}_purged\t" $fdep | awk '{printf("-d %.0f -D %.0f\n",$2/3,$2*2)}'`
sbatch -A Coregonus -t 48:00:00 --job-name sequences --wrap\
 "bcftools mpileup -Ou -q 30 -Q 25 -f $ref ${from}/${individual}_purged.bam\
 | bcftools call -c - | vcfutils.pl vcf2fq $dep > ${to}/${individual}.fq"
done

#----------Include only long scaffolds (>1 Mb)
folder="/faststorage/project/Coregonus/Aja/Dragonfly/SEQ"

cd $folder

for NameRoot in `ls | grep "fq$" | cut -d . -f 1 | sort | uniq`
do
sbatch -A Coregonus -t 12:00:00 --job-name seqtk --wrap\
 "seqtk subseq ${NameRoot}.fq ../scaffold_names.1Mb_noSexChr.txt > ${NameRoot}.1Mb_noSexChr.fq"
done

#===========Make PSMC input files
input="/faststorage/project/Coregonus/Aja/Dragonfly/SEQ"
folder="/faststorage/project/Coregonus/Aja/Dragonfly/PSMC"

mkdir ${folder}

for file in `ls ${input} | grep "1Mb_noSexChr.fq$"`
do
sbatch -A Coregonus --job-name fq2psmcfa --wrap "fq2psmcfa -q20 ${input}/${file} > ${folder}/${file/.fq/.psmcfa}"
done

#-----------Run PSMC (bin size 100)
folder="/faststorage/project/Coregonus/Aja/Dragonfly/PSMC"

for file in ${folder}/*.psmcfa
do
sbatch -A Coregonus -t 12:00:00 --mem 16G --job-name psmc --wrap\
 "psmc -N25 -t15 -r5 -p 1+1+1+1+25*2+4+6 -o ${file/psmcfa/1+1+1+1.psmc} $file"
done

#----------Plot PSMC with outgroup, bin size 100 (default) and point estimate mutation rate (R)
setwd("PSMC")
source("plot_PSMC.r")

#read popmap
popmap<-read.table("popmap2_region.txt",sep="\t",stringsAsFactors=F)

# Create a list of vectors - one for each pop
vector_list <- lapply(split(popmap$V1, popmap$V2), function(x) as.character(x))

#make plot colorful
pops<-c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")
x<-unique(pops)
n<-length(x)
col.palette<-c("#9b5fe0","#16a4d8","#60dbe8","#8bd346","#efdf48","#f9a52c","#d64e12")
col.palette<-t(as.data.frame(col.palette))
colnames(col.palette)<-x

#create file
pdf(file="plots/PSMC2.2.9e-9.2+2.pdf", width = 169/25.4, height = 150/25.4)

#create layout
nf<-layout(matrix(c(1,2,3,4,5,6,7,0,0), nrow=3, ncol=3, byrow=T))

poplist<-c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")

for (pop in poplist){
files <- vector_list[[pop]]  #assign vector from list to variable "files"
col<-col.palette[, pop]  #assign color
print(files)
print(col)

par(mgp=c(1.5,0.5,0), mar=c(2.5,3,0.5,0.5)+0.1)
plot(1,1,type="n",log="x",xlab="Years ago",ylab="Effective population size",xlim=c(1000,10e6),ylim=c(1e2,6.5e5),xaxt="n",yaxt="n")
axis(side=1, at=c(10,100,1000,10000,100000,1000000,10000000),labels=c("1e+01","1e+02","1e+03","1e+04","1e+05","1e+06","1e+07"))
axis(side=2, at=c(0,200000,400000,600000),labels=c("0","2e+05","4e+05","6e+05"))
abline(v=20000,col="grey",lwd=3)
for(file in files)lines(psmc.result(paste(file, ".1Mb_noSexChr.2+2.psmc", sep=""),i.iteration=25,mu=2.9e-9,s=100,g=2.5),col=col,lwd=2)
legend("topleft",pop,col=col,lty=1,lwd=3,cex=0.75,bg="white")
}

dev.off()
