#For producing site frequency spectra using Lin's script vcf2sfs.r (https://github.com/shenglin-liu/vcf2sfs)

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#vcftools
#r-base
#java

#Requires installation of Stairway Plot 2: https://github.com/xiaoming-liu/stairway-plot-v2

#==============Subset VCF to include only NEZ individuals and keep only variant sites with no missing sites for each
VCF="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.vcf"
folder="/faststorage/project/Coregonus/Aja/Dragonfly/StairwayPlot/"

mkdir $folder/VCF

sbatch -A Coregonus -t 12:00:00 --wrap\
 "vcftools --vcf ${VCF} --keep ${folder}/NEZ.txt --max-missing 1 --mac 1 --recode --recode-INFO-all --out $folder/VCF/NEZ_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.mac1"
 
#3,807,320 SNPs (NEZ)
#4,000,942 SNPs (NEZ_22-23)

#==============Create SFS (run in R)
setwd("/faststorage/project/Coregonus/Aja/Dragonfly/SFS/")
source("/faststorage/home/anbt/software/VCF2SFS/vcf2sfs.R")   #https://github.com/shenglin-liu/vcf2sfs

## Read VCF file and popmap file"," and create a gt (genotype) object.
gt<-vcf2gt("/faststorage/project/Coregonus/Aja/Dragonfly/StairwayPlot/VCF/NEZ_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.mac1.recode.vcf","/faststorage/project/Coregonus/Aja/Dragonfly/StairwayPlot/popmap_NEZ.txt")

#produce folded 1D-SFS from gt object
folded.sfs<-fold.sfs(gt2sfs.raw(gt,c("NEZ")))
write.table(folded.sfs,"NEZ_folded.obs",quote=F,col.names=F,row.names=F)

#Now create the blueprint file with SFS info, path to software, etc. The easiest way is to modify the example blueprint file that comes with the software.
#Refer to the manual: https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/READMEv2.1.pdf

#=============Create batch file
cd "/faststorage/project/Coregonus/Aja/Dragonfly/StairwayPlot/"

for pop in NEZ NEZ22-23
do
java -cp "/faststorage/home/anbt/software/stairway_plot_v2.1.1/stairway_plot_es/" Stairbuilder ${pop}_fold.blueprint
done

#=============Run batch file
for pop in NEZ NEZ22-23
do
sbatch -A Coregonus -t 12:00:00 --job-name blueprint.sh --wrap\
 "bash ${pop}_fold.blueprint.sh"
done

#=============Plot stairway plot including confidence intervals for each population
pops <- c("NEZ")

cols <- c("#efdf48")
names(cols) <- pops

tiff("Stairway_NEZ.tiff",units="cm", width = 16, height = 15, res=200)

par(mar=c(5,4,2,2)+0.1)
plot(1,1, type="n", xlim=c(5,1e7), ylim=c(5,1e6), log="x", xlab="Years ago", ylab="Effective population size",xaxt='n')
axis(side=1, at=c(10,100,1000,10000,100000,1000000,10000000))
abline(v=20000,col="grey",lwd=3)
for(pop in pops)
{
	a<-read.table(sprintf("%s/%s.final.summary",pop,pop),header=T,sep="\t")
	lines(a$year,a$Ne_median,col=cols[pop],lwd=5)
  lines(a$year,a$Ne_2.5.,col=cols[pop],lty=3)
  lines(a$year,a$Ne_97.5.,col=cols[pop],lty=3)
}
legend("topright",pops,col=cols,lty=1,lwd=3,bty="n")
dev.off()

#=============Plot stairway plot including confidence intervals for each population
pops <- c("NEZ")

cols <- c("#efdf48")
names(cols) <- pops

tiff("Stairway_NEZ22-23.tiff",units="cm", width = 16, height = 15, res=200)

par(mar=c(5,4,2,2)+0.1)
plot(1,1, type="n", xlim=c(5,1e7), ylim=c(5,1e6), log="x", xlab="Years ago", ylab="Effective population size",xaxt='n')
axis(side=1, at=c(10,100,1000,10000,100000,1000000,10000000))
abline(v=20000,col="grey",lwd=3)
for(pop in pops)
{
	a<-read.table(sprintf("%s/%s.final.summary",paste(pop, "22-23", sep=""),pop),header=T,sep="\t")
	lines(a$year,a$Ne_median,col=cols[pop],lwd=5)
  lines(a$year,a$Ne_2.5.,col=cols[pop],lty=3)
  lines(a$year,a$Ne_97.5.,col=cols[pop],lty=3)
}
legend("topright",pops,col=cols,lty=1,lwd=3,bty="n")
dev.off()
