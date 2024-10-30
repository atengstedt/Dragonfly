#Script for running preliminary PCA on dragonfly data

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#--------Filter for missing data in preparation for preliminary PCA
folder="/faststorage/project/Coregonus/Aja/Dragonfly/"

sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name vcftools --wrap\
 "vcftools --gzvcf ${folder}/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.vcf.gz --max-missing 1 --recode --recode-INFO-all --out ${folder}/PCA_prelim/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.miss1"

#create random subset of VCF file to reduce computation time
sbatch -A Eels -t 1:00:00 --job-name subset --wrap \
 "awk '/^#/{print;next}rand()<0.02{print}' /faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.miss1.recode.vcf > /faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.miss1.subset0.05.vcf"
  
#line count (how many variants in subset file?)
file="/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.miss1.subset0.05.vcf"

sbatch -A Eels -t 1:00:00 --job-name line_count --wrap \
 "bcftools view -H ${file} | wc -l"

### PCA (R)
f.vcf<-"/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.miss1.subset0.05.vcf"	# Give the absolute pathname of the VCF subset.
f.popmap<-"/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/popmap_region.txt"	# Give the absolute pathname of the popmap file.
#popmap file contains individual names in column 1 (in the same order as they appear in the VCF!) and population name in column 2.

#-----------Read the input files.
a<-read.table(f.vcf,sep="\t",stringsAsFactors=F)
popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=F)

#-----------Parse the VCF file.
vcf<-as.matrix(a[,-c(1:9)])
nrow.vcf<-nrow(vcf)
ncol.vcf<-ncol(vcf)
chrom1<-as.integer(substring(vcf,1,1))
chrom2<-as.integer(substring(vcf,3,3))
chrom<-matrix(chrom1+chrom2,nrow.vcf,ncol.vcf)
colnames(chrom)<-popmap[,1]
rm(vcf,chrom1,chrom2)

#-----------Run PCA.
b<-prcomp(t(na.omit(chrom)))
e<-b$sdev^2
indis<-rownames(b$x)

#-----------Function for generating PCA biplot.
plot.pca<-function(pca,x,y,...)
{
	e<-pca$sdev^2
	plot(pca$x[,c(x,y)],
		xlab=sprintf("PC%d (%.1f%%)",x,e[x]/sum(e)*100),
		ylab=sprintf("PC%d (%.1f%%)",y,e[y]/sum(e)*100),...)
}

#-----------Make the plot colorful.
pops<-c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")  #pop names as they appear in the popmap file
colp<-c("#9b5fe0","#16a4d8","#60dbe8","#8bd346","#efdf48","#f9a52c","#d64e12")  #assign pop colours
names(colp)<-pops
cols<-colp[popmap[,2]]

#PCA biplot
png(file="/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/PCA_region.png", units="in", width = 2, height = 4.1, res=300)

layout(matrix(c(4,4,4,1,2,3), 3, 2, byrow=FALSE),
        heights=c(1,1), widths=c(0.4,1))        
        
par(mar=c(2.5,2.3,0.5,0.2)+0.1, mgp=c(1,0.1,0), xpd=TRUE,cex=0.6)
plot.pca(b,1,2,type="n",tck=-0.01)
text(b$x[,c(1,2)],"x",col=cols,xpd=NA,cex=1.5)

par(mar=c(2.5,2.3,0.5,0.2)+0.1, mgp=c(1,0.1,0), xpd=TRUE,cex=0.6)
plot.pca(b,1,3,type="n",tck=-0.01)
text(b$x[,c(1,3)],"x",col=cols,xpd=NA,cex=1.5)

par(mar=c(2.5,2.3,0.5,0.2)+0.1, mgp=c(1,0.1,0), xpd=TRUE,cex=0.6)
plot.pca(b,2,3,type="n",tck=-0.01)
text(b$x[,c(2,3)],"x",col=cols,xpd=NA,cex=1.5)

par(mar=c(0,0,0,0)+0.1)
plot.new()

#legend("left",horiz=FALSE,legend=pops, col=colp, bty="n",pch=19,cex=0.8,y.intersp=1.2)  #automatic legend
legend("left",horiz=FALSE,legend=c("NWJ","EJ","WJ","SJ","NEZ","SZ","B"), col=c("#9b5fe0","#16a4d8","#60dbe8","#8bd346","#efdf48","#f9a52c","#d64e12"), bty="n", pch=4,cex=0.8,y.intersp=1.2)  #legend in custom order

dev.off()

#plot the eigenvalues.
dev.new() 
barplot(e,main="Eigenvalues")
savePlot("/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/eigenvalues.png","png")

#Population structure in preliminary PCA (PC axis 1) suggests that there are sex chromosomes/scaffolds in the data set
#Print loadings from PC1
write.table(b$x[,"PC1"], "/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/PC1.txt",col.names=F,quote=F,sep="\t")

#use this info to create a new popmap file dividing individuals into two groups based on sex (popmap_sex.txt)



###PCA for each scaffold using the popmap where individuals are grouped by sex

#split VCF
input="/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.miss1.recode.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/scaffolds"

mkdir $output

cat scaffold_names.100kb.txt | while read scaffold
do
sbatch -A Eels -t 12:00:00 --job-name subset --wrap \
 "vcftools --vcf $input --chr $scaffold --recode --recode-INFO-all --out ${output}/${scaffold}_filtered"
done

#run PCA for each scaffold
input="/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/scaffolds"
output="/faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/per_scaffold/"

mkdir $output

cat scaffold_names.100kb.txt | while read scaffold
do
sbatch -A Eels -t 12:00:00 --job-name plot_PCA --wrap \
 "Rscript /faststorage/project/Coregonus/Aja/Dragonfly/PCA_prelim/PCA.r ${input}/${scaffold}_filtered.recode.vcf ${output}/${scaffold}.PCA.png"
done

#Check the resulting individual plots to identify scaffolds where individuals are clearly split by sex on either PC1 or PC2
