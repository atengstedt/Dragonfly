#Script for running PCA on dragonfly data

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#r-base

#===========PCA (run in R)
f.vcf<-"/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.vcf"	# Give the absolute pathname of the VCF.
f.popmap<-"/faststorage/project/Coregonus/Aja/Dragonfly/PCA_LD-thinned/popmap2_region.txt"	# Give the absolute pathname of the popmap file.

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
pops<-c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")
colp<-c("#9b5fe0","#16a4d8","#60dbe8","#8bd346","#efdf48","#f9a52c","#d64e12")
names(colp)<-pops
cols<-colp[popmap[,2]]

#PCA biplot
png(file="/faststorage/project/Coregonus/Aja/Dragonfly/PCA_LD-thinned/PCA1-3_region_noSexChr.thinned2.png", units="in", width = 2, height = 4.1, res=300)

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
legend("left",horiz=FALSE,legend=c("NWJ","EJ","WJ","SJ","NEZ","SZ","B"), col=c("#9b5fe0","#16a4d8","#60dbe8","#8bd346","#efdf48","#f9a52c","#d64e12"), bty="n", pch=4,cex=0.8,y.intersp=1.2)  #custom legend

dev.off()

# Plot the eigenvalues.
dev.new() 
barplot(e,main="Eigenvalues")
savePlot("/faststorage/project/Coregonus/Aja/Dragonfly/PCA_LD-thinned/eigenvalues2_region_noSexChr.png","png")
