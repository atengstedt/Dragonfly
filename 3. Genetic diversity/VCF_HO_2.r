input<-commandArgs(trailingOnly=T)
f.in<-input[1]

vcf<-read.table(f.in,sep="\t",stringsAsFactors=F)
vcf<-as.matrix(vcf[,-c(1:9)])
nrow.vcf<-nrow(vcf)
ncol.vcf<-ncol(vcf)
chrom1<-matrix(as.integer(substring(vcf,1,1)),nrow.vcf,ncol.vcf)
chrom2<-matrix(as.integer(substring(vcf,3,3)),nrow.vcf,ncol.vcf)

cat(sub("VCF/for_HO_2/(.*).vcf","\\1",f.in),nrow.vcf,colMeans(chrom1!=chrom2,na.rm=T),sep="\t")
cat("\n")
