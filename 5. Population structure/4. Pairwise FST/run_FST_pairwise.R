
### Pairwise Fst (R)
library(snpR)
library(ggplot2)
library(ggrepel)

vcf<-read_vcf("/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.vcf",
              sample.meta = "/faststorage/project/Coregonus/Aja/Dragonfly/FST_pairwise/popmap2.csv")
              
vcf <- calc_pairwise_fst(vcf, facets = "POP", boot = 100)
result <- get.snpR.stats(vcf, "POP", "fst") #get results

write.table(result$fst.matrix$POP$fst, "/faststorage/project/Coregonus/Aja/Dragonfly/FST_pairwise/Fst.thinned2.txt", row.names = F) #save pairwise fst values
write.table(result$fst.matrix$POP$p, "/faststorage/project/Coregonus/Aja/Dragonfly/FST_pairwise/pval.thinned2.txt", row.names = F) #save p values
