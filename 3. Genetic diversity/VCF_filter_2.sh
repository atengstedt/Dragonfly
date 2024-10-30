#!/bin/bash
#SBATCH --time=100:00:00

ref="/faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"
from="/faststorage/project/Coregonus/Aja/Dragonfly/mem"
to="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/for_HO_2"

bcftools mpileup -Ou -r $1 -q 30 -Q 25 -a FORMAT/DP,FORMAT/AD -f $ref ${from}/*_purged.bam\
 | bcftools call -f GQ -m\
 | bcftools filter -e'INFO/DP<600 || INFO/DP>1900 || INDEL==1'\
 > ${to}/${1}.vcf
