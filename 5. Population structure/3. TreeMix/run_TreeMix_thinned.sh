#TreeMix analysis using scripts from https://github.com/carolindahms/TreeMix

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#bcftools
#r-base

#Requires installation of a number of R packages. Check https://github.com/carolindahms/TreeMix for original scripts and instructions.

cp "/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.vcf" "/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.vcf"

#-----------make clust file
input="/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.clust"

bcftools query -l ${input} | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[2]}' > ${output}
#manually add a third column with population name

###activate environment
conda activate treemix

#----------convert VCF to TreeMix file
vcf="/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.vcf"
clust="/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.clust"

"/faststorage/project/Coregonus/Aja/Coregonus/WGS/treemix/vcf2treemix.sh" ${vcf} ${clust}

####============Run TreeMix
####following tutorial on https://github.com/carolindahms/TreeMix

#Build consensus tree with multiple migration events
input="/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.treemix.frq.gz"
output="/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/noRoot2_k100_bs500_m1-10_rep10/"

mkdir ${output}
cd ${output}

sbatch -A Coregonus -t 12:00:00 --mem 8G -c 2 --job-name TreeMix_step1 --wrap \
 "sh ../Step1_TreeMix.sh ${input} 2 100 noRoot 500 consense mytree 1 10 10"

#Test migration edges with OptM
#run part A of Step2+4_TreeMix.R
#conda environment TreeMix

#Final runs with optimum number of migration edges
input="/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.treemix.frq.gz"

cd "/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/noRoot2_k100_bs500_m1-10_rep10/"

#1 migration event
#sbatch -A Coregonus -t 48:00:00 --mem 25G -c 8 --job-name TreeMix_step3 --wrap \
# "sh ../Step3_TreeMix.sh ${input} 8 100 noRoot 500 1 mytree 10 mytree2_constree.newick consense"
#mv "/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/noRoot2_k100_bs500_m1-10_rep10/final_runs/" "/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/noRoot2_k100_bs500_m1-10_rep10/final_runs_m1/"

#2 migration events
sbatch -A Coregonus -t 48:00:00 --mem 8G -c 2 --job-name TreeMix_step3 --wrap \
 "sh ../Step3_TreeMix.sh ${input} 2 100 noRoot 500 2 mytree 10 mytree2_constree.newick consense"

#Tree visualization + Migration stats and support
#run part B of Step2+4_TreeMix.R
#conda environment TreeMix


#calculate variance explained (R)
Rscript "/faststorage/project/Coregonus/Aja/Coregonus/WGS/treemix/VarianceExplained.R" mytree_1m_finalrun_1

#1 migration event
#Standard error for all entries in the covariance matrix estimated from the data 0.000465826510204082
#Variance of relatedness between populations explained by the model      0.980932712095605

Rscript "/faststorage/project/Coregonus/Aja/Coregonus/WGS/treemix/VarianceExplained.R" mytree_2m_finalrun_1

#2 migration events
#Standard error for all entries in the covariance matrix estimated from the data 0.000465826510204082
#Variance of relatedness between populations explained by the model      0.991050664326613

             
