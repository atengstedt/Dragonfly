##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
########################## 

#Requires the following packages
#bcftools
#vcftools
#plink

#Requires installation of GONE: https://github.com/esrud/GONE

#fix VCF header (list only the 100 scaffolds that are present in VCF)
sbatch -A Coregonus -t 12:00:00 --wrap\
 "bcftools reheader --fai /faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.GONE.fa.fai /faststorage/project/Coregonus/Aja/Dragonfly/GONE/VCF/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.vcf -o /faststorage/project/Coregonus/Aja/Dragonfly/GONE/VCF/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.fai.vcf"
 
#=========Rename scaffolds and annotate SNP IDs in VCF file
#scaffolds need to be numbered (ONLY numbers!) and occur in correct order with no missing scaffolds/numbers
folder="/faststorage/project/Coregonus/Aja/Dragonfly/GONE/VCF/"

sbatch -A Coregonus -t 12:00:00 --mem 8G --wrap\
 "bcftools annotate --rename-chrs /faststorage/project/Coregonus/Aja/Dragonfly/plink/chr_conv.txt ${folder}/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.fai.vcf --set-id '%CHROM\_%POS' -O v -o ${folder}/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.fai.ann.vcf"

#=========split vcf by population and keep only variant sites with no missing sites for each
folder="/faststorage/project/Coregonus/Aja/Dragonfly/GONE/"

for population in NEZALL B2018 NEZ2018 NEZ2023 NWJ2021 SJ2021
do

mkdir $folder/$population

sbatch -A Coregonus -t 12:00:00 --wrap\
 "vcftools --vcf ${folder}/VCF/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.fai.ann.vcf --keep $folder/$population.txt --max-missing 1 --mac 1 --recode --recode-INFO-all --out $folder/${population}/${population}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1"
done

#=========convert to ped/map format
ID_file="/faststorage/project/Coregonus/Aja/Dragonfly/plink/FID.txt"
folder="/faststorage/project/Coregonus/Aja/Dragonfly/GONE/"

for population in NEZALL B2018 NEZ2018 NEZ2023 NWJ2021 SJ2021
do
sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name conv_pedmap --wrap\
 "plink --vcf ${folder}/${population}/${population}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1.recode.vcf --const-fid --allow-extra-chr --update-ids $ID_file --recode --out ${folder}/${population}/${population}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1"
done

#remove prefix from scaffold names
for population in NEZALL B2018 NEZ2018 NEZ2023 NWJ2021 SJ2021
do
sed -i 's/Contig//g' /faststorage/project/Coregonus/Aja/Dragonfly/GONE/${population}/${population}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1.map
done

#remove unnecessary files
for population in NEZALL B2018 NEZ2018 NEZ2023 NWJ2021 SJ2021
do
rm /faststorage/project/Coregonus/Aja/Dragonfly/GONE/${population}/*.log
rm /faststorage/project/Coregonus/Aja/Dragonfly/GONE/${population}/*.nosex
done

#=========run GONE for each population
programme="/faststorage/home/anbt/software/GONE-master/Linux/"
folder="/faststorage/project/Coregonus/Aja/Dragonfly/GONE"

for population in NEZALL B2018 NEZ2018 NEZ2023 NWJ2021 SJ2021
do
mkdir $folder/$population/GONE_hc0.05_7.37cM
cd $folder/$population/GONE_hc0.05_7.37cM

cp -r $programme/* .
cp $folder/$population/${population}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1.map .
cp $folder/$population/${population}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1.ped .

sbatch -A Coregonus -t 12:00:00 --mem 8G -c 16 --job-name GONE --wrap\
 "bash script_GONE.sh ${population}_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1"
done

#remove temporary files
for pop in NEZALL B2018 NEZ2018 NEZ2023 NWJ2021 SJ2021
do
rm -r /faststorage/project/Coregonus/Aja/Dragonfly/GONE/${pop}/GONE*/PROGRAMMES
rm -r /faststorage/project/Coregonus/Aja/Dragonfly/GONE/${pop}/GONE*/TEMPORARY_FILES
rm /faststorage/project/Coregonus/Aja/Dragonfly/GONE/${pop}/GONE*/script_GONE.sh
done
