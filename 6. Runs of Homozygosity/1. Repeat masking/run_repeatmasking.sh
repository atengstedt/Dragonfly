### Perform repeat masking on draft assembly

### Using a combination of protocols
### https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
### https://star-protocols.cell.com/protocols/1799#fig2
### https://github.com/jeppebayer/EcoGenetics/tree/master/scripts/05_repeat_masking/modules
### https://www.bioinformatics.uni-muenster.de/publication_data/P.californicus_annotation/repeat_masking.hbi

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

### 0. Prepare environment

#download environment yml file from Jeppes github page 
# https://github.com/jeppebayer/EcoGenetics/tree/master/scripts/05_repeat_masking/env

#make environment
conda env create -f repeatmasking.yml

#activate environment
conda activate repeatmasking

### 1. Run RepeatModeler on Assembly

# build new RepeatModeler BLAST database
folder="/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/library/RepeatModelerDB/"
ref="/faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"

cd $folder

BuildDatabase -name AesVir_Aeshna-viridis -engine ncbi $ref

#now run RepeatModeler
folder="/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/"
cd $folder/repeatmodeler

sbatch -A Coregonus -t 48:00:00 --job-name RepeatModeler --mem 100G -c 32 --wrap\
 "RepeatModeler -pa 32 -engine ncbi -database $folder/library/RepeatModelerDB/AesVir_Aeshna-viridis -LTRStruct"
 
#move results from database folder
mv "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/library/RepeatModelerDB/AesVir_Aeshna-viridis-families.fa" "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmodeler/"
mv "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/library/RepeatModelerDB/AesVir_Aeshna-viridis-families.stk" "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmodeler/"
mv "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/library/RepeatModelerDB/AesVir_Aeshna-viridis-rmod.log" "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmodeler/"

#add species code to output fasta file
cat ./AesVir_Aeshna-viridis-families.fa | \
seqkit fx2tab | \
awk '{ print "EupAur_"$0 }' | \
seqkit tab2fx \
> AesVir_Aeshna-viridis-families.prefix.fa

#split fasta file into classified and unclassified elements
cat ./AesVir_Aeshna-viridis-families.prefix.fa | \
seqkit fx2tab | \
rg -j 10 -v "Unknown" | \
seqkit tab2fx \
> AesVir_Aeshna-viridis-families.prefix.fa.known

cat ./AesVir_Aeshna-viridis-families.prefix.fa | \
seqkit fx2tab | \
rg -j 10 "Unknown" | \
seqkit tab2fx \
> AesVir_Aeshna-viridis-families.prefix.fa.unknown

# quantify number of classified elements
grep -c ">" AesVir_Aeshna-viridis-families.prefix.fa.known #934
# quantify number of unknown elements
grep -c ">" AesVir_Aeshna-viridis-families.prefix.fa.unknown #1872

### 2. Run RepeatMasker on Assembly using RepBase DB
sbatch -A Coregonus -t 12:00:00 --job-name RepeatMasker --mem 50G -c 16 --wrap\
 "RepeatMasker -e rmblast -pa 16 -dir /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepBaseRun -xsmall -lib /faststorage/project/Coregonus/Aja/Butterflies/assembly/repeatMasking/library/RepBase/Arthropoda.rep /faststorage/project/Coregonus/Aja/Dragonfly/ref/dragonfly.contig.fa"
 
#make GFF 3 of repeat regions from RepeatMasker output file
"/faststorage/project/Coregonus/Aja/Butterflies/assembly/repeatMasking/rm2gff3.sh" "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepBaseRun/dragonfly.contig.fa.out" > "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepBaseRun/AesVir_genomic.repeats_giri.gff"

### 3. Run RepeatMasker on Assembly using RepeatModeler data
sbatch -A Coregonus -t 12:00:00 --job-name RepeatMasker --mem 50G -c 16 --wrap\
 "RepeatMasker -e rmblast -pa 8 -dir /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepModRun/ -xsmall -lib /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmodeler/AesVir_Aeshna-viridis-families.prefix.fa /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepBaseRun/dragonfly.contig.fa.masked"

#make GFF 3 of repeat regions from RepeatMasker output file
"/faststorage/project/Coregonus/Aja/Butterflies/assembly/repeatMasking/rm2gff3.sh" "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepModRun/dragonfly.contig.fa.masked.out" > "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepModRun/AesVir_genomic.repeats_giri.gff"

### 4. Combine result of 2 and 3 and Process them to create new files

# Combine full RepeatMasker result files
cat \
"/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepBaseRun/dragonfly.contig.fa.cat.gz" \
"/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepModRun/dragonfly.contig.fa.masked.cat.gz" \
> "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/AesVir_genomic.full_mask.cat.gz"

# Combine RepeatMasker tabular files
cat \
"/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepBaseRun/dragonfly.contig.fa.out" \
<(tail -n +4 "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/RepModRun/dragonfly.contig.fa.masked.out") \
> "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/AesVir_genomic.full_mask.out"

sbatch -A Coregonus -t 12:00:00 --job-name ProcessRepeats --mem 32G --wrap\
 "ProcessRepeats -lib /faststorage/project/Coregonus/Aja/Butterflies/assembly/repeatMasking/library/RepBase/Arthropoda.rep /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/AesVir_genomic.full_mask.cat.gz"

# Make GFF 3 of repeat regions from RepeatMasker output file
"/faststorage/project/Coregonus/Aja/Butterflies/assembly/repeatMasking/rm2gff3.sh" "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/repeatmasker/AesVir_genomic.full_mask.out" > "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/AesVir_genomic.repeats.gff"

# Make BED3 file of repeat regions
awk \
    -F "\t" \
    'BEGIN {OFS = "\t"}
    {if ($0 ~ /^[^#]/)
        {print $1, ($4 - 1), $5}
    }' \
    /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/AesVir_genomic.repeats.gff \
    > /faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/AesVir_genomic.repeats.bed

##merge book-ended features in bed file (increases speed of subsequent filtering process)
bedtools merge -i "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/AesVir_genomic.repeats.bed" > "/faststorage/project/Coregonus/Aja/Dragonfly/repeatmasking/AesVir_genomic.repeats.merged.bed"
