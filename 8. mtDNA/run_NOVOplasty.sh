#de novo assembly of mitochondrial genomes using NOVOplasty

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#Requires the following packages
#novoplasty

for file in config1 config2 config3
do
sbatch -A Coregonus -t 12:00:00 --mem 300G --wrap\
 "NOVOPlasty4.3.5.pl -c /faststorage/project/Coregonus/Aja/Dragonfly/NOVOplasty/${file}.txt"
done


#Concatenate results
# Create or empty the output file
> mtDNA_sequences.fasta

# Loop through each .fasta file in the directory
for file in Circularized_assembly_*.fasta; do
    # Extract the unique identifier from the filename
    identifier=$(echo "$file" | grep -oP '(?<=NHMA-ENT-)\d+')
    
    # Replace ">Contig1" with the identifier and append to the output file
    sed "s/>Contig1/>NHMA-ENT-${identifier}/" "$file" >> concatenated_sequences.fasta
done
