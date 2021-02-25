#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/Stringtie_output
#SBATCH	-e /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/Stringtie_error

#This script takes bam files from HISAT (processed by SAMtools) and performs StringTie assembly and quantification and converts

echo "START" $(date)

module load StringTie/2.1.1-GCCcore-7.3.0 # new version of Stringtie
module load gffcompare/0.11.5-foss-2018b # new version of gffcompare

# Create variable for each path to trimmed and quality filtered data folder
D=/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/RAW_DATA
CV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/Cvir_Genome_and_Indexes

# assemble transcripts for each sample with the GFF3 annotation file
array1=($(ls $D/*.bam))
for i in ${array1[@]}; do
	stringtie -G $CV/ref_C_virginica-3.0_top_level.gff3 -o ${i}.gtf ${i}
	echo "${i} assembled"
	echo "${i}.gtf" >> $D/Dermo_2020_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CV/ref_C_virginica-3.0_top_level.gff3 -o $D/Dermo_2020_stringtie_merged.gtf $D/Dermo_2020_mergelist.txt
echo "Dermo merged"

#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CV/ref_C_virginica-3.0_top_level.gff3 -G -o $D/Dermo_2020_stringtie_merged $D/Dermo_2020_stringtie_merged.gtf
echo "Dermo gffcompared"

#Re-estimate transcript abundance after merge step
for i in ${array1[@]}; do
		stringtie -A $(echo ${i}|sed "s/\..*//").abd.tab -e -G $D/Dermo_2020_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} Dermo transcript abundance re-estimated"
done

echo "Dermo 2020 Stringtie complete $(date)"


echo "DONE ALL $(date)"
