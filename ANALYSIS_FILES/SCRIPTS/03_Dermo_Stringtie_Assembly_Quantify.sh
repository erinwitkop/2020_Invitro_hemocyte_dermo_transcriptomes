#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/Stringtie_output_dermo_3_10_21
#SBATCH	-e /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/Stringtie_error_dermo_3_10_21

#This script takes bam files from HISAT (processed by SAMtools) and performs StringTie assembly and quantification and converts

echo "START" $(date)

module load StringTie/2.1.1-GCCcore-7.3.0 # new version of Stringtie
module load gffcompare/0.11.5-foss-2018b # new version of gffcompare

# Create variable for each path to trimmed and quality filtered data folder
PM=/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/ANALYSIS_FILES/Dermo_HISAT_STRINGTIE


# assemble transcripts for each sample with the GFF3 annotation file
array1=($(ls $PM/*.bam))
for i in ${array1[@]}; do
	stringtie -G $PM/GCF_000006405.1_JCVI_PMG_1.0_genomic.gff -o ${i}.gtf ${i}
	echo "${i} assembled"
	echo "${i}.gtf" >> $PM/Dermo_2020_PERK_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $PM/GCF_000006405.1_JCVI_PMG_1.0_genomic.gff -o $PM/Dermo_2020_PERK_stringtie_merged.gtf $PM/Dermo_2020_PERK_mergelist.txt
echo "Dermo PERK merged"

#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $PM/ref_C_virginica-3.0_top_level.gff3 -G -o $PM/Dermo_2020_PERK_stringtie_merged $PM/Dermo_2020_PERK_stringtie_merged.gtf
echo "Dermo PERK gffcompared"

#Re-estimate transcript abundance after merge step
for i in ${array1[@]}; do
		stringtie -A $(echo ${i}|sed "s/\..*//").abd.tab -e -G $PM/Dermo_2020_PERK_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} Dermo PERK transcript abundance re-estimated"
done

echo "Dermo PERK 2020 Stringtie complete $(date)"


echo "DONE ALL $(date)"
