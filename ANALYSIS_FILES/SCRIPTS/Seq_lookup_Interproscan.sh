#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH -o /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/extract_sig_seq_Interproscan_out_3_16_21
#SBATCH -e /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/extract_sig_seq_Interproscan_error_3_16_21

echo "START $(date)"

# load modules
module load InterProScan/5.44-79.0-foss-2018b

# Set paths needed
O=/data/marine_diseases_lab/erin/OrthoFinder_2020
I=/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/ANALYSIS_FILES/Interproscan

# fetch protein sequences from C. virginica protein sequeces to run through Inteproscan
array1=($(cat $I/hemo_dds_deseq_sig_XP_df_lookup.txt))
for i in ${array1[@]}; do
	sed -n "/${i}/,/^>/p" $O/GCF_002022765.2_C_virginica-3.0_protein.faa | sed '$d' >> $I/hemo_dds_deseq_sig_XP_df_lookup_seq.fa
	echo "done ${i}"
done

# Run InterProScan
interproscan.sh -i $I/hemo_dds_deseq_sig_XP_df_lookup_seq.fa -d $I/ -f GFF3,TSV

# -i is the input data
# -b is the output file base
# -d is the output directory, the output filenames are the same as the input filename
# -f is formats

echo "DONE $(date)"
