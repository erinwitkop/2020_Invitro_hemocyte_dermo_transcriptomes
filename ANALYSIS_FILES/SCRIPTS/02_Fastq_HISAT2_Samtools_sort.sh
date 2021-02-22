#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/HISAT_out
#SBATCH	-e /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/HISAT_error_

echo "START" $(date)

module load HISAT2/2.1.0-foss-2018b
module load SAMtools/1.9-foss-2018b

# Create variable for each path to trimmed and quality filtered data folder
D=/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/RAW_DATA
CV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/Cvir_Genome_and_Indexes

############## BUILDING HISAT2 INDEXES ###################
#Index will be made with reference genome and no annotation file (allowing for novel transcript discovery)
	# create new directory for the HISAT index called genome, and put the genome inside it
	# copy all reads files into this directory as well to ensure easy access by commands

# C. virginica genome index - already built
#Build HISAT index with cvir_edited (this file has extra spaces in header removed so that genome and annotation don't conflict)
#hisat2-build -f $CV/cvir_edited.fa $CV/cvir_edited_index
  # -f indicates that the reference input files are FASTA files

############# USE HISAT TO ALIGN RNA-READS TO GENOME ##############
# reads are paired end

# C_vir_Probiotic_SRA_ID
array1=($(ls $D/*_1.fastq.gz.clean.trim.filter.gz))
for i in ${array1[@]}; do
  # outputs a single bam file
	hisat2 --dta -x $CV/cvir_edited_index  -1 ${i} -2 $(echo ${i}|sed s/_1/_2/) -S ${i}.sam
	echo "HISAT2 PE ${i} $(date)"
  #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
  samtools sort ${i}.sam > ${i}.bam
  #Get bam file statistics for percentage aligned with flagstat
  samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
 echo "${i} sorted bam done"
done

echo "HISAT DONE $(date)"

# -x before the index
# -U before the file to be aligned
# --dta : Report alignments tailored for transcript assemblers including StringTie.
     # With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites.
     #This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

#reference: Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown
#https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
#http://www.htslib.org/doc/samtools.html
