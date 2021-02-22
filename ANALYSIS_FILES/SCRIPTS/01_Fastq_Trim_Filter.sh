#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH	-o /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/trim_out_2_22_21
#SBATCH	-e /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/SCRIPTS/SCRIPT_out_error_files/trim_error_2_22_21
#SBATCH	--mail-user=erin_roberts@my.uri.edu

echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

# Create variable for each path to raw data folder
D=/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/RAW_DATA

##### Paired End Read Trimming + Filtering  ######
# Create array variables

array1=($(ls $D/*R1_001.fastq.gz))
for i in ${array1[@]}; do
	bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/R1/R2/) out2=$(echo ${i}|sed s/R1/R2/).clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "adapter trimming ${i} $(date)"
  #Quality trimmming, of both the left and the right sides to get rid of reads that are less than quality 20
	bbduk.sh in1=${i}.clean out1=${i}.clean.trim in2=$(echo ${i}|sed s/R1/R2/).clean out2=$(echo ${i}|sed s/R1/R2/).clean.trim qtrim=rl trimq=20
	echo "quality trimming ${i} $(date)"
	#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/R1/R2/).clean.trim out2=$(echo ${i}|sed s/R1/R2/).clean.trim.filter maq=10
	rm ${i}.clean
  rm $(echo ${i}|sed s/R1/R2/).clean
  rm ${i}.clean.trim
  rm $(echo ${i}|sed s/R1/R2/).clean.trim
	echo "quality filtering ${i} $(date)"
done

	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge
	#which does not require known adapter sequences)
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)

#Histogram generation, only generating for one of the pair (assuming that similar stats will be present).
#All histogram output contents are combined into one file
for i in ${array1[@]}; do
  	 bbduk.sh in1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/R1/R2/).clean.trim.filter  bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
     echo ${i} > ${i}.hist.all
     echo "bhist" >> ${i}.hist.all
		 cat ${i}.b.hist >> ${i}.hist.all
		 echo "qhist" >> ${i}.hist.all
     cat ${i}.q.hist >> ${i}.hist.all
     echo "gchist" >> ${i}.hist.all
     cat ${i}.gc.hist >> ${i}.hist.all
     echo "lhist" >> ${i}.hist.all
     cat ${i}.l.hist >> ${i}.hist.all
	   echo "histogram DONE $(date)"
		 rm ${i}.*.hist
	   gzip ${i}.clean.trim.filter
		 gzip $(echo ${i}|sed s/R1/R2/).clean.trim.filter
done
		#lhist = output a read length histogram
        #qhist = per base average quality
        #bhist = output a per-base composition histogram
        #gchist = output a gc content histogram


echo "full trim complete $(date)"
