#!bin/bash

##############
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Shell script that downloads data from ncbi, performs QC of the reads before and after the trimming procedure. 
#Than script maps trimmed only paired reads to the reference genome, prepares
#sorted and indexed BAM-files suitable for visualization with IGV

#Requirements: sra toolkit, factqc, trimmomatic, bwa mem, samtools 
#This variables should be in the path (or replace them with the path to the particular program)
##############



#######
#Variables to be defined.
#######

#Path to the working directory, contains /Raw_data folder with raw reads files.
PWD='/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/NUMT_masked'
RAW='/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/Raw_data'

echo $PWD
cd $PWD

#Path to the file containing sequencing adapters sequences for trimmomatic uses. Typically in the Trimmomatic-0.36/adapters/XXX.fa
Adapters='All_TruSeq.fa'
trimmomatic='trimmomatic-0.39.jar'
#Path to the reference genome.
Ref_genome=$PWD/Ref_genome/GCA_000001735.2_TAIR10.1_genomic_NUMT_masked.fna




#######
#Quality control and sequencing data preparation.
#######
echo '
######################
Initial quality control is in progress...
######################
'

#Initial quality control
mkdir $PWD/Fastqc_analysis/
mkdir $PWD/Fastqc_analysis/Initial
fastqc -t 20 -o $PWD/Fastqc_analysis/Initial $RAW/*



#######
#Reads trimming
#######
echo '
######################
Reads trimming...
######################
'

mkdir $PWD/Trimmed/
for i in `ls -a $RAW/ | grep 'fastq' | sed -r "s/(.+)_R[1,2]_001\.fastq\.gz/\1/g" | uniq | sort -d`; do
echo $i
java -jar $trimmomatic PE -threads 15 -phred33 $RAW/${i}_R1_001.fastq.gz $RAW/${i}_R2_001.fastq.gz $PWD/Trimmed/${i}_paired_R1.fastq.gz $PWD/Trimmed/${i}_unpaired_R1.fastq.gz $PWD/Trimmed/${i}_paired_R2.fastq.gz $PWD/Trimmed/${i}_unpaired_R2.fastq.gz ILLUMINACLIP:$Adapters:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:30 ; done



#######
#Quality control after the trimming procedure
#######
echo '
######################
Quality control after trimming...
######################
'
mkdir $PWD/Fastqc_analysis/Trimmed/
fastqc -t 20 -o $PWD/Fastqc_analysis/Trimmed/ $PWD/Trimmed/*



#######
#Prepare index for reference genome.
#######
echo '
########################
Reference genome indexing...
########################
'

bwa index $Ref_genome

#######
#Reads mapping, alignment conversion to IGV-compatible format (sorted indexed BAM).
#######
echo '
######################
Reads mapping, SAM files generation...
######################
'

#Reads mapping to the reference genome: make SAM-files
mkdir $PWD/SAM/
for i in `ls -a /home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/Trimmed/ | grep '_paired_' | sed -r "s/(.+)_paired_R[1,2]\.fastq\.gz/\1/g" | uniq | sort -d`; do 
bwa mem -t 15 $Ref_genome /home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/Trimmed/${i}_paired_R1.fastq.gz /home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/Trimmed/${i}_paired_R2.fastq.gz > $PWD/SAM/$i.sam; done

#######
#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
#######
echo '
######################
BAM files preparation...
######################
'
mkdir $PWD/BAM/ChIP/
#Makes BAM-files
for i in `ls -a $PWD/SAM/ChIP/ | grep '.sam' | sed -r "s/(.+).sam/\1/g"`; do 
samtools view -S -b $PWD/SAM/ChIP/${i}.sam > $PWD/BAM/ChIP/${i}.bam ; done

#Sorts BAM-files
echo '
######################
BAM files sorting...
######################
'
mkdir $PWD/BAM_sorted/ChIP/
for i in `ls -a $PWD/BAM/ChIP/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools sort $PWD/BAM/ChIP/${i}.bam -o $PWD/BAM_sorted/ChIP/${i}_sorted.bam ; done

#Converts bam to bed.
echo '
#######################
BAM to bed conversion...
#######################
'

mkdir $PWD/Cov_depth/
mkdir $PWD/Cov_depth/ChIP/

for i in `ls -a $PWD/BAM_sorted/ChIP/ | grep '.bam$' | sed -r "s/(.+)_sorted.bam/\1/g"`; do
samtools depth -a -d 0 $PWD/BAM_sorted/ChIP/${i}_sorted.bam > $PWD/Cov_depth/ChIP/${i}.bed; done


echo '
Script ended its work succesfully!
'

