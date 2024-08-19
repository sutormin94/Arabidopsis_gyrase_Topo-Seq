#!bin/bash

##############
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Shell script that makes QC of the reads before and after the trimming procedure. 
#Than script maps trimmed only paired reads to the reference genome, prepares
#sorted and indexed BAM-files suitable for visualization with IGV

#Requirements: factqc, trimmomatic, bwa mem, samtools 
#This variables should be in the path (or replace them with the path to the particular program)
##############


#######
#Variables to be defined.
#######

#Path to the working directory, contains /Raw_data folder with raw reads files.
PWD='/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/'
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

#Initial quality control
echo '
######################
Initial quality control is in progress...
######################
'
mkdir $PWD/Fastqc_analysis/
fastqc -t 20 -o $PWD/Fastqc_analysis/ $PWD/Raw_data/*

#######
#Reads mapping, alignment conversion to IGV-compatible format (sorted indexed BAM).
#######

#Reads mapping to the reference genome: make SAM-files
echo '
######################
Reads mapping, SAM files generation...
######################
'
mkdir $PWD/SAM/
for i in `ls -a $PWD/Raw_data/ | grep 'fastq.gz' | sed -r "s/(.+)_R[1,2]\.fastq\.gz/\1/g" | uniq | sort -d`; do 
bwa mem -t 20 $Ref_genome $PWD/Raw_data/${i}_R1.fastq.gz $PWD/Raw_data/${i}_R2.fastq.gz > $PWD/SAM/$i.sam; done

#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
mkdir $PWD/BAM_sorted/
mkdir $PWD/BAM/
#Makes BAM-files
echo '
######################
BAM files preparation...
######################
'
for i in `ls -a $PWD/SAM/ | grep '.sam' | sed -r "s/(.+).sam/\1/g"`; do 
samtools view -S -b $PWD/SAM/${i}.sam > $PWD/BAM/${i}.bam ; done

#Sorts BAM-files
echo '
######################
BAM files sorting...
######################
'
for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools sort $PWD/BAM/${i}.bam $PWD/BAM_sorted/${i}_sorted.bam ; done

#Makes index files
echo '
######################
BAM files indexing...
######################
'
for i in `ls -a $PWD/BAM_sorted/`; do 
samtools index $PWD/BAM_sorted/${i} ; done

echo '
Script ended its work succesfully!
'

