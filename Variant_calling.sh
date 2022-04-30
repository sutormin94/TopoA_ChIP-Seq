#!/bin/bash

#Prepare folders for data processing.
PWD=/home/cls01/Data_Spongy_metagenomes/E_coli_BW25113_topA_mutants_set3/

mkdir $PWD/Initial_QC
Init_QC=$PWD/Initial_QC

mkdir $PWD/Trimmed_data
Trim_data=$PWD/Trimmed_data

mkdir $PWD/Trimmed_QC
Trim_QC=$PWD/Trimmed_QC

Genome_data=$PWD/Genome/

mkdir $PWD/SAM
SAM_data=$PWD/SAM

mkdir $PWD/SAM_mapped/
SAM_map=$PWD/SAM_mapped

mkdir $PWD/SAM_mapped_header
SAM_map_head=$PWD/SAM_mapped_header

mkdir $PWD/BAM
BAM_data=$PWD/BAM

mkdir $PWD/BCF
BCF_data=$PWD/BCF

mkdir $PWD/VCF
VCF_data=$PWD/VCF

mkdir $PWD/VCF_called
VCF_data_call=$PWD/VCF_called

mkdir $PWD/Cov_depth
Cov_depth_data=$PWD/Cov_depth

mkdir $PWD/BAM_regions
BAM_extracts=$PWD/BAM_regions



echo '
###########
#Data initial QC
###########
'

for i in `ls -a $PWD/Raw_data/`
do
	echo $i        
	fastqc -t 30 -o $Init_QC $PWD/Raw_data/$i
done


echo '
###########
#Data trimming
###########
'

for i in `ls -a $PWD/Raw_data/ | sed -r "s/(.+)_R[1,2]_001.fastq.gz/\1/" | uniq | sort -d`
do
        echo $i
        java -jar /home/cls01/Prog/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 30 -phred33 $PWD/Raw_data/${i}_R1_001.fastq.gz $PWD/Raw_data/${i}_R2_001.fastq.gz  $Trim_data/paired_${i}_R1_001.fastq.gz  $Trim_data/unpaired_${i}_R1_001.fastq.gz  $Trim_data/paired_${i}_R2_001.fastq.gz  $Trim_data/unpaired_${i}_R2_001.fastq.gz ILLUMINACLIP:/home/cls01/Prog/Trimmomatic-0.38/adapters/All_TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done


echo '
###########
#Data QC after trimmomatic
###########
'
for i in `ls -a $Trim_data/`
do
	echo $i
        fastqc -t 30 -o $Trim_QC $Trim_data/$i
done


echo '
###########
#Prepare reference genome, make index.
###########
'

bwa index $Genome_data/BW25113_NZ_CP009273.fasta


echo '
###########
#Reads mapping to the reference genome.
###########
'

for i in `ls -a $Trim_data/ | grep '^paired' | sed -r "s/paired_(.+)_R[1,2]_001\.fastq\.gz/\1/g" | uniq | sort -d`; do 
	echo $i
	bwa mem -t 30 $Genome_data/BW25113_NZ_CP009273.fasta $Trim_data/paired_${i}_R1_001.fastq.gz $Trim_data/paired_${i}_R2_001.fastq.gz > $SAM_data/$i.sam; 
done


echo '
###########
#Keep only mapped reads.
###########
'

for i in `ls -a $SAM_data/ | sed -r "s/(.*).sam/\1/g" | uniq | sort -d`;
do
	echo $i
	samtools view -S -F 4 $SAM_data/${i}.sam > $SAM_map/${i}_mapped.sam
done


echo ' 
###########
#Add headers to sam files with mapped reads only.
###########
'

for i in `ls -a $SAM_data | sed -r "s/(.*).sam/\1/g" | uniq | sort -d`;
do
	echo $i	
	head -2 $SAM_data/${i}.sam | cat - $SAM_map/${i}_mapped.sam > $SAM_map_head/${i}_mapped_header.sam
done



echo '
###########
#Convert SAM to BAM, sort BAM.
###########
'

for i in `ls -a $SAM_data/ | grep ".sam" | sed -r "s/(.*).sam/\1/g" | uniq | sort -d`;
do
	echo $i
	samtools view -S $SAM_data/${i}.sam -b > $BAM_data/${i}.bam
	samtools sort -@ 24 $BAM_data/${i}.bam -o $BAM_data/${i}.sorted.bam
done


echo '
###########
#BAM to bed conversion
###########
'

for i in `ls -a $BAM_data/ | grep '.sorted.bam' | sed -r "s/(.+).sorted.bam/\1/g"`; 
do
	echo $i
	samtools depth -d 0 -a $BAM_data/${i}.sorted.bam > $Cov_depth_data/${i}.bed;
done



echo '
###########
#Prepare VCF from sorted BAM.
###########
'

for i in `ls -a $BAM_data/ | grep ".sorted.bam" | sed -r "s/(.*).sorted.bam/\1/g" | uniq | sort -d`;
do
	echo $i
	bcftools mpileup -f $Genome_data/BW25113_NZ_CP009273.fasta $BAM_data/${i}.sorted.bam -Oz -o $VCF_data/${i}.vcf.gz
done


echo '
###########
#Call variants from VCF.
###########
'

for i in `ls -a $VCF_data/ | grep ".vcf.gz" | sed -r "s/(.*).vcf.gz/\1/g" | uniq | sort -d`;
do
	echo $i
	bcftools call -cv -Ov $VCF_data/${i}.vcf.gz > $VCF_data_call/${i}.vcf
done


