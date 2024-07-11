#!/bin/bash

declare -a names=(
    Co_20X_1
    Co_20X_2
    Co_20X_3
    Co_1000X_1
    Co_1000X_2
    Co_1000X_3
)

for i in ${names[@]};do
	R1=$i"_1.fq.gz"
    R2=$i"_2.fq.gz"
    PREFIX=$i

#1. trimmomatic and fastQC were used for quality control and evaluation

	java -jar trimmomatic-0.39.jar PE -phred33 $R1 $R2  $PREFIX"_1.clean.fq.gz"  $PREFIX"_1.unpaired.fq.gz" $PREFIX"_2.clean.fq.gz" $PREFIX"_2.unpaired.fq.gz" ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:1:true  LEADING:25 TRAILING:25  MINLEN:30 -threads 5 
	fastqc -t 16 -o ./ $PREFIX"_1.clean.fq.gz" $PREFIX"_2.clean.fq.gz"



#2. The reference genomes of GA1_B2776 and GA5_B3788 were used to establish the index, respectively, and the filtered sequences were analyzed by Bowtie2

    bowtie2-build -f ref.fa Gilliamella
    bowtie2 -q --phred33 -p 8 -x Gilliamella -1 $PREFIX"_1.clean.fq.gz" -2 $PREFIX"_2.clean.fq.gz" -S  $PREFIX".sam"
    samtools view -@ 16 -b -S $PREFIX".sam" -o $PREFIX".bam"
    rm $PREFIX".sam"

#3. featurecounts was used to summarize the reads in file *.bam

    featureCounts -T 10 -a ref.gtf -o $PREFIX".count" -p -B -C -f -t exon -g gene_id $PREFIX".bam"
done



