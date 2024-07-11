#Overview

In this directory, the approach used to sequencing processing and differentially expressed genes analysis


#01.rawdata_to_read_summarization

Commands are detailed in the bash-script:

workflow.sh


#This script contains the flowing processes:

#a. Rawdata quality control and evaluation

1. Using trimmomatic to remove: (1) reads with adapter; (2) reads with N; (3) Low-quality reads(reads with a mass value of Qpred <= 20 accounted for more than 50% of the total read)

2. Using fastQC to evaluate the filetered read quality  

#b. Mapping to reference database:

1. The reference genomes of GA1_B2776 and GA5_B3788 were used to establish the index, respectively (bowtie2-build)

2. The filtered clean reads are mapped to the referance database by Bowtie2

3. Extract mapped reads (using samtools)

#c. Counting reads is perfomed by featurecounts



#02.DesSeq2_and_enrichment


R commands for DesSeq2 and GO enrichment analysis are given in:

"RNA_seq.R"

