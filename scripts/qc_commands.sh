#!/bin/bash
# Run pre-trim QC
fastqc raw_reads/*.fastq.gz -o qc/pre_trim/

# MultiQC summary
multiqc qc/pre_trim/ -o qc/pre_trim/

# Trimming
trimmomatic PE -threads 4 -phred33 \
raw_reads/sample_R1.fastq.gz raw_reads/sample_R2.fastq.gz \
trimmed_reads/sample_R1_paired.fq.gz trimmed_reads/sample_R1_unpaired.fq.gz \
trimmed_reads/sample_R2_paired.fq.gz trimmed_reads/sample_R2_unpaired.fq.gz \
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Post-trim QC
fastqc trimmed_reads/*_paired.fq.gz -o qc/post_trim/
multiqc qc/post_trim/ -o qc/post_trim/
