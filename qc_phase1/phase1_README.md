MDR-TB Genomics Pipeline (Phase 1 - QC)

This repository documents the step-by-step replication of the pipeline from

##Comparative genomic analysis of multi-drug resistant Mycobacterium tuberculosis clinical isolates from Nigeria (PLOS ONE, 2021)
link.. journals.plos.org/plosone/article?id=10.1371/journal.pone.0258774

## Environment Setup
bash

#1. Check conda is installed
   conda --version

#2. Configure channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict

# 3. Create project environment
     conda create -n tbqc fastqc multiqc trimmomatic bwa samtools bcftools -y

# 4. Activate the environment
    conda activate tbqc

# 5. Test tools
    fastqc -h           #Quality control
    multiqc -h          #Quality control
    trimmomatic -h      #Quality control
    bwa	                #Alignment tool
    samtools            #Alignment file operations
    bcftools            #Variant calling

# 6. Tools Versions
    FastQC v0.12.1
    multiqc, version 1.30
    Trimmomatic v0.39
    conda 25.5.1


## Quality Control Workflow 
The raw Illumina paired-end reads were subjected to quality control using FastQC and summarized with MultiQC. Adapter trimming and low-quality base removal were performed with Trimmomatic. The pipeline ensured that downstream analysis used only high-quality reads.


After activating the environment:

```bash

# 1. Run FastQC on raw reads
  fastqc raw_data/*.fastq.gz -o qc/pre_trim/

# 2. Summarize results with MultiQC
  multiqc qc/pre_trim/ -o qc/pre_trim/

# 3. Trim adapters and low-quality reads
  trimmomatic PE -threads 4 -phred33 \
  raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz \
  qc/sample_R1_paired.fq.gz qc/sample_R1_unpaired.fq.gz \
  qc/sample_R2_paired.fq.gz qc/sample_R2_unpaired.fq.gz \
  ILLUMINACLIP:/path/to/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# 4. Run FastQC again on trimmed reads
  fastqc qc/*_paired.fq.gz -o qc/post_trim/

# 5. Summarize post-trim QC with MultiQC
  multiqc qc/post_trim/ -o qc/post_trim/


#QC Summary (Biological Interpretation)

Based on the quality control analyses of the Mycobacterium tuberculosis SRR31065062 sequencing data:

#Per-base quality (Phred scores)
*Post-trimming reads exhibit median Phred scores ≥ 30, indicating very low base-calling errors (<0.1%).
*This ensures high confidence in nucleotide calls, which is essential for alignment, variant detection, and downstream classification.

#Adapter sequences
*Raw reads contained Illumina adapter contamination, particularly at the 3′ ends.
*Trimming successfully removed these adapters, preventing false alignments and improving the accuracy of downstream analyses.

#Read length distribution
*Trimming removed low-quality bases from read ends, resulting in slightly shorter reads.
*The retained read lengths remain sufficient for genome coverage and robust taxonomic profiling.

#GC content
*Raw reads had a slightly lower GC content (~61%) due to low-quality AT-rich tails.
*Post-trimming reads show GC content ~64%, closer to the expected M. tuberculosis genome (~65%).
*This confirms that trimming enriches for true biological sequences and reduces sequencing artifacts.

#Overall data integrity
*No excessive duplication or abnormal sequence patterns detected after trimming.
*The dataset maintains biological fidelity, reflecting the characteristics of M. tuberculosis without evidence of major contamination.

##All QC and trimming commands are saved in `scripts/qc_commands.sh` for reproducibility.

