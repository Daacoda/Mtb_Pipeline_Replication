#MDR-TB Genomics Pipeline (Phase 2 - Alignment, filtering and taxonomic classification)
This phase focuses on removing human and viral reads from metagenomic data and classifying bacterial sequences using Centrifuge for taxonomic profiling. All visualizations were performed using RStudio Server.

#This repository documents the step-by-step replication of Phase 2 from:

#Comparative genomic analysis of multi-drug resistant Mycobacterium tuberculosis clinical isolates from Nigeria (PLOS ONE, 2021)



# 1. Environment Setup
Phase 2 builds upon the tbqc Conda environment created in Phase 1, with additional R and visualization dependencies.

# Activate environment
conda activate tbqc

# bwa (for mapping reads to host/viral genomes)
conda install -c bioconda bwa -y

# seqtk (for subsampling reads)
conda install -c bioconda seqtk -y

# centrifuge (for taxonomic classification)
conda install -c bioconda centrifuge -y

# ncbi-datasets-cli (download bacterial genomes)
conda install -c conda-forge ncbi-datasets-cli -y


#2. Test tools
bwa          #Aignment tool
seqtk        #subsampling tool
Centrifuge   #taxonomic classification
datasets     #download genomes


#3. Tools Versions
bwa v0.7.19-r1273
seqtk Version: 1.5-r133
Centrifuge version 1.0.4
datasets version: 18.7.0


# 4. Subsampling Reads with Seqtk
#Subsample 200,000 paired-end reads
seqtk sample -s42 <(zcat trimmed_reads/*_1_p*.fastq.gz) 200000 | gzip > sub_R1.fastq.gz
seqtk sample -s42 <(zcat trimmed_reads/*_2_p*.fastq.gz) 200000 | gzip > sub_R2.fastq.gz


# 5. Host Genome Removal (Human GRCh38) 
# Download the latest GRCh38.p14 human reference genome from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# Unzip the genome file
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz# 3. Host Genome Removal (Human GRCh38)

#Index the human genome
bwa index GCF_000001405.40_GRCh38.p14_genomic.fna

# Map subsampled reads to human genome
bwa mem GCF_000001405.40_GRCh38.p14_genomic.fna sub_R1.fastq.gz sub_R2.fastq.gz > host_aln.sam

# Convert SAM → BAM
samtools view -S -b host_aln.sam > host_aln.bam

# Alignment stats
samtools flagstat host_aln.bam

# Sort and extract unmapped pairs
samtools sort -n -@2 -o host_aln_namesorted.bam host_aln.bam
samtools view -b -f 12 -F 256 host_aln_namesorted.bam > host_aln_unmapped_pairs.bam
samtools fastq -1 host_cleaned_R1.fastq.gz -2 host_cleaned_R2.fastq.gz -0 /dev/null -s /dev/null -n host_aln_unmapped_pairs.bam


# 6. Viral Decontamination
For hardware constraints, 10 common human viruses were used for practice.

# Run the script
./vir_gen.sh

# Index the 10 viral genomes
bwa index sub_vir10.fasta

# Map reads to viral reference
bwa mem viral_human_sub10/sub_vir10.fasta host_cleaned_R1.fastq.gz host_cleaned_R2.fastq.gz > vir_aln.sam
samtools view -bS vir_aln.sam > vir_aln.bam
samtools flagstat vir_aln.bam

# Extract host+virus-unmapped reads
samtools view -b -f 12 -F 256 vir_aln.bam > vir+host_unmapped.bam
samtools fastq -1 bac_input_R1.fastq.gz -2 bac_input_R2.fastq.gz -0 /dev/null -s /dev/null -n vir+host_unmapped.bam


# 7. Bacterial Classification Using Centrifuge
i chose the the top 25 bacteria that were discovered from the journal due to hardware constraints
 
# Run the bacterial genome download script
./bac_downl.sh

# Extract accession numbers
grep "^>" bacteria_com25/bacteria_25.fna | cut -d' ' -f1 | sed 's/^>//' > bacteria_accessions.txt

# Run the conversion table generator script
./conv_table.sh

# Download and extract taxonomy files
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir taxdump
tar -xvzf taxdump.tar.gz -C taxdump

# Build Centrifuge index
centrifuge-build \
    -p 4 \
    --conversion-table conversion_table.txt \
    --taxonomy-tree taxdump/nodes.dmp \
    --name-table taxdump/names.dmp \
    bacteria_com25/bacteria_25.fna \
    bacteria_index

# Run Centrifuge classification
centrifuge -x bacteria_index \
    -1 bac_input_R1.fastq.gz \
    -2 bac_input_R2.fastq.gz \
    -S centrifuge_bac_class_output.txt \
    --report-file centrifuge_bac_class_report.tsv \
    -p 4


# 8. Alignment Summary Report

# Create alignment summary (from terminal)
echo "Sample,Total_Reads,Mapped_Reads,Mapped_Percent,Properly_paired,Singletons
host_aln,400000,12199,3.05,11058,27
host_aln_unmapped_pairs,387774,0,0.00,0,0
vir_aln,387774,0,0.00,0,0" > alignment_summary.csv



# 9. Visualization in RStudio Server

# RStudio Server Setup
Visualization was performed using RStudio Server running inside Ubuntu (WSL).

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)


# Alignment Summary (Mapped vs Unmapped Reads)
mydata <- read.csv("alignment_summary.csv")

# Data cleaning
str(mydata)
summary(mydata)
head(mydata)

# Add Unmapped Reads column
mydata_long <- mydata %>%
  mutate(Unmapped_Reads = Total_Reads - Mapped_Reads) %>%
  pivot_longer(cols = c(Mapped_Reads, Unmapped_Reads),
               names_to = "Category",
               values_to = "Reads")

# Stacked Bar Chart
ggplot(mydata_long, aes(x = Sample, y = Reads, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Mapped vs Unmapped Reads",
       y = "Number of Reads",
       x = "Sample")

# Percentage Pie Chart
mydata_long_pct <- mydata_long %>%
  group_by(Sample) %>%
  mutate(Percent = Reads / sum(Reads) * 100)

ggplot(mydata_long_pct, aes(x = "", y = Percent, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~Sample) +
  theme_void() +
  labs(title = "Mapped vs Unmapped Reads (%)")



# Bacterial Taxonomic Abundance
class_data <- read_tsv("centrifuge_bac_class_report.tsv")

# Clean and sort
df_class_data <- class_data %>%
  select(name, numReads, abundance) %>%
  arrange(desc(abundance))

# Horizontal Bar Chart
ggplot(df_class_data, aes(x = reorder(name, abundance), y = abundance, fill = name)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Taxonomic Abundance",
    x = "Organism",
    y = "Abundance"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


#Results

The sample is dominated by M.tuberculosis reads (>99% of classified bacterial reads), consistent with the dataset and the selected reference panel.

No viral contamination** was detected against the small 10-virus DB (0 mapped reads); this likely reflects either absence of those viruses or the limited viral index used for this test.

Low taxonomic diversity as expected here because 
	(1) subsampling to 200k pairs. 
	(2) the intentionally limited bacterial reference set (25 genomes) used for replication.

The pipeline successfully removed host and viral reads and produced a reproducible bacterial classification result ready for downstream variant calling or further analyses.


#Limitations 

 Subsampling and a small reference DB reduce sensitivity to low-abundance taxa — results are appropriate for pipeline demonstration and CV/portfolio use but not for publishing large-scale microbiome claims.
 SURPI-style edit-distance sweeps were not performed here due to index scope and compute constraints (not required to demonstrate the pipeline logic).
