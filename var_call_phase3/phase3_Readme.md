# Phase 3 – Variant Calling and Annotation of Mycobacterium tuberculosis

This directory contains scripts, commands, and reports for variant calling, filtering, and annotation of **Mycobacterium tuberculosis** reads extracted after Phase 2 (host and viral removal).

## Environment Setup

Install required tools:

```bash
sudo apt install freebayes tabix bcftools samtools bwa

# Verify installations:
freebayes --version
tabix --version
bcftools --version
samtools --version
bwa

# Reference Preparation
Extract TB reference from bacterial genomes:
samtools faidx aligment_taxclass_phase2/bacteria_com25/bacteria_25.fna NC_000962.3 > H37Rv_ref.fna

Build BWA index:
bwa index H37Rv_ref.fna
samtools faidx H37Rv_ref.fna

# Mapping Reads to TB Reference
bwa mem H37Rv_ref.fna bac_input_R1.fastq.gz bac_input_R2.fastq.gz > tb_mapped.sam
samtools view -bS tb_mapped.sam | samtools sort -o tb_mapped_sorted.bam
samtools index tb_mapped_sorted.bam

# Variant Calling
freebayes -f H37Rv_ref.fna tb_mapped_sorted.bam > tb_variants.vcf
bgzip -c tb_variants.vcf > tb_variants.vcf.gz
tabix -p vcf tb_variants.vcf.gz

# Filtering and Classification
Calculate Allele Frequency and Depth
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AD]\n' tb_variants.vcf.gz \
| awk -F'\t' 'BEGIN{OFS="\t"; print "CHROM","POS","REF","ALT","QUAL","REF_DEPTH","ALT_DEPTH","AF"}
{
  split($6,a,",");
  if(a[1]+a[2]>0) af=a[2]/(a[1]+a[2]);
  else af=0;
  print $1,$2,$3,$4,$5,a[1],a[2],af;
}' > tb_variants_with_AF.tsv

# Classify Variants by Confidence
awk 'BEGIN{OFS="\t"; print "CHROM","POS","REF","ALT","QUAL","REF_DP","ALT_DP","DP","AF","CLASS"}
NR>1{
  dp=$6+$7; alt=$7; ref=$6; af=$8; qual=$5;
  cls="LOWCONF";
  if(dp>=10 && af>=0.90 && alt>=5) cls="HOMO";
  else if(dp>=8 && af>=0.85 && alt>=3) cls="HOMO";
  else if(dp>=6 && alt>=3 && af>=0.05) cls="HETERO";
  else if(dp>=10 && alt>=2 && af>=0.05) cls="HETERO";
  print $1,$2,$3,$4,qual,ref,alt,dp,af,cls
}' tb_variants_with_AF.tsv > tb_variants_classified.tsv

# Quick Counts by Class
awk -F'\t' 'NR>1{c[$10]++} END{for(k in c) print k,c[k]}' tb_variants_classified.tsv

# SnpEff Annotation
Install and Setup
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
mkdir -p data/Mycobacterium_tuberculosis_h37rv
curl -L "https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/references/NC_000962.3.fasta" -o data/Mycobacterium_tuberculosis_h37rv/genome.fa
curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz" -o data/Mycobacterium_tuberculosis_h37rv/genes.gff.gz
gunzip data/Mycobacterium_tuberculosis_h37rv/genes.gff.gz
echo "Mycobacterium_tuberculosis_h37rv.genome : Mycobacterium_tuberculosis_h37rv" >> snpEff.config
java -Xmx2g -jar snpEff.jar build -gff3 -v Mycobacterium_tuberculosis_h37rv

# Prepare VCF and Annotate
gunzip tb_variants.HOMO.vcf.gz
sed -i 's/NC_000962.3/Chromosome/g' tb_variants.HOMO.vcf
bgzip tb_variants.HOMO.vcf
tabix tb_variants.HOMO.vcf.gz
java -Xmx2g -jar snpEff.jar -v Mycobacterium_tuberculosis_h37rv tb_variants.HOMO.vcf.gz > tb_variants.HOMO.ann.vcf

# Filter High and Moderate Impact Variants
java -jar SnpSift.jar filter "ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE'" tb_variants.HOMO.ann.vcf > tb_variants.HOMO.highModerate.vcf

# Extract TSV Report
java -jar SnpSift.jar extractFields tb_variants.HOMO.highModerate.vcf \
CHROM POS REF ALT AF \
GEN[unknown].DP GEN[unknown].AD \
"ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].HGVS_P" \
> report_highModerate.tsv

# Summary Reports
Gene Counts
awk '{print $10}' report_highModerate.tsv | sort | uniq -c | sort -nr | awk '{print $2 "\t" $1}' > top_gene_counts.tsv
printf "GENE\tCOUNT\n" | cat - top_gene_counts.tsv > tmp && mv tmp top_gene_counts.tsv

Impact Summary
awk '{print $9}' report_highModerate.tsv | sort | uniq -c | awk '{print $2 "\t"$1}' > impact_summary.tsv

Variant Effect Counts
awk 'NR>1 {print $8}' report_highModerate.tsv | sort | uniq -c | sort -nr | awk '{print $2"\t"$1}' > effect_counts.tsv

High Impact Variants Only
awk 'NR==1 || $9=="HIGH"' report_highModerate.tsv > report_HIGH.tsv
awk 'NR>1 {print $8}' report_HIGH.tsv | sort | uniq -c | sort -nr > high_effect_counts.tsv

#R plots

#plot of Effect Counts by Type

mut = read_tsv("snpEff/effect_counts.tsv")
colnames(mut) <- c("VariantType", "Count")
ggplot(mut, aes(x = reorder(VariantType, Count), y = Count, fill = VariantType)) +
     geom_col(show.legend = FALSE) +
     coord_flip() +
     labs(title = "Effect Counts by Type", 
          y = "Number of Variants", 
          x = "Effect Type") +
     theme_minimal()

# Notes
All scripts are reproducible and correspond to the PLOS ONE TB variant calling workflow.

Full tb_variants_classified.tsv and report_highModerate.tsv files can be provided as supplementary tables.

For visualization, plots can be generated in R using the TSV outputs.

