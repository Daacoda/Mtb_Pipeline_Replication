#install freebayes tabix bcftools
	sudo apt install freebayes
	sudo apt install tabix
	sudo apt install bcftools

#sort the file post alignment
	samtools sort -o vir+host_unmapped+sorted.bam vir+host_unmapped.bam

#index the sorted 
	samtools index vir+host_unmapped+sorted.bam

#extract the tb reads 
	samtools faidx aligment_taxclass_phase2/bacteria_com25/bacteria_25.fna NC_000962.3 > H37Rv_ref.fna

#build bwa index for tb
	bwa index H37Rv_ref.fna
  
#map reads to tb ref
	bwa mem H37Rv_ref.fna bac_input_R1.fastq.gz bac_input_R2.fastq.gz > tb_mapped.sam

#sam to bam to sort and index
	samtools view -bS tb_mapped.sam | samtools sort -o tb_mapped_sorted.bam
	samtools index tb_mapped_sorted.bam
 
#index the fna file
	samtools faidx H37Rv_ref.fna

#variant calling 
	freebayes -f H37Rv_ref.fna tb_mapped_sorted.bam > tb_variants.vcf
	bgzip -c tb_variants.vcf > tb_variants.vcf.gz
	tabix -p vcf tb_variants.vcf.gz

#to filter the vcf variant file 

#filter by quality
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AD]\n' tb_variants.vcf.gz \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "CHROM","POS","REF","ALT","QUAL","REF_DEPTH","ALT_DEPTH","AF"}
    {
      split($6,a,",");
      if (a[1]+a[2]>0) af=a[2]/(a[1]+a[2]);
      else af=0;
      print $1,$2,$3,$4,$5,a[1],a[2],af;
    }' > tb_variants_with_AF.tsv


#Create classified TSV (HOMO / HETERO / LOWCONF)
	awk 'BEGIN{OFS="\t"; print "CHROM","POS","REF","ALT","QUAL","REF_DP","ALT_DP","DP","AF","CLASS"}
NR>1{
  dp=$6+$7; alt=$7; ref=$6; af=$8; qual=$5;
  cls="LOWCONF";
  # High-confidence homozygous-like (dominant strain)
  if( dp>=10 && af>=0.90 && alt>=5 ) cls="HOMO";
  # Moderate-confidence homozygous-like (relax alt requirement if AF very high)
  else if( dp>=8 && af>=0.85 && alt>=3 ) cls="HOMO";
  # Heterogeneous / minor strain candidates
  else if( dp>=6 && alt>=3 && af>=0.05 ) cls="HETERO";
  # Sensitive capture of plausible minor variants (lower alt but AF>=5%)
  else if( dp>=10 && alt>=2 && af>=0.05 ) cls="HETERO";
  # otherwise LOWCONF (likely noise)
  print $1,$2,$3,$4,qual,ref,alt,dp,af,cls
}' tb_variants_with_AF.tsv > tb_variants_classified.tsv


#Quick counts by class
	awk -F'\t' 'NR>1{c[$10]++} END{for(k in c) print k,c[k]}' tb_variants_classified.tsv


#annotate high confident variants 

#download and rename gff file
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz
gunzip GCF_000195955.2_ASM19595v2_genomic.gff.gz
mv GCF_000195955.2_ASM19595v2_genomic.gff H37Rv.gff

#annotate Homologous(high cofidence variants)
 awk -F'\t' '$10=="HOMO"{print $1"\t"$2}' tb_variants_classified.tsv > homo.txt
 bcftools view -R homo.txt tb_variants.vcf.gz -Oz -o tb_variants.HOMO.vcf.gz
 tabix -p vcf tb_variants.HOMO.vcf.gz


#SnpEff installation and setup
Download SnpEff:

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff

#Create SnpEff data directory for H37Rv
	mkdir -p data/Mycobacterium_tuberculosis_h37rv

#Download H37Rv reference genome and annotation
#Genome fasta:
	curl -L "https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/references/NC_000962.3.fasta" -o data/Mycobacterium_tuberculosis_h37rv/genome.fa

#Annotation (GFF):
	curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000

#Unzip GFF:
	gunzip data/Mycobacterium_tuberculosis_h37rv/genes.gff.gz

#Build SnpEff database for H37Rv

#Edit config to add genome entry:
	echo "Mycobacterium_tuberculosis_h37rv.genome : Mycobacterium_tuberculosis_h37rv" >> snpEff.config

#Build database:
	java -Xmx2g -jar snpEff.jar build -gff3 -v Mycobacterium_tuberculosis_h37rv

#Prepare VCF (chromosome rename fix)
H37Rv uses "Chromosome" as chromosome name, so modify VCF header accordingly:
gunzip tb_variants.HOMO.vcf.gz
sed -i 's/NC_000962.3/Chromosome/g' tb_variants.HOMO.vcf
bgzip tb_variants.HOMO.vcf
mv tb_variants.HOMO.vcf.gz tb_variants.HOMO.chrom.vcf.gz

#Index the fixed VCF:
	tabix tb_variants.HOMO.chrom.vcf.gz

#Run SnpEff annotation
	cd ~/snpEff
	java -Xmx2g -jar snpEff.jar -v Mycobacterium_tuberculosis_h37rv ~/tb_variants.HOMO.chrom.vcf.gz > tb_variants.HOMO.ann.

#Extract HIGH + MODERATE impact variants:
	java -jar SnpSift.jar filter "ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE'" tb_variants.HOMO.ann.vcf > tb_variants.HOMO.highModerate.vcf

#Make a final CSV report for HIGH/MODERATE:
java -jar SnpSift.jar extractFields tb_variants.HOMO.highModerate.vcf \
CHROM POS REF ALT AF \
GEN[unknown].DP GEN[unknown].AD \
"ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].HGVS_P" \
> report_highModerate.tsv


#Gene counts TSV
	awk '{print $10}' report_highModerate.tsv | sort | uniq -c | sort -nr | awk '{print $2 "\t" $1}' > top_gene_counts.tsv

	printf "GENE\tCOUNT\n" | cat - top_gene_counts.tsv > tmp && mv tmp top_gene_counts.tsv

#impact summary tsv 
	awk '{print $9}' report_highModerate.tsv | sort | uniq -c | \
awk '{print $2 "\t" $1}' > impact_summary.tsv

	
#Variant Effect Summary
	awk 'NR>1 {print $8}' report_highModerate.tsv \
  | sort \
  | uniq -c \
  | sort -nr \
  | awk '{print $2"\t"$1}' \
  > effect_counts.tsv

#high impact alone + count
	awk 'NR==1 || $9=="HIGH"' report_highModerate.tsv > report_HIGH.tsv
	awk 'NR>1 {print $8}' report_HIGH.tsv | sort | uniq -c | sort -nr > high_effect_counts.tsv

#r plots

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
