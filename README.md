# A Reproducible Low-Resource Workflow for NGS Variant Detection in *Mycobacterium tuberculosis*

**Author:** Olawuyi Samuel Babatunde  
**Location:** Ilorin, Nigeria  
**Repository:** https://github.com/Daacoda/Mtb_Pipeline_Replication  
**Date:** November 2025  

---

## Overview

This repository contains a lightweight, fully reproducible workflow for processing next-generation sequencing (NGS) data and detecting genomic variants in *Mycobacterium tuberculosis*.  
It is designed for **low-resource environments**, using standard, accessible tools and subsampling strategies to reduce computational load.  
It serves as a **training-focused pipeline** for students developing practical skills in microbial genomics and bioinformatics.

---

## Objective

To implement a reproducible and resource-efficient pipeline for:

- Quality control  
- Read filtering  
- Taxonomic profiling  
- Variant calling and annotation  

…using accessible tools suitable for low-power systems.

---

## Methods (Summary)

- **Dataset:** Public MTB sequencing dataset (SRR31065062)  
- **QC:** FastQC + Trimmomatic (Phred ≥30)  
- **Subsampling:** 400,000 total reads using Seqtk  
- **Filtering & Alignment:** Host removal + alignment to MTB reference genome using BWA  
- **Taxonomic Profiling:** Centrifuge classification  
- **Variant Calling:** FreeBayes + bcftools  
- **Annotation:** SnpEff classification by impact (HIGH, MODERATE, LOW)  
- **Reproducibility:** All phases scripted and environment documented

---

## Results (Summary)

- High-quality reads retained, GC content ~64%  
- >99% reads classified as *Mycobacterium tuberculosis*  
- **706 functional variants detected**, including:  
  - 72 HIGH-impact variants  
  - 634 MODERATE-impact variants  
- Includes frameshift, stop-gained, and missense mutations  
- All output tables (TSV, HTML) and scripts are included

---

## Significance

- Demonstrates a fully reproducible, low-resource MTB workflow  
- Provides practical experience in QC, alignment, filtering, annotation, and Bash automation  
- Useful for bioinformatics training and postgraduate research preparation  
- Works on modest hardware

---

## Limitations

- Subsampling reduces sensitivity to low-frequency variants  
- Small bacterial/viral reference panels limit broad microbial detection  
- Training-oriented, not intended for clinical or large-scale discovery

---

## Repository Structure

```text
Mtb_Pipeline_Replication/
├── README.md
├── environment.yml
├── project_report.txt
├── qc_phase1/
│   ├── multiqc_report/
│   ├── phase1_README.md
│   └── post_trim/
├── aligment_taxclass_phase2/
│   ├── alignment_summary.csv
│   ├── centrifuge_bac_class_report.tsv
│   ├── phase2_README.md
│   ├── Rplots/
│   │   ├── Rplot_alignment_pie.pdf
│   │   ├── Rplot_alignment_pie.png
│   │   ├── Rplot_bar_ch_map_eff.pdf
│   │   ├── Rplot_bar_ch_mapp_eff.png
│   │   ├── tax_class_bar_plot.png
│   │   └── tax_class_plot.pdf
│   └── bacteria_com25/
│       ├── bac_downl.sh
│       └── bacteria_25.fna
├── var_call_phase3/
│   (Final VCFs, SnpEff annotation files, TSV extracts, etc.)
├── plots/
│   ├── Rplot_bar_ch_mapp_eff.png
│   ├── effect_count.png
│   ├── fastqc_per_sequence_gc_content_plot.png
│   ├── fastqc_per_sequence_quality_scores_plot.png
│   ├── snpEff_summary.html
│   └── tax_abun.png
└── scripts/
    (All shell scripts used in the 3 phases)

---

## Reproducible Environment

conda env create -f environment.yml
conda activate tbqc
Tools Used
FastQC
Trimmomatic
Seqtk
BWA
Samtools / Bcftools
FreeBayes
Centrifuge
SnpEff / SnpSift
MultiQC

---

## License

This project is open for educational and research purposes.
