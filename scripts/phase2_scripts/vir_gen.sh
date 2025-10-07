#!/bin/bash

# Empty the output file first (so reruns don't duplicate entries)
> sub_vir10.fasta

for acc in \
NC_001405.1 \
NC_007605.1 \
NC_001806.2 \
NC_006273.2 \
NC_001348.1 \
NC_003977.2 \
NC_001526.4 \
NC_001802.1 \
NC_026431.1 \
NC_045512.2
do
  echo "Downloading $acc ..."
  wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" >> sub_vir10.fasta
done
