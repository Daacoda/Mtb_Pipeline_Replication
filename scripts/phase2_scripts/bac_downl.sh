#!/bin/bash

# Output combined FASTA file
OUTFILE="bacteria_25.fna"

# Empty the output file first (so reruns don’t duplicate entries)
> $OUTFILE

# List of bacterial species
for species in \
"Streptococcus pneumoniae" \
"Mycobacterium tuberculosis" \
"Escherichia coli" \
"Salmonella enterica" \
"Bacillus subtilis" \
"Clostridium difficile" \
"Staphylococcus aureus" \
"Listeria monocytogenes" \
"Helicobacter pylori" \
"Vibrio cholerae" \
"Shigella dysenteriae" \
"Klebsiella pneumoniae" \
"Neisseria meningitidis" \
"Yersinia pestis" \
"Campylobacter jejuni" \
"Haemophilus influenzae" \
"Legionella pneumophila" \
"Treponema pallidum" \
"Chlamydia trachomatis" \
"Borrelia burgdorferi" \
"Enterococcus faecalis" \
"Corynebacterium diphtheriae" \
"Pseudomonas aeruginosa" \
"Acinetobacter baumannii" \
"Proteus mirabilis"
do
	echo "Downloading $species ..."
	fname=$(echo $species | tr  ' ' '_') #input underscore in filename

	#Download genome dataset(ref only)
	datasets download genome taxon "$species" --reference --filename ${fname}.zip

	# Unzip quietly into a temp directory
	 unzip -q -o ${fname}.zip -d ${fname}_dir

	# Append only the FASTA sequences (.fna) to the master file
        cat ${fname}_dir/ncbi_dataset/data/*/*.fna >> $OUTFILE

	# Delete zips and temp folder
	rm -rf ${fname}.zip ${fname}_dir
done

echo "All FASTAs combined into $OUTFILE"
