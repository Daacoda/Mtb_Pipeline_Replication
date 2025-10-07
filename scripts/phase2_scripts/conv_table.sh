#!/bin/bash
# input: bacteria_accessions.txt
# output: conversion_table.txt

# Empty output file before writing
> conversion_table.txt

for acc in $(cat bacteria_accessions.txt); do
    echo "Processing $acc..."
    taxid=$(esearch -db nuccore -query "${acc}[Accession]" \
            | elink -target taxonomy \
            | efetch -format uid \
            | head -n 1)   # take only the first taxid
    echo -e "$acc\t$taxid" >> conversion_table.txt
done

