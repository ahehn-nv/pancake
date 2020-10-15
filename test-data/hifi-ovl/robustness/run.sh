#! /bin/bash

# Extracts the reads from the dataset.

seqdb=/pbi/flash/isovic/large-data-output-volatile/asm-regression-testing-2020_09/out-asm-ipa/02-dros/RUN/02-build_db/reads.seqdb

for f in reads.robustness.*
do
    echo $f
    pancake seqfetch ${f%.*}.fasta ${f} ${seqdb}
done
