The pair of reads is taken from the synthetic dataset ivan-200k-t1.

```
function pair1 {
    local out_prefix=reads1.weird_masking_issue
    echo "000000551" > ${out_prefix}.txt
    echo "000000571" >> ${out_prefix}.txt
    SEQ_DB_PREFIX="../RUN/02-build_db/reads"
    pancake seqfetch --write-ids --alias ${SEQ_DB_PREFIX}.seqdb ${out_prefix}.fasta ${out_prefix}.txt ${SEQ_DB_PREFIX}.seqdb
}
```
