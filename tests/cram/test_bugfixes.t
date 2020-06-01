Labeling of overlaps. There was a bug where one of the overlaps used to be marked as 3 prime instead of 5 prime.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/bugfixes/overlap_labelling/reads.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --skip-sym --write-rev --num-threads 1 reads reads 0 0 0
  pat_m64011_190605_003147/160434972/ccs pat_m64011_190605_003147/32900288/ccs -11690 99.76 0 0 11718 13995 1 0 11736 13966 5 * * * u
  pat_m64011_190605_003147/32900288/ccs pat_m64011_190605_003147/160434972/ccs -11690 99.76 0 0 11736 13966 1 0 11718 13995 5 * * * u

Labeling of overlaps, same as before but write IDs instead of headers.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/bugfixes/overlap_labelling/reads.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --skip-sym --write-rev --write-ids --num-threads 1 reads reads 0 0 0
  000000001 000000000 -11690 99.76 0 0 11718 13995 1 0 11736 13966 5 * * * u
  000000000 000000001 -11690 99.76 0 0 11736 13966 1 0 11718 13995 5 * * * u
