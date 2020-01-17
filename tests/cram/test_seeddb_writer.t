Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ ${BIN_DIR}/pancake seeddb -k 30 -w 80 ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb test-1a
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1a.seeddb test-1a.seeddb
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1a.seeddb.0.seeds test-1a.seeddb.0.seeds
  > ls -1 test-1a.seeddb*
  test-1a.seeddb
  test-1a.seeddb.0.seeds

Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ ${BIN_DIR}/pancake seeddb -k 30 -w 80 --split-blocks ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb test-1b
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1b.seeddb test-1b.seeddb
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1b.seeddb.0.seeds test-1b.seeddb.0.seeds
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1b.seeddb.1.seeds test-1b.seeddb.1.seeds
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1b.seeddb.2.seeds test-1b.seeddb.2.seeds
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1b.seeddb.3.seeds test-1b.seeddb.3.seeds
  > diff ${PROJECT_DIR}/test-data/seeddb-writer/test-1b.seeddb.4.seeds test-1b.seeddb.4.seeds
  > ls -1 test-1b.seeddb*
  test-1b.seeddb
  test-1b.seeddb.0.seeds
  test-1b.seeddb.1.seeds
  test-1b.seeddb.2.seeds
  test-1b.seeddb.3.seeds
  test-1b.seeddb.4.seeds
