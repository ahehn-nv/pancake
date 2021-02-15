Tests dumping of the entire DB to FASTA.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb out.fasta
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Tests dumping of the entire DB to stdout.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb - > out.stdout.fasta
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta out.stdout.fasta | wc -l | awk '{ print $1 }'
  0

Test dumping of a particular block ID.
Note - the Awk command here concatenates the multiline FASTA into a single line because Samtools automatically wraps.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb out.fasta --block-id 2
  > samtools faidx ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" | awk '/^>/ {printf("%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test fetching the homopolymer compressed sequences.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-10a-compressed-in-2-small.seqdb - --use-hpc | sed 's/>//g'
  Seq1
  ATCGTCAGCGCG
  Seq2
  AGCGCAGCAGACATACGATATGCA
  Seq3
  C
  Seq4
  A
  Seq5
  ACTCGCGCAGTGACATCGCAGAGTATAGCAGCGTACAGCGTACGTG
