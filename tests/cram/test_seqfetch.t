Test fetching of sequences from two input files, one FASTA and one FASTQ.
Output in FASTA.
  $ rm -f out.* test.in.*
  > head -n 6 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta > test.in.1.fasta
  > tail -n 8 ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq > test.in.2.fastq
  > samtools faidx test.in.1.fasta
  > samtools fqidx test.in.2.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt test.in.1.fasta test.in.2.fastq
  > head -n 6 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta | tail -n 2 > expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from the SeqDB.
Output in FASTA.
  $ rm -f out.* test.in.*
  > ${BIN_DIR}/pancake seqdb test.in.1 ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt test.in.1.seqdb
  > head -n 6 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta | tail -n 2 > expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from a BAM file.
Output in FASTA.
  $ rm -f out.* test.in.*
  > echo "ref1/3/0_42176" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt ${PROJECT_DIR}/test-data/seqdb-writer/bam/subreads1.bam
  > samtools fasta ${PROJECT_DIR}/test-data/seqdb-writer/bam/subreads1.bam | head -n 6 | tail -n 2 > expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from two input files, one FASTA and one FASTQ.
Output in FASTA.
  $ rm -f out.* test.in.*
  > head -n 6 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta > test.in.1.fasta
  > tail -n 8 ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq > test.in.2.fastq
  > samtools faidx test.in.1.fasta
  > samtools fqidx test.in.2.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852" > test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt test.in.1.fasta test.in.2.fastq
  > grep ">" ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta | sort > expected.headers
  > grep ">" out.fasta | sort > results.headers
  > diff expected.headers results.headers | wc -l | awk '{ print $1 }'
  0

Output in FASTQ.
  $ rm -f out.* test.in.*
  > cp ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fastq out.fastq test.in.list.txt test.in.1.fastq
  > head -n 12 ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq | tail -n 4 > expected.fastq
  > diff expected.fastq out.fastq | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from a FOFN.
  $ rm -f out.* test.in.*
  > head -n 6 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta > test.in.1.fasta
  > tail -n 8 ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq > test.in.2.fastq
  > samtools faidx test.in.1.fasta
  > samtools fqidx test.in.2.fastq
  > echo "test.in.1.fasta" > input.fofn
  > echo "test.in.2.fastq" >> input.fofn
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852" > test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt input.fofn
  > grep ">" ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta | sort > expected.headers
  > grep ">" out.fasta | sort > results.headers
  > diff expected.headers results.headers | wc -l | awk '{ print $1 }'
  0

Test fetching using alias. Output in FASTA.
  $ rm -f out.* test.in.*
  > cp ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > ${BIN_DIR}/pancake seqdb test.in.1 test.in.1.fastq
  > echo "3" > test.in.list.txt
  > echo "4" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta --alias test.in.1.seqdb out.fasta test.in.list.txt test.in.1.fastq
  > echo ">m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105" > expected.headers
  > echo ">m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> expected.headers
  > cat expected.headers | sort > expected.headers.sorted
  > grep ">" out.fasta | sort > results.headers
  > diff expected.headers.sorted results.headers | wc -l | awk '{ print $1 }'
  0

Nonexistent sequence, but do not fail.
  $ rm -f out.* test.in.*
  > cp ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > echo "some_read_which_does_not_exist" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fastq out.fastq test.in.list.txt test.in.1.fastq
  > wc -l out.fastq | awk '{ print $1 }'
  0

Nonexistent sequence, FAIL when it cannot be found.
  $ rm -f out.* test.in.*
  > cp ${PROJECT_DIR}/test-data/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > echo "some_read_which_does_not_exist" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --fail --out-fmt fastq out.fastq test.in.list.txt test.in.1.fastq 2>&1 | sed 's/.*pancake //g'
  seqfetch ERROR: Not all queries were found in the provided input files! Found sequences: 0 / 1 .
