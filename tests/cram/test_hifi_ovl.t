Dovetail overlap 5prime test. Read "m64030_190330_071939/101844710/ccs" should have only 5prime overlaps with other reads.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile1-5prime.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/101844710/ccs"
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9078 99.72 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60884410/ccs -4958 99.60 0 0 4946 11811 0 4181 9139 9139 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61276980/ccs -3845 99.77 0 0 3845 11811 1 0 3838 10150 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61737670/ccs -6031 99.82 0 0 6031 11811 0 2564 8590 8590 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/62523930/ccs -7895 99.94 0 0 7895 11811 1 0 7894 11307 5

Dovetail overlap 3prime test. Read "m64030_190330_071939/101909220/ccs" should have only 3prime overlaps with other reads.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile2-3prime.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/101909220/ccs"
  m64030_190330_071939/101909220/ccs m64030_190330_071939/22611030/ccs -5435 99.10 0 4096 9530 9530 0 0 5435 9487 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/25626700/ccs -3286 99.91 0 6247 9530 9530 0 0 3286 10208 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/154862090/ccs -3003 99.53 0 6527 9530 9530 1 5043 8044 8044 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/157549180/ccs -1834 98.96 0 7696 9530 9530 1 8299 10128 10128 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/165937510/ccs -7423 99.89 0 2111 9530 9530 1 2475 9898 9898 3

Dovetail overlap all fwd oriented.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -6813 99.91 0 1793 8602 8602 0 0 6813 8913 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -6996 99.21 0 0 6996 8602 0 656 7651 7651 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -1328 99.85 0 7276 8602 8602 0 0 1328 9090 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -6836 99.84 0 0 6829 8602 0 1987 8823 8823 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1789 99.78 0 0 1787 8602 0 7139 8928 8928 5

Dovetail overlap all rev oriented.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile4-rev.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/102172020/ccs"
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -6182 99.84 0 4459 10635 10635 1 3541 9723 9723 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -9425 99.86 0 1217 10635 10635 1 946 10371 10371 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -9851 99.90 0 0 9847 10635 1 0 9851 10977 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -4547 99.67 0 0 4547 10635 1 0 4542 9121 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -3660 99.73 0 6981 10635 10635 1 7165 10825 10825 3

Perfect matching overlap.
One sequence is a real CCS read ("m64030_190330_071939/101844710/ccs"), one is a fake read composed of 5000bp copy/pasted
from the real read, and the third sequence is the reverse complement of the fake read.
All should have perfect 100.00% identity matches.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained

Test symmetry in the aligned identity when A->B and B->A.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile6-short-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --min-idt 96 reads reads 0 0 0
  m64030_190330_071939/39060280/ccs m64030_190330_071939/101975180/ccs -1962 96.53 0 5293 7255 9907 1 0 1961 10984 u
  m64030_190330_071939/101975180/ccs m64030_190330_071939/39060280/ccs -1962 96.53 0 0 1962 10984 1 5294 7255 9907 u

Contained overlap.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile7-contained.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0
  m64030_190330_071939/100008040/ccs m64030_190330_071939/117113640/ccs -9847 99.88 0 0 9847 9847 0 388 10231 11163 contained
  m64030_190330_071939/117113640/ccs m64030_190330_071939/100008040/ccs -9847 99.88 0 388 10231 11163 0 0 9847 9847 contains

Internal overlap.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile8-internal.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --min-idt 94.0 reads reads 0 0 0 | grep "^m64030_190330_071939/101451080/ccs"
  m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -3347 96.86 0 2140 5474 8022 1 1716 5063 12064 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -3370 95.37 0 2116 5474 8022 1 3400 6770 10056 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -3368 95.72 0 2116 5474 8022 0 5649 9017 11350 u
### Note: this differs from Raptor, Raptor is more sensitive here:
# m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -2633 96.7319 0 2256 7060 8022 1 126 4930 12064
# m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -2505 96.2606 0 2227 7415 8022 0 5753 10965 11350
# m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -2314 95.9715 0 2256 7444 8022 1 1428 6637 10056
# m64030_190330_071939/101451080/ccs m64030_190330_071939/165414370/ccs -1392 95.4594 0 4513 7390 8022 0 0 2797 10706
# m64030_190330_071939/101451080/ccs m64030_190330_071939/160039920/ccs -500 98.9011 0 6261 7357 8022 1 7586 8678 8678
### Also, to-do: larger '--aln-bw' actually can reduce sensitivity. Look into that.

Overlapping a single full pile.
(Awk is a hack to support cross-platform 'wc'. On some systems it outputs a prefix whitespace, and on some it does not.)
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile9-single-full-pile.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/111673510/ccs" | wc -l | awk '{ print $1 }'
  69

##################################################################
### Test removing symmetric arcs and writing reverse overlaps. ###
##################################################################
This input consists of a single dovetail overlap. By default, both arcs should be output.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9078 99.72 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9078 99.72 0 981 10059 10059 0 0 9077 11811 3

Skip symmetric overlaps. The pair where Aid > Bid should be removed. (The ID is given in the DB, not the Aname and Bname.)
This input consists of a single dovetail overlap. By default, both arcs should be output.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym reads reads 0 0 0
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9078 99.72 0 981 10059 10059 0 0 9077 11811 3

Writing symmetric overlaps should provide the other arc of the pair of reads, even though only one of them was
actually computed. Order may change though.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym --write-rev reads reads 0 0 0
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9078 99.72 0 981 10059 10059 0 0 9077 11811 3
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9078 99.72 0 0 9077 11811 0 981 10059 10059 5
##################################################################

Test writing IDs instead of headers. Other than that, the same as the Perfect match overlap test.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --write-ids reads reads 0 0 0
  000000000 000000001 -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  000000000 000000002 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  000000001 000000000 -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  000000001 000000002 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  000000002 000000000 -5000 100.00 0 0 5000 15000 1 0 5000 11811 5
  000000002 000000001 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained

