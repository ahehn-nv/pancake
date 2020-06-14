Dovetail overlap 5prime test. Read "m64030_190330_071939/101844710/ccs" should have only 5prime overlaps with other reads.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile1-5prime.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/101844710/ccs"
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/62523930/ccs -7889 99.94 0 0 7895 11811 1 0 7894 11307 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61737670/ccs -6015 99.82 0 0 6031 11811 0 2564 8590 8590 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60884410/ccs -4942 99.68 0 0 4959 11811 0 4181 9139 9139 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61276980/ccs -3829 99.77 0 0 3845 11811 1 0 3838 10150 5

Dovetail overlap 3prime test. Read "m64030_190330_071939/101909220/ccs" should have only 3prime overlaps with other reads.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile2-3prime.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/101909220/ccs"
  m64030_190330_071939/101909220/ccs m64030_190330_071939/165937510/ccs -7412 99.91 0 2111 9530 9530 1 2475 9898 9898 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/22611030/ccs -5388 99.15 0 4096 9530 9530 0 0 5435 9487 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/25626700/ccs -3280 99.91 0 6247 9530 9530 0 0 3286 10208 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/154862090/ccs -2988 99.57 0 6527 9530 9530 1 5043 8044 8044 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/157549180/ccs -1811 99.02 0 7696 9530 9530 1 8299 10128 10128 3

Dovetail overlap all fwd oriented.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -6944 99.27 0 0 6996 8602 0 656 7651 7651 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -6818 99.84 0 0 6829 8602 0 1987 8823 8823 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -6803 99.91 0 1793 8602 8602 0 0 6813 8913 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1783 99.78 0 0 1787 8602 0 7139 8928 8928 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -1324 99.85 0 7276 8602 8602 0 0 1328 9090 3

Dovetail overlap all rev oriented.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile4-rev.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/102172020/ccs"
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -9838 99.91 0 0 9847 10635 1 0 9851 10977 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -9406 99.87 0 1217 10635 10635 1 946 10371 10371 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -6167 99.85 0 4459 10635 10635 1 3541 9723 9723 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -4527 99.67 0 0 4547 10635 1 0 4542 9121 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -3645 99.75 0 6981 10635 10635 1 7165 10825 10825 3

Perfect matching overlap.
One sequence is a real CCS read ("m64030_190330_071939/101844710/ccs"), one is a fake read composed of 5000bp copy/pasted
from the real read, and the third sequence is the reverse complement of the fake read.
All should have perfect 100.00% identity matches.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

Test symmetry in the aligned identity when A->B and B->A.
TODO: This overlap is not symmetric because we could use trimming of the local alignments. It is context dependent how far
the alignment is extended in each direction.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile6-short-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --min-idt 0 reads reads 0 0 0
  m64030_190330_071939/39060280/ccs m64030_190330_071939/101975180/ccs -1914 94.88 0 5184 7255 9907 1 0 2020 10984 u
  m64030_190330_071939/101975180/ccs m64030_190330_071939/39060280/ccs -1929 94.38 0 0 2099 10984 1 5208 7255 9907 u

Contained overlap.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile7-contained.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/100008040/ccs m64030_190330_071939/117113640/ccs -9831 99.88 0 0 9847 9847 0 388 10231 11163 contained
  m64030_190330_071939/117113640/ccs m64030_190330_071939/100008040/ccs -9831 99.88 0 388 10231 11163 0 0 9847 9847 contains

Internal overlap.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile8-internal.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --min-idt 94.0 reads reads 0 0 0 | grep "^m64030_190330_071939/101451080/ccs"
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -5266 95.51 0 2168 7686 8022 0 5698 11212 11350 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -5232 95.48 0 2168 7648 8022 1 1238 6721 10056 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -4865 95.35 0 2043 7186 8022 1 0 5104 12064 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/160039920/ccs -1347 94.10 0 6263 7737 8022 1 7244 8678 8678 u
### Note2: these were the results with the SES-1986 algorithm. Internal overlaps don't have a clear boundary, so it's
### arbitrary how much we extend them beyond their "endpoint".
#  m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -3229 96.86 0 2140 5474 8022 1 1716 5063 12064 u
#  m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -3214 95.72 0 2116 5474 8022 0 5649 9017 11350 u
#  m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -3202 95.37 0 2116 5474 8022 1 3400 6770 10056 u
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
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Skip symmetric overlaps. The pair where Aid > Bid should be removed. (The ID is given in the DB, not the Aname and Bname.)
This input consists of a single dovetail overlap. By default, both arcs should be output.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym reads reads 0 0 0
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Writing symmetric overlaps should provide the other arc of the pair of reads, even though only one of them was
actually computed. Order may change though.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym --write-rev reads reads 0 0 0
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
##################################################################

Test writing IDs instead of headers. Other than that, the same as the Perfect match overlap test.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --write-ids reads reads 0 0 0 --out-fmt m4
  000000000 000000001 -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  000000000 000000002 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  000000001 000000002 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  000000001 000000000 -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  000000002 000000001 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  000000002 000000000 -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

##############
### Best N ###
##############
Bestn testing on a set of perfect matching overlaps. Value "0" should deactivate the filter.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn 0 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

Bestn testing on a set of perfect matching overlaps. Value < "0" should deactivate the filter.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn -1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

Bestn testing on a set of perfect matching overlaps. Positive value should be applied normally.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn 1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained

Bestn testing on a set of perfect matching overlaps. A large positive value, to test for potential edge cases.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn 10 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5
##############

###########################
### Traceback alignment ###
###########################
Traceback alignment.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -6968 99.27 0 0 6996 8602 0 656 7651 7651 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -6827 99.84 0 0 6829 8602 0 1987 8823 8823 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -6808 99.91 0 1793 8602 8602 0 0 6813 8913 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1786 99.78 0 0 1787 8602 0 7139 8928 8928 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -1326 99.85 0 7276 8602 8602 0 0 1328 9090 3

Traceback alignment.
Dovetail overlap all rev oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile4-rev.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --out-fmt m4 | grep "^m64030_190330_071939/102172020/ccs"
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -9844 99.91 0 0 9847 10635 1 0 9851 10977 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -9415 99.87 0 1217 10635 10635 1 946 10371 10371 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -6174 99.85 0 4459 10635 10635 1 3541 9723 9723 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -4537 99.67 0 0 4547 10635 1 0 4542 9121 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -3652 99.75 0 6981 10635 10635 1 7165 10825 10825 3

Traceback and HPC turned on at the same time.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --use-hpc --min-map-len 0 --min-anchor-span 100 --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 99.90 0 0 5145 6298 0 496 5640 5640 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 99.98 0 0 5026 6298 0 1470 6497 6497 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3

Traceback and HPC turned on at the same time.
Dovetail overlap all rev oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile4-rev.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --use-hpc --min-map-len 0 --min-anchor-span 100 --out-fmt m4 | grep "^m64030_190330_071939/102172020/ccs"
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -7275 99.99 0 0 7276 7845 1 0 7276 8160 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -6944 99.99 0 900 7845 7845 1 705 7650 7650 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -4530 99.98 0 3314 7845 7845 1 2617 7148 7148 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -3379 99.97 0 0 3379 7845 1 0 3380 6838 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -2680 99.96 0 5164 7845 7845 1 5338 8019 8019 3

Writing the CIGAR string produced by the traceback. HPC is turned on at the same time just so that the CIGAR string is shorter.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --use-hpc --min-map-len 0 --min-anchor-span 100 --write-cigar --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 99.90 0 0 5145 6298 0 496 5640 5640 5 221=1D649=1D296=1I169=1I167=1I3640=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 99.98 0 0 5026 6298 0 1470 6497 6497 5 5026=1D
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3 5008=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5 1288=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3 952=

Ignore indels when computing the alignment identity.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --no-indels --use-hpc --min-map-len 0 --min-anchor-span 100 --write-cigar --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 100.00 0 0 5145 6298 0 496 5640 5640 5 221=1D649=1D296=1I169=1I167=1I3640=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 100.00 0 0 5026 6298 0 1470 6497 6497 5 5026=1D
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3 5008=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5 1288=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3 952=

Ignore SNPs when computing the alignment identity.
This is a synthetic test case.
Note: The mismatches introduced into the fake read are all near to homopolymers, so the SES aligner can easily
mistake them for an indel combination with interspersed matching bases.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile11-3_mismatches.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --min-map-len 0 --min-anchor-span 100 --write-cigar --no-snps
  m64030_190330_071939/101844710/ccs fake_read_with_3_snps/1/ccs -4997 100.00 0 0 5000 11811 0 10000 15000 15000 5 3000=1X999=1X898=1X100=
  fake_read_with_3_snps/1/ccs m64030_190330_071939/101844710/ccs -4997 100.00 0 10000 15000 15000 0 0 5000 11811 3 3000=1X999=1X898=1X100=

Ignore indels when computing the alignment identity.
This is a synthetic test case.
Note: The mismatches introduced into the fake read are all near to homopolymers, so the SES aligner can easily
mistake them for an indel combination with interspersed matching bases.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile11-3_mismatches.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --min-map-len 0 --min-anchor-span 100 --write-cigar --no-indels --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read_with_3_snps/1/ccs -4997 99.94 0 0 5000 11811 0 10000 15000 15000 5 3000=1X999=1X898=1X100=
  fake_read_with_3_snps/1/ccs m64030_190330_071939/101844710/ccs -4997 99.94 0 10000 15000 15000 0 0 5000 11811 3 3000=1X999=1X898=1X100=

#############################
### Test variant strings. ###
#############################
Overlap two reads where one should be reverse complemented.
This is a synthetic test case, where there is only one mismatches and the sequences are short.
The same alignment should be reported in fwd and rev directions, and it should be the same.
The variant strings should always be in the fwd strand of the sequence, so when the reads are swapped, the
actual variant sequences should in this case be the same (just swapped).
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt ipa
  read1-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  read4-rev read1-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *

An expansion of the previous test, but with more than 2 sequences.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/hifi-ovl/reads.pile13-simple_errors_7_reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt ipa
  read1-fwd read2-fwd -180 100.0000 0 0 180 180 0 0 180 180 c c u 180= * * *
  read1-fwd read1-rev -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read1-fwd read3-fwd -179 99.4444 0 0 180 180 0 0 179 179 c c u 54=1I125= G * *
  read1-fwd read4-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= T C *
  read1-fwd read3-rev -179 99.4444 0 0 180 180 1 0 179 179 c c u 54=1I125= G * *
  read1-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  read2-fwd read1-fwd -180 100.0000 0 0 180 180 0 0 180 180 c c u 180= * * *
  read2-fwd read1-rev -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read2-fwd read3-fwd -179 99.4444 0 0 180 180 0 0 179 179 c c u 54=1I125= G * *
  read2-fwd read4-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= T C *
  read2-fwd read3-rev -179 99.4444 0 0 180 180 1 0 179 179 c c u 54=1I125= G * *
  read2-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  read3-fwd read1-fwd -179 99.4444 0 0 179 179 0 0 180 180 c c u 54=1D125= * G *
  read3-fwd read2-fwd -179 99.4444 0 0 179 179 0 0 180 180 c c u 54=1D125= * G *
  read3-fwd read1-rev -179 99.4444 0 0 179 179 1 0 180 180 c c u 54=1D125= * C *
  read3-fwd read3-rev -179 100.0000 0 0 179 179 1 0 179 179 c c u 179= * * *
  read3-fwd read4-fwd -178 98.8889 0 0 179 179 0 0 180 180 c c u 54=1D97=1X27= T GC *
  read3-fwd read4-rev -178 98.8889 0 0 179 179 1 0 180 180 c c u 54=1D97=1X27= T GC *
  read4-fwd read4-rev -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read4-fwd read1-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= C T *
  read4-fwd read2-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= C T *
  read4-fwd read1-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= C A *
  read4-fwd read3-fwd -178 98.8889 0 0 180 180 0 0 179 179 c c u 54=1I97=1X27= GC T *
  read4-fwd read3-rev -178 98.8889 0 0 180 180 1 0 179 179 c c u 54=1I97=1X27= GC A *
  read1-rev read1-fwd -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read1-rev read2-fwd -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read1-rev read3-fwd -179 99.4444 0 0 180 180 1 0 179 179 c c u 125=1I54= C * *
  read1-rev read4-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= A C *
  read1-rev read3-rev -179 99.4444 0 0 180 180 0 0 179 179 c c u 125=1I54= C * *
  read1-rev read4-rev -179 99.4444 0 0 180 180 0 0 180 180 c c u 27=1X152= A G *
  read3-rev read1-fwd -179 99.4444 0 0 179 179 1 0 180 180 c c u 125=1D54= * G *
  read3-rev read2-fwd -179 99.4444 0 0 179 179 1 0 180 180 c c u 125=1D54= * G *
  read3-rev read3-fwd -179 100.0000 0 0 179 179 1 0 179 179 c c u 179= * * *
  read3-rev read1-rev -179 99.4444 0 0 179 179 0 0 180 180 c c u 125=1D54= * C *
  read3-rev read4-fwd -178 98.8889 0 0 179 179 1 0 180 180 c c u 27=1X97=1D54= A GC *
  read3-rev read4-rev -178 98.8889 0 0 179 179 0 0 180 180 c c u 27=1X97=1D54= A GC *
  read4-rev read4-fwd -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read4-rev read1-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *
  read4-rev read2-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *
  read4-rev read1-rev -179 99.4444 0 0 180 180 0 0 180 180 c c u 27=1X152= G A *
  read4-rev read3-fwd -178 98.8889 0 0 180 180 1 0 179 179 c c u 27=1X97=1I54= GC T *
  read4-rev read3-rev -178 98.8889 0 0 180 180 0 0 179 179 c c u 27=1X97=1I54= GC A *
