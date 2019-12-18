Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ ${BIN_DIR}/pancake seqdb test-1 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --block-size 0 --buffer-size 1024
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.0.seq test-1.seqdb.0.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.1.seq test-1.seqdb.1.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.2.seq test-1.seqdb.2.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.3.seq test-1.seqdb.3.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.4.seq test-1.seqdb.4.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb test-1.seqdb

Same as before, but test writing to a different folder. The file paths should be local.
  $ mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/test-1 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --block-size 0 --buffer-size 1024
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.0.seq out/test-1.seqdb.0.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.1.seq out/test-1.seqdb.1.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.2.seq out/test-1.seqdb.2.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.3.seq out/test-1.seqdb.3.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.4.seq out/test-1.seqdb.4.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb test-1.seqdb

Write the SeqDB like before, but store all sequences into the same single 2-bit compressed file.
  $ ${BIN_DIR}/pancake seqdb out ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --block-size 1024 --buffer-size 1024
  > cat ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.*.seq > expected.seq
  > cat out.seqdb
  > diff expected.seq out.seqdb.0.seq
  V	0.1.0
  C	1
  F	0	out.seqdb.0.seq	5	16560	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001

Using any size buffer should not impact the results, only perhaps the speed of writing. This will write out one read at a time.
  $ ${BIN_DIR}/pancake seqdb out ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --block-size 1024 --buffer-size 0
  > cat ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.*.seq > expected.seq
  > cat out.seqdb
  > diff expected.seq out.seqdb.0.seq
  V	0.1.0
  C	1
  F	0	out.seqdb.0.seq	5	16560	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001

Uncompressed SeqDB construction, each sequence in one block.
  $ ${BIN_DIR}/pancake seqdb test-3 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --compression 0 --block-size 0 --buffer-size 1024
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-3.seqdb.0.seq test-3.seqdb.0.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-3.seqdb.1.seq test-3.seqdb.1.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-3.seqdb.2.seq test-3.seqdb.2.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-3.seqdb.3.seq test-3.seqdb.3.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-3.seqdb.4.seq test-3.seqdb.4.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-3.seqdb test-3.seqdb

Uncompressed SeqDB construction, single block.
  $ ${BIN_DIR}/pancake seqdb out ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --compression 0 --block-size 1024 --buffer-size 1024
  > cat ${PROJECT_DIR}/test-data/seqdb-writer/test-3.seqdb.*.seq > expected.seq
  > cat out.seqdb
  > diff expected.seq out.seqdb.0.seq
  V	0.1.0
  C	0
  F	0	out.seqdb.0.seq	5	66233	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	5852	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	5852	11983	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	17835	24292	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	42127	5105	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	47232	19001	19001	1	0	19001
