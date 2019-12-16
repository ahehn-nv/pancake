Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ ${BIN_DIR}/pancake seqdb test-1 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --block-size 0 --buffer-size 1024
  > cat test-1.seqdb
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.0.seq test-1.seqdb.0.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.1.seq test-1.seqdb.1.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.2.seq test-1.seqdb.2.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.3.seq test-1.seqdb.3.seq
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.4.seq test-1.seqdb.4.seq
  V	0.1.0
  C	1
  F	0	test-1.seqdb.0.seq	1	1463	5852
  F	1	test-1.seqdb.1.seq	1	2996	11983
  F	2	test-1.seqdb.2.seq	1	6073	24292
  F	3	test-1.seqdb.3.seq	1	1277	5105
  F	4	test-1.seqdb.4.seq	1	4751	19001
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	1	0	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	2	0	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	3	0	1277	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	4	0	4751	19001	1	0	19001

Write the SeqDB like before, but store all sequences into the same single 2-bit compressed file.
  $ ${BIN_DIR}/pancake seqdb test-1 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --block-size 1024 --buffer-size 1024
  > cat ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.*.seq > expected.seq
  > cat test-1.seqdb
  > diff expected.seq test-1.seqdb.0.seq
  V	0.1.0
  C	1
  F	0	test-1.seqdb.0.seq	5	16560	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001

Using any size buffer should not impact the results, only perhaps the speed of writing. This will write out one read at a time.
  $ ${BIN_DIR}/pancake seqdb test-1 ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta --block-size 1024 --buffer-size 0
  > cat ${PROJECT_DIR}/test-data/seqdb-writer/test-1.seqdb.*.seq > expected.seq
  > cat test-1.seqdb
  > diff expected.seq test-1.seqdb.0.seq
  V	0.1.0
  C	1
  F	0	test-1.seqdb.0.seq	5	16560	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
