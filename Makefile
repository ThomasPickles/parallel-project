CC=mpicc
CFLAGS=-Wall -O3 -lm

#   make TEST_FILES='TP03/**/bad*.c' tests
ifdef TEST_FILES
export TEST_FILES
endif

all: tests mst

tests: mst test_06

test_06: test_%:
	mpirun -np 1 ./mst tests/test_06.txt kruskal-seq
	# mpirun -np 3 ./mst tests/test_$*.txt prim-par

single: mst
	mpirun -np 1 ./mst tests/test_06.txt kruskal-seq > tests/krus_seq_out_06.txt
	diff -y tests/krus_seq_out_06.txt tests/exp_06.txt

all_tests: mst
	mpirun -np 1 ./mst tests/test_06.txt kruskal-seq
	# mpirun -np 1 ./mst tests/test_$*.txt prim-seq > tests/prim_seq_out_$*.txt
	# mpirun -np 1 ./mst tests/test_$*.txt kruskal-seq > tests/krus_seq_out_$*.txt
	# mpirun -np 3 ./mst tests/test_$*.txt prim-par > tests/prim_par_out_$*.txt
	# diff -y tests/prim_seq_out_$*.txt tests/exp_$*.txt
	# # diff command returns 1 if error found, so fails at first wrong one
	# diff -y tests/krus_seq_out_$*.txt tests/exp_$*.txt
	# diff -y tests/prim_par_out_$*.txt tests/exp_$*.txt

mst: mst-skeleton.c mst-solution.c
	$(CC) $(CFLAGS) -o $@ $< -lm

graph: create-graph.py
	python2.7 $< 6 14 5 graph.txt

run: mst
	mpirun -np 1 ./mst graph.txt kruskal-seq

clean:
	rm -rf *.o mst

.PHONY: clean test_06

