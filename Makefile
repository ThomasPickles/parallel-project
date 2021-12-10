CC=mpicc
CFLAGS=-Wall -O3 -lm

TEST_DIR := ./tests
TEST_FILES := $(shell find $(TEST_DIR) -name 'data-*.txt')
size = 05 06
algos = prim-seq kruskal-seq prim-par kruskal-par
procs = 3 4

tests: mst
	mpirun -np 1 ./mst tests/data-06.txt prim-seq > tests/out-06-prim-seq.txt
	diff -w tests/out-06-prim-seq.txt tests/exp-06-prim.txt
	mpirun -np 1 ./mst tests/data-06.txt kruskal-seq > tests/out-06-krus-seq.txt
	diff -w tests/out-06-krus-seq.txt tests/exp-06-krus.txt
	mpirun -np 3 ./mst tests/data-06.txt prim-par > tests/out-06-prim-par.txt
	diff -w tests/out-06-prim-par.txt tests/exp-06-prim.txt

mst: mst-skeleton.c mst-solution.c
	$(CC) $(CFLAGS) -o $@ $< -lm

graph: create-graph.py
	python2.7 $< 6 14 5 graph.txt

run: mst
	mpirun -np 1 ./mst graph.txt kruskal-seq

clean:
	rm -rf *.o mst

.PHONY: clean

