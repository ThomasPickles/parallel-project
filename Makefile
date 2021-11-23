CC=mpicc
CFLAGS=-Wall -O3 -lm

all: mst

mst: mst-skeleton.c mst-solution.c
	$(CC) $(CFLAGS) -o $@ $< -lm

graph: create-graph.py
	python2.7 $< 5 10 5 graph.txt

run: mst
	mpirun -np 1 ./mst graph.txt prim-seq

clean:
	rm -rf *.o mst

.PHONY: clean

