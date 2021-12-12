CC=mpicc
CFLAGS=-Wall -O3 -lm

tests: mst run_cases.sh
	grep -ri 'TODO' .
	./run_cases.sh

charts: mst
	python3 make_charts.py

mst: mst-skeleton.c mst-solution.c
	$(CC) $(CFLAGS) -o $@ $< -lm

graph: create-graph.py
	python2.7 $< 6 14 5 graph.txt

run: mst
	mpirun -np 1 ./mst graph.txt kruskal-seq

clean:
	rm -rf *.o mst
	rm tests/out*

.PHONY: clean run_cases.sh
