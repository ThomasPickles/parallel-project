CC=mpicc
CFLAGS=-Wall -O3 -lm

tests: mst run_cases.sh
	grep -ri 'TODO' .
	./run_cases.sh

charts: mst gen_perf.py generate-data.sh
	python3 gen_perf.py > perf-cases.txt
	cat perf-cases.txt | generate-data.sh
	python3 make_charts.py

run: mst
	mpirun -np 3 ./mst ./tests/data-06.txt kruskal-par

mst: mst-skeleton.c mst-solution.c
	$(CC) $(CFLAGS) -o $@ $< -lm

graph: create-graph.py
	python2.7 $< 6 14 5 graph.txt


clean:
	rm -rf *.o mst
	rm tests/out*

.PHONY: clean run_cases.sh
