CC=mpicc
CFLAGS=-Wall -O3 -lm

tests: mst run_cases.sh
	$(CC) $(CFLAGS) -DTEST -o mst mst-skeleton.c -lm
	grep -ri 'TODO' .
	./run_cases.sh
	touch mst-skeleton.c

performance: mst run_perf.sh
	smpicc $(CFLAGS) -o mst mst-skeleton.c
	./run_perf.sh


charts: mst perf/gen_perf.py perf/generate-data.sh
	python3 perf/gen_perf.py > perf/perf-cases.txt
	cat perf/perf-cases.txt | perf/generate-data.sh
	python3 charts/make_charts.py

run: mst
	mpirun -np 2 ./mst ./tests/data-06.txt prim-par

mst: mst-skeleton.c mst-solution.c
	$(CC) $(CFLAGS) -o $@ $< -lm

graph: create-graph.py
	python2.7 $< 6 14 5 graph.txt


clean:
	rm -rf *.o mst
	rm tests/out*

.PHONY: clean run_cases.sh
