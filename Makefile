CC=mpicc
CFLAGS=-Wall -O3 -lm

all: mst

mst: mst-skeleton.c
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -rf *.o mst

.PHONY: clean

