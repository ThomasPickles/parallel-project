#!/usr/bin/env bash

# todo: expand these test cases
algos=('prim' 'kruskal')
flavour=('seq')
nodes=('06' '100')
# proc_number = 3 4

for alg in "${algos[@]}"; do
    for f in "${flavour[@]}"; do
    for n in "${nodes[@]}"; do
        data="./tests/data-$n.txt"
        expected="./tests/exp-$n-$alg.txt"
        actual="./tests/out-$n-$alg-$f.txt"
        touch $actual
        mpirun -np 1 ./mst $data $alg-$f > $actual
        diff -w $actual $expected
    done
    done
done
echo "Done"