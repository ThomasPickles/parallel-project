#!/usr/bin/env bash

# todo: expand these test cases
algos=('prim' 'kruskal')
nodes=('05' '06' '100')
proc_number=('1' '2' '3')

for alg in "${algos[@]}"; do
    for n in "${nodes[@]}"; do
    for p in "${proc_number[@]}"; do
        data="./tests/data-$n.txt"
        expected="./tests/exp-$n-$alg.txt"
        if [[ $p == '1' ]]; then
            f='seq'
        else
            f='par'
        fi
        actual="./tests/out-$n-$alg-$f-$p.txt"
        touch $actual
        mpirun -np $p ./mst $data $alg-$f > $actual
        diff -wq $actual $expected
    done
    done
done
echo "'diff -wy actual expected' to compare side-by-side"
echo "Done"