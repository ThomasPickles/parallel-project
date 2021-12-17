#!/usr/bin/env bash
N=1000
algos='prim kruskal'
procs='1 2 4 8 16 32'
for algo in $algos; do
    outfile="./perf/results_${algo}.txt"
    touch $outfile
    header="writing data for N=$N"
    for p in $procs; do
        echo "running $algo on $p procs"
        o_procs="${p}"
        if [[ $p == '1' ]]; then
            f='seq'
        else
            f='par'
            # python2 smpi-generate-ring.py $p 100 100 100Gbps 1us
        fi
        echo "smpirun -hostfile smpi/ring-$p-hostfile.txt -platform smpi/ring-$p-platform.xml ./mst tests/data-$N.txt ${algo}-$f | tail -1"
        o_time=$(smpirun -hostfile smpi/ring-$p-hostfile.txt -platform smpi/ring-$p-platform.xml ./mst tests/data-$N.txt ${algo}-$f | tail -1)
        printf '%s\n' $o_procs $o_time | paste -sd ',' >> $outfile
    done
done
