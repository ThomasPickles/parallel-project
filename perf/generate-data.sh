while read line
do eval "python2.7 create-graph.py $line"
	# mpirun -np 2 ./mst ./tests/data-06.txt prim-par
done