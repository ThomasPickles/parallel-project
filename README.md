todo: python script to produce charts
todo: shell script to loop over prim and seq

# Qs


problem statement says we have to output equal weight edges in lexicographical order,
but it's not clear whether the edges should be listed globally from smallest to largest, and Prim could start by adding large weight edges.  For example, for Prim, since we select a random vertex to start, to there's no guarantee on the ordering.
do we have to output edges in particular order?  Using lexicographic order to break ties does seem to guarantee a particular output set though.

Parallel, does p always divide n exactly?