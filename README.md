todo: simgrid processor array, find this in earlier TP
todo: python script to produce charts



## Speedup
### Data structures
- Use queue rather than array of status, then no need to check for each vertices.  Yields speedups when the vast majority of vertices are rejected as avoids repeated checking
- Union-find: follow pointers, discuss pointer jumping to give logarithmic behaviour.  Can I present a chart to illustrate this?
- Heaps to keep track of mins per vertex?  Is this too much data to store?

## Scaling
### Kruskal vs Prim
Should (todo: check) be same in theory.  My Prim throws away lots of useful calcs though

At what size do we start seeing asymptotic behaviour?  What constant are hidden in the small N behaviour?  Does one method overtake the other?

### Benefits of parallelism
what speedups do we see?  do we get to a regime where the inherent sequentialism of one method means that the T_seq/p*T_par (what is this called?)  starts to become a problem?  Can we study this for different values of p, L, \beta?

### Limits of parallelism
What is optimum p, q?  Is there any reason to favour are particular order of vertices for each processor?  I don't think so, as we're not passing to a previous, next.  When does communication cost become too high?  As a function of bandwidth, latency etc.


todo: vertex rather than node
Parallel, does p always divide n exactly?