todo: simgrid processor array, find this in earlier TP
todo: python script to produce charts

todo: remember, this is a parallel algos course, so report should focus more on the asymptotics of parallelism, and not so much on the scaling as function of system size


## Speedup
### Data structures

- Use queue rather than array of status, then no need to check for each vertices.  Yields speedups when the vast majority of vertices are rejected as avoids repeated checking
- Union-find: follow pointers, discuss pointer jumping to give logarithmic behaviour.  Can I present a chart to illustrate this?
- Heaps to keep track of mins per vertex?  Is this too much data to store?

## Scaling

### As a function of graph size

### Kruskal vs Prim
we can't do better than O(N^2) since we're in an adjacency matrix representation.  Accoridingly, it's just the time it takes to read all the values.  If we had an adjacency list represnetation, we could (TODO: how?) improve our implementation as we wouldn't have the same limits.  Note on sparse vs dense graphs?

Should (todo: check) be same in theory.  My Prim throws away lots of useful calcs though

At what size do we start seeing asymptotic behaviour?  What constant are hidden in the small N behaviour?  Does one method overtake the other?

todo: do we reach a point where sparse vs dense is important?  I think that, in principle, since we're using adjacency matrix this should not change things todo: demo with case of m~O(n) [sparse] and m~O(n^2) [dense] - we're using a bad represenation if we're in a sparse regime though.

### Benefits of parallelism
what speedups do we see?  do we get to a regime where the inherent sequentialism of one method means that the T_seq/p*T_par (what is this called?)  starts to become a problem?  Can we study this for different values of p, L, \beta?

### Limits of parallelism
What is optimum p, q?  Is there any reason to favour are particular order of vertices for each processor?  I don't think so, as we're not passing to a previous, next.  When does communication cost become too high?  As a function of bandwidth, latency etc.

- Communication vs computation
TODO: chart showing when we reach a case where communication exceeds computuation

todo: vertex rather than node
Parallel, does p always divide n exactly?