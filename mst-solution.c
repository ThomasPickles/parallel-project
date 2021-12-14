/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algo_name the name of the algorithm to be executed
 */

int VERBOSE = 0;
void print_edges(int *edges, int count)
{
  printf("edge count %d :[", count);
  for (int j = 0; j < count; j++)
  {
    printf("(%d, %d-%d) ", edges[3 * j], edges[3 * j + 1], edges[3 * j + 2]);
  }
  printf("]\n");
}
void output_edges(int *edges, int count)
{
  for (int j = 0; j < count; j++)
  {
    printf("%d %d\n", edges[3 * j + 1], edges[3 * j + 2]);
  }
}
int has_edge(int weight)
{
  return weight > 0;
}
int previous_candidate(int min_weight)
{
  return min_weight != 10;
}

int find_set(int *parent_array, int elem)
{
  int parent = parent_array[elem];
  if (parent == elem)
  {
    return elem;
  }
  else
  {
    // Optim: CAN MAKE MORE EFFICIENT BY UPDATING POINTER
    // EACH TIME WE GET LEADER.  MAKES TREE MUCH SHALLOWER
    parent = find_set(parent_array, parent);
    parent_array[elem] = parent;
    return parent;
  }
}
void union_sets(int *parents_array, int t1, int t2)
{
  // usually have t1 < t2, so join sets to lower numbers
  int l1 = find_set(parents_array, t1);
  int l2 = find_set(parents_array, t2);
  parents_array[l2] = l1;
}
void print_vertices(int *arr, int count, int log_level)
{
  if (log_level > 0)
  {
    printf("stat:[");
    for (int i = 0; i < count - 3; i = i + 3)
    {
      printf("%d,", arr[i]);
    }
    printf("%d]\n", arr[count - 3]);
    printf("weig:[");
    for (int i = 1; i < count - 3; i = i + 3)
    {
      printf("%d,", arr[i]);
    }
    printf("%d]\n", arr[count - 2]);
    printf("edge:[");
    for (int i = 2; i < count - 1; i = i + 3)
    {
      printf("%d,", arr[i]);
    }
    printf("%d]\n", arr[count - 1]);
  }
}
void print_edge(int e1, int e2)
{
  int min, max;
  min = e1 < e2 ? e1 : e2;
  max = e1 < e2 ? e2 : e1;
  printf("%d %d\n", min, max);
}

int compare_elem_at_address(const void *v1, const void *v2)
{
  const int i1 = **(const int **)v1;
  const int i2 = **(const int **)v2;
  return i1 < i2 ? -1 : (i1 > i2);
}
void sort_edges(int *sorted_edges, int *edges, int count)
{

  int **loc = calloc(count, sizeof(int *)); // initialise array of pointers
  for (int i = 0; i < count; i++)
  {
    loc[i] = &edges[3 * i];
  }

  // sort edge list
  // loc -> iteration order
  if (VERBOSE > 0)
  {
    print_vertices(edges, 3 * count, 1);
    puts("Before sorting:");
    for (int i = 0; i < count; i++)
    {
      printf("value: %-6d adress: %p\n", *loc[i], loc[i]);
    }
  }

  // edge order
  qsort(loc, count, sizeof loc, compare_elem_at_address);

  if (VERBOSE > 0)
  {
    puts("After sorting:");
    for (int i = 0; i < count; i++)
    {
      printf("value: %-6d posn: %d\n", *loc[i], (int)(loc[i] - edges));
    }
  }

  int idx;
  for (int i = 0; i < count; i++)
  {
    idx = (int)(loc[i] - edges);
    sorted_edges[3 * i] = edges[idx];
    sorted_edges[3 * i + 1] = edges[idx + 1];
    sorted_edges[3 * i + 2] = edges[idx + 2];
  }
  free(loc);
}
void sort_array(int **pointer_array, int count)
{
}
int compare_edges(int *e1, int *e2)
{
  int e1_w = *e1, e1_1 = *(e1 + 1), e1_2 = *(e1 + 2);
  int e2_w = *e2, e2_1 = *(e2 + 1), e2_2 = *(e2 + 2);
  e1_1 = e1_1 < e1_2 ? e1_1 : e1_2;
  e1_2 = e1_1 < e1_2 ? e1_2 : e1_1;
  e2_1 = e2_1 < e2_2 ? e2_1 : e2_2;
  e2_2 = e2_1 < e2_2 ? e2_2 : e2_1;
  if (e1_w < e2_w)
    return 1; // first better
  if (e1_w > e2_w)
    return -1; // second better
  if (e1_1 < e2_1)
    return 1;
  if (e1_1 > e2_1)
    return -1;
  if (e1_2 < e2_2)
    return 1;
  if (e1_2 > e2_2)
    return -1;
  return 0; // match
}

int is_better_edge(int candidate_w, int best_w, int candidate_1, int candidate_2, int best_1, int best_2)
{
  if (candidate_w == 0)
    return 0; // not an edge
  if (best_w == 0)
    return 1; // no competition
  if (candidate_w > best_w)
    return 0;
  if (candidate_w < best_w)
    return 1;
  // same weight, so break ties with lexicographical order
  int best_low, best_high, candidate_low, candidate_high;
  best_low = best_1 < best_2 ? best_1 : best_2;
  best_high = best_1 < best_2 ? best_2 : best_1;
  candidate_low = candidate_1 < candidate_2 ? candidate_1 : candidate_2;
  candidate_high = candidate_1 < candidate_2 ? candidate_2 : candidate_1;
  if (candidate_low < best_low)
    return 1;
  if (candidate_low > best_low)
    return 0;
  if (candidate_high < best_high)
    return 1;
  return 0;
}

void create_edge_list(int *ret, int *edges, int *adj, int N, int q, int parallel, int proc_rank)
{
  // adjancency matrix is symmetric, so only consider upper
  // triangular part
  int *ptr = edges; //, *ptr1 = v1, *ptr2 = v2;
  int i_local, i_global, j, w;
  int count = 0;
  int start;
  for (i_local = 0; i_local < q; i_local++)
  {
    if (parallel)
    {
      i_global = i_local + q * proc_rank;
    }
    else
    {
      i_global = i_local;
    }
    // start at i+1 to avoid self loops
    start = parallel ? 0 : i_local + 1;
    for (j = start; j < N; j++)
    {
      // not a great implementation of this conditional
      // right inside a double for-loop :(
      if (i_global == j)
        continue;
      w = adj[i_local * N + j];
      if (w > 0)
      {
        // update and increment
        *ptr++ = w;                           // weights
        *ptr++ = i_global < j ? i_global : j; // v1
        *ptr++ = i_global < j ? j : i_global; // v2
        count++;
      }
    }
  }
  *ret = count;
}

int merge(int *out_set, int *set1, int *set2, int M1, int M2)
{
  // both sets should be sorted in weight order, so merge in turn

  int *p1 = set1, *p2 = set2;
  int *end1 = p1 + 3 * M1, *end2 = p2 + 3 * M2;
  int i = 0;
  int cmp = 0;

  while (p1 < end1 || p2 < end2)
  {
    // elems in both lists
    if (p1 < end1 && p2 < end2)
    {

      cmp = compare_edges(p1, p2);
      if (cmp == 1)
      {
        out_set[i++] = *p1++;
        out_set[i++] = *p1++;
        out_set[i++] = *p1++;
      }
      else if (cmp == -1)
      {
        out_set[i++] = *p2++;
        out_set[i++] = *p2++;
        out_set[i++] = *p2++;
      }
      else
      {
        out_set[i++] = *p1++;
        out_set[i++] = *p1++;
        out_set[i++] = *p1++;
        p2 = p2 + 3;
      }
    }
    else if (p1 < end1)
    {
      out_set[i++] = *p1++;
      out_set[i++] = *p1++;
      out_set[i++] = *p1++;
    }
    else
    {
      out_set[i++] = *p2++;
      out_set[i++] = *p2++;
      out_set[i++] = *p2++;
    }
  }
  return i / 3;
}

int kruskal(int *out_buf, int edge_count, int *sorted_edges, int N, int *vertex_list)
{
  int *leaders = calloc(N, sizeof(int));
  int idx, w;
  int mst = 0;
  int count = 0;
  // todo: vertices to include

  for (int i = 0; i < N; i++)
  {
    leaders[i] = i;
  }
  // todo: remove added variable if it doesn't do anything
  int *added = calloc(N, sizeof(int));
  int t1, t2, l1, l2;
  for (idx = 0; idx < edge_count; idx++)
  {
    w = sorted_edges[3 * idx];
    t1 = sorted_edges[3 * idx + 1]; // t1 < t2
    t2 = sorted_edges[3 * idx + 2];

    l1 = find_set(leaders, t1);
    l2 = find_set(leaders, t2);

    if (l1 != l2)
    {
      // different set, add edge to combine sets

      // combine sets
      union_sets(leaders, t1, t2);
      // add vertices
      added[t1] = 1;
      added[t2] = 1;
      out_buf[3 * count] = w;
      out_buf[3 * count + 1] = t1;
      out_buf[3 * count + 2] = t2;
      mst += w;
      count++;
    }
    else
    {
      // Edge from %d to %d in same set with leader %d.  Ignore
    }
  }

  free(added);
  free(leaders);
  return count;
}

void update_mins_for_vertices(int *vertices, int *adj, int last_added, int N)
{

  int w, min_weight, status;
  // Update mins for each vertex
  for (int v = 0; v < N; v++)
  {
    min_weight = 0;
    status = vertices[3 * v];
    if (status != 2)
    {
      // no point adding a frond edge
      w = adj[last_added + v * N];
      min_weight = vertices[3 * v + 1];
      if (is_better_edge(w, min_weight, v, last_added, v, vertices[3 * v + 2]))
      {
        if (v != last_added)
        { // no self-loops
          vertices[3 * v] = 1;
          vertices[3 * v + 1] = w;
          vertices[3 * v + 2] = last_added;
        }
      }
    }
  }
}

void get_abs_min(int *last_added, int *vertices, int N)
{
  // Get absolute min
  // find the absolute lowest weight among candidate edges
  // TODO: possible optimisation:
  // might be better to keep an "adjacency set to iterate over", rahter than reading
  // the vertices array again.  Although this is only really a O(N) step, so probably
  // not a big deal

  int u, v, w, w_min = 0, u_min = 0, v_min = 0, status;
  for (v = 0; v < N; v++)
  {
    status = vertices[3 * v];
    if (status != 2)
    {
      w = vertices[3 * v + 1];
      u = vertices[3 * v + 2];
      if (is_better_edge(w, w_min, u, v, u_min, v_min))
      {
        // update min
        v_min = v;
        u_min = u;
        w_min = w;
      }
    }
  }

  print_edge(u_min, v_min);

  // remove edge as candidate
  vertices[3 * v_min] = 2;
  vertices[3 * v_min + 1] = 0;
  *last_added = v_min;
}

int get_global_id(int local_id, int vertices_per_proc, int proc_rank)
{
  return local_id + vertices_per_proc * proc_rank;
}
void prim(int *vertices, int N, int *adj, int start_vertex)
{
  int last_added = start_vertex;
  vertices[3 * last_added] = 2;

  for (int vertex_count = 1; vertex_count < N; vertex_count++)
  {
    update_mins_for_vertices(vertices, adj, last_added, N);
    get_abs_min(&last_added, vertices, N);
  }
}
void compute_mst(
    int N,
    int M,
    int *adj,
    char *algo_name)
{
  int proc_rank = 0, nb_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);

  // TODO: HANDLE CASE WHERE N MOD P != 0

  // todo: interesting idea is assigning "levels" for the edges,
  // BRANCH, REJECTED, BASIC, so we can mark the rejected eges
  // and the unseen ones

  // TODO: CONVERT INTs TO UNSIGNED INT
  if (strcmp(algo_name, "prim-seq") == 0)
  { // Sequential Prim's algorithm
    if (proc_rank == 0 && nb_procs != 1)
    {
      printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", nb_procs);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // BEGIN IMPLEMENTATION HERE

    // TODO: SPEEDUP USING HEAP?
    // adjacency matrix size N**2
    // accordingly might not be feasible to convert to binary heap as this
    // requires adjacency list representation
    int *vertices = calloc(3 * N, sizeof(int));

    // (added, min_weight_edge, to_vertex)
    for (int i = 0; i < N; i++)
    {
      // 0 - unseen; 1 - adjacent; 2 - added
      vertices[3 * i] = 0;      // not added
      vertices[3 * i + 1] = 0;  // min weight
      vertices[3 * i + 2] = -1; // key
    }
    prim(vertices, N, adj, 0);
  }
  else if (strcmp(algo_name, "kruskal-seq") == 0)
  { // Sequential Kruskal's algorithm
    if (proc_rank == 0 && nb_procs != 1)
    {
      printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", nb_procs);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // BEGIN IMPLEMENTATION HERE
    // TODO: DO I WANT TO USE CONST FOR SOME OF THESE ARRAY POINTERS?

    // To check for set membership for vertices, just keep an array of length N
    //  initialised to 0, and set value to 1 if vertex is added to set.

    int *edges = calloc(3 * M, sizeof(int));
    int *sorted_edges = calloc(3 * M, sizeof(int));
    int count = 0; // NOT THE SAME AS M, SINCE WE IGNORE SELF-LOOPS

    // todo: check all pointers are allocated properly
    // todo: free all memory, particularly anything allocated within a loop
    if (edges == NULL)
    {
      fprintf(stderr, "malloc failed\n");
    }

    int *vertices_to_include = calloc(N, sizeof(int));

    // do for all vertices
    for (int i = 0; i < N; i++)
    {
      vertices_to_include[i] = 1;
    }

    create_edge_list(&count, edges, adj, N, N, 0, 0);
    sort_edges(sorted_edges, edges, count);
    count = kruskal(edges, count, sorted_edges, N, vertices_to_include);
    output_edges(edges, count);

    free(edges);

    //////
  }
  else if (strcmp(algo_name, "prim-par") == 0)
  { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE
    if (N % nb_procs != 0)
    {
      printf("Procs doesn't exactly divide vertices, exiting...");
      return;
    }
    int root = 0;
    int q = N / nb_procs;
    int start_vertex = 0; // less than q
    // shared among all processes
    int *added = calloc(N, sizeof(int));
    int *min_weight = calloc(q, sizeof(int));

    int *vertices = calloc(3 * nb_procs, sizeof(int));
    int proc_min[3] = {0};
    int i, u, w, next, v; //

    if (proc_rank == root)
    {
      next = start_vertex;
    }

    for (int vertex_count = 0; vertex_count < N - 1; vertex_count++)
    {
      MPI_Bcast(&next, 1, MPI_INT, root, MPI_COMM_WORLD);

      // add to set
      added[next] = 1;

      // get minimum for process, do on all processes
      // I'm concerned that we can't just keep the minimum, since when this
      // gets sucked up into the tree, how do we know what is the second
      // smallest? For now we can just do the double for-loop, and perhaps
      // look at some optimatisions later (make it run, make it right, make it fast)
      proc_min[0] = 0;
      proc_min[1] = 0;
      proc_min[2] = 0;

      for (i = 0; i < q; i++)
      {
        // yuck, this is the slow step, we're throwing away all the previous
        // calculations.  can we keep a heap somehow?
        min_weight[i] = 0;

        v = get_global_id(i, q, proc_rank);

        // todo: use prim seq subroutine to help with this - lots of duplication here
        if (added[v] == 1)
        {
          min_weight[i] = 0;
          continue;
        }

        for (u = 0; u < N; u++)
        {

          if (added[u] == 0)
          {
            continue;
          }

          if (v == u)
            // don't include self loops,
            continue;

          w = adj[u + i * N];
          if (w == 0)
            continue; // only want edges

          if (w < min_weight[i] || min_weight[i] == 0) // and edge is lower weight
          {
            // update its key
            min_weight[i] = w;
            // TODO: feels costly to check whether we've set the proc_min
            // to zero every time, since we just want to do it once
            // for the first edge we find
            // then, if it's the lowest weight yet seen by process
            if (w < proc_min[0] || proc_min[0] == 0)
            {
              // printf("lower weight edge R%dP%d %d:%d-%d\n", vertex_count, proc_rank, edge, v, u);
              proc_min[0] = w;
              // Put into lexicographic order
              proc_min[1] = v < u ? v : u;
              proc_min[2] = v < u ? u : v;
            }
          }
        }
      }
      // send mins back to root
      // TODO: can just do a reduce here, as we can discard any
      // values that are less than the minimum
      MPI_Gather(&proc_min, 3, MPI_INT, vertices, 3, MPI_INT, root, MPI_COMM_WORLD);

      if (proc_rank == root)
      {
        // todo: take this decl out of loop
        int w_min = 0, u_min = 0, v_min = 0;
        int w, u, v;
        for (int i = 0; i < nb_procs; i++)
        {
          w = vertices[3 * i];
          u = vertices[3 * i + 1];
          v = vertices[3 * i + 2];
          if (is_better_edge(w, w_min, u, v, u_min, v_min))
          {
            w_min = w;
            u_min = u;
            v_min = v;
          }
        }
        if (added[u_min])
        {
          next = v_min;
        }
        else
        {
          next = u_min;
        }
        print_edge(u_min, v_min);
      }
    }
  }
  else if (strcmp(algo_name, "kruskal-par") == 0)
  { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

    if (N % nb_procs != 0)
    {
      printf("Procs doesn't exactly divide vertices, exiting...");
      return;
    }
    int root = 0;
    int q = N / nb_procs;
    int MAX_EDGES = N * q; // in case fully connected
    int *include = calloc(N, sizeof(int));
    for (int i = proc_rank * q; i < (proc_rank + 1) * q; i++)
    {
      include[i] = 1;
    }

    // don't actually know how many edges we'll have on each processor,
    // and it many be more than M / q, so allocate
    // enough space for a fully-connected graph in the worst case
    int *edges = calloc(3 * MAX_EDGES, sizeof(int));
    int *sorted_edges = calloc(3 * MAX_EDGES, sizeof(int));
    int *counts = calloc(nb_procs, sizeof(int));
    int count = 0;

    // TODO: this is a huge waste of memory, since we only need it
    // on root proc
    int *edge_buf = calloc(3 * nb_procs * MAX_EDGES, sizeof(int));

    create_edge_list(&count, edges, adj, N, q, 1, proc_rank);
    sort_edges(sorted_edges, edges, count);

    MPI_Gather(&count, 1, MPI_INT, counts, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Gather(sorted_edges, 3 * MAX_EDGES, MPI_INT, edge_buf, 3 * MAX_EDGES, MPI_INT, root, MPI_COMM_WORLD);

    if (proc_rank == root)
    {

      int *out = calloc(3 * MAX_EDGES, sizeof(int));      // MST is N-1 edges
      int *merged = calloc(3 * MAX_EDGES, sizeof(int));   // MST is N-1 edges
      int *to_merge = calloc(3 * MAX_EDGES, sizeof(int)); // MST is N-1 edges
      int edge_count = 0, next_count;
      // todo: speed up with memcopy?

      for (int proc = 0; proc < nb_procs; proc++)
      {
        int start = proc * (3 * MAX_EDGES);
        for (int i = 0; i < 3 * counts[proc] + 3; i++)
        {
          to_merge[i] = edge_buf[i + start];
        }
        next_count = counts[proc];
        // merge all procs into root, and do kruskal
        // swap buffers between two processes
        edge_count = merge(out, merged, to_merge, edge_count, next_count);
        edge_count = kruskal(merged, edge_count, out, N, NULL);
      }
      output_edges(merged, edge_count);
    }
  }
  else
  { // Invalid algorithm name
    if (proc_rank == 0)
    {
      printf("ERROR: Invalid algorithm name: %s.\n", algo_name);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
