/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algo_name the name of the algorithm to be executed
 */

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
  // TODO: CAN MAKE MORE EFFICIENT BY UPDATING POINTER
  // EACH TIME WE GET LEADER.  MAKES TREE MUCH SHALLOWER
  int parent = parent_array[elem];
  if (parent == elem)
  {
    return elem;
  }
  else
  {
    return find_set(parent_array, parent);
  }
}
void union_sets(int *parents_array, int t1, int t2)
{
  // usually have t1 < t2, so join sets to lower numbers
  int l1 = find_set(parents_array, t1);
  int l2 = find_set(parents_array, t2);
  parents_array[l2] = l1;
}
void print_array(int *arr, int count, int log_level)
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
void add_edge(int e1, int e2)
{
  int min, max;
  min = e1 < e2 ? e1 : e2;
  max = e1 < e2 ? e2 : e1;
  printf("%d %d\n", min, max);
}

int deref_pointer(const void *v1, const void *v2)
{
  const int i1 = **(const int **)v1;
  const int i2 = **(const int **)v2;
  return i1 < i2 ? -1 : (i1 > i2);
}

void sort_array(int **loc, int count)
{
  qsort(loc, count, sizeof loc, deref_pointer);
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

int VERBOSE = 0;
void create_edge_list(int *ret, int *edges, int **loc, int *adj, int N)
{
  // adjancency matrix is symmetric, so only consider upper
  // triangular part
  int *ptr = edges; //, *ptr1 = v1, *ptr2 = v2;
  int i, j, w;
  int count = 0;
  for (i = 0; i < N; i++)
  {
    for (j = i + 1; j < N; j++)
    {
      w = adj[i * N + j];
      if (has_edge(w))
      {
        // update and increment
        loc[count] = ptr; // edge address
        *ptr++ = w;       // weights
        *ptr++ = i;       // v1
        *ptr++ = j;       // v2
        count++;
      }
    }
  }

  print_array(edges, 3 * count, 0);

  if (VERBOSE > 0)
  {
    puts("Before sorting:");
    for (i = 0; i < count; i++)
    {
      printf("value: %-6d adress: %p\n", *loc[i], loc[i]);
    }
  }
  // read through code
  // pointer to pointer is useful to give location of a result,
  // since we can only pass by value

  // edge order
  sort_array(loc, count);

  if (VERBOSE > 0)
  {
    puts("After sorting:");
    for (i = 0; i < count; i++)
    {
      printf("value: %-6d posn: %d\n", *loc[i], (int)(loc[i] - edges));
    }
  }

  *ret = count;
}

void kruskal(int *out_buf, int edge_count, int *edges, int **loc, int N, int *vertices_to_include)
{
  int *leaders = calloc(N, sizeof(int));
  int idx, i, w;
  int mst = 0;

  for (int i = 0; i < N; i++)
  {
    leaders[i] = i;
  }
  int *added = calloc(N, sizeof(int));
  int t1, t2, l1, l2;
  for (i = 0; i < edge_count; i++)
  {
    idx = loc[i] - edges;
    w = edges[idx];
    t1 = edges[idx + 1]; // t1 < t2
    t2 = edges[idx + 2];

    l1 = find_set(leaders, t1);
    l2 = find_set(leaders, t2);

    if (l1 != l2)
    {
      // different set, add edge to combine sets
      add_edge(t1, t2);
      // combine sets
      union_sets(leaders, t1, t2);
      // add vertices
      added[t1] = 1;
      added[t2] = 1;
      mst += w;
    }
    else
    {
      // Edge from %d to %d in same set with leader %d.  Ignore
    }
  }

  free(added);
  free(leaders);
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

  add_edge(u_min, v_min);

  // remove edge as candidate
  vertices[3 * v_min] = 2;
  vertices[3 * v_min + 1] = 0;
  *last_added = v_min;
}

int get_global_id(int local_id, int vertices_per_proc, int proc_rank)
{
  return local_id + vertices_per_proc * proc_rank;
}

void compute_mst(
    int N,
    int M,
    int *adj,
    char *algo_name)
{
  int proc_rank = 0, nb_procs = 0;
  int mst = 0;
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
    int start_vertex = 0;
    int *vertices = calloc(3 * N, sizeof(int));
    int last_added = start_vertex;

    // (added, min_weight_edge, to_vertex)
    for (int i = 0; i < N; i++)
    {
      // 0 - unseen; 1 - adjacent; 2 - added
      vertices[3 * i] = 0;      // not added
      vertices[3 * i + 1] = 0;  // min weight
      vertices[3 * i + 2] = -1; // key
    }

    vertices[2 * last_added] = 2;

    for (int vertex_count = 1; vertex_count < N; vertex_count++)
    {

      update_mins_for_vertices(vertices, adj, last_added, N);
      get_abs_min(&last_added, vertices, N);
    }
  }
  else if (strcmp(algo_name, "kruskal-seq") == 0)
  { // Sequential Kruskal's algorithm
    if (proc_rank == 0 && nb_procs != 1)
    {
      printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", nb_procs);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // BEGIN IMPLEMENTATION HERE
    // TODO: HOW TO WRITE TESTS IN C?
    // TODO: HOW TO PASS ARRAY SIZE AS A PARAM?
    // TODO: DO I WANT TO USE CONST FOR SOME OF THESE ARRAY POINTERS?

    // To check for set membership for vertices, just keep an array of length N
    //  initialised to 0, and set value to 1 if vertex is added to set.

    int *edges = calloc(3 * M, sizeof(int));
    int **loc = calloc(M, sizeof(int *)); // initialise array of pointers
    int count = 0;                        // NOT THE SAME AS M, SINCE WE IGNORE SELF-LOOPS

    // todo: check all pointers are allocated properly
    if (edges == NULL)
    {
      fprintf(stderr, "malloc failed\n");
    }

    create_edge_list(&count, edges, loc, adj, N);

    int *vertices_to_include = calloc(N, sizeof(int));

    // do for all vertices
    for (int i = 0; i < N; i++)
    {
      vertices_to_include[i] = 1;
    }

    kruskal(NULL, count, edges, loc, N, vertices_to_include);

    free(edges);
    free(loc);

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
    int nb_loc = N / nb_procs;
    int start_vertex = 0; // less than nb_loc
    int *added = calloc(N, sizeof(int));
    int *min_weight = calloc(nb_loc, sizeof(int));

    int *vertices = calloc(3 * nb_procs, sizeof(int));
    int proc_min[3] = {0};
    int i, w, next, v_i; //

    if (proc_rank == root)
    {
      next = start_vertex;
    }

    for (int vertex_count = 0; vertex_count < N - 1; vertex_count++)
    {
      MPI_Bcast(&next, 1, MPI_INT, root, MPI_COMM_WORLD);

      // add to set
      // printf("Adding %d to set\n", next);
      added[next] = 1;

      // get minimum for process, do on all processes
      // I'm concerned that we can't just keep the minimum, since when this
      // gets sucked up into the tree, how do we know what is the second
      // smallest? For now we can just do the double for-loop, and perhaps
      // look at some optimatisions later (make it run, make it right, make it fast)
      proc_min[0] = 0;
      proc_min[1] = 0;
      proc_min[2] = 0;

      for (i = 0; i < nb_loc; i++)
      {
        // yuck, this is the slow step, we're throwing away all the previous
        // calculations.  can we keep a heap somehow?
        min_weight[i] = 0;

        v_i = get_global_id(i, nb_loc, proc_rank);

        if (added[v_i] == 1)
        {
          min_weight[i] = 0;
          continue;
        }

        for (int v_j = 0; v_j < N; v_j++)
        {

          if (added[v_j] == 0)
          {
            continue;
          }

          if (v_i == v_j)
            // don't include self loops,
            continue;

          w = adj[v_j + i * N];
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
              // printf("lower weight edge R%dP%d %d:%d-%d\n", vertex_count, proc_rank, edge, v_i, v_j);
              proc_min[0] = w;
              // Put into lexicographic order
              proc_min[1] = v_i < v_j ? v_i : v_j;
              proc_min[2] = v_i < v_j ? v_j : v_i;
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
        add_edge(u_min, v_min);
      }
    }
  }
  else if (strcmp(algo_name, "kruskal-par") == 0)
  { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE
    // TEST merge
    // N = 5, M=5
    // int V[5] = {0, 1, 2, 3, 4};
    int e[15] = {0, 1, 2, 0, 3, 1, 1, 2, 3, 1, 3, 4, 3, 4, 2}; // tuple (v1,v2,w)
    // in-place modification of e:
    // kruskal(e, 5, V, 5);
    for (int i = 0; i < 15; i++)
    {
      printf("%d ", e[i]);
    }
    // V1 = [ 1, 1, 1, 0, 0 ]; // vertices {0,1,2};
    // V2 = [ 0, 1, 0, 1, 1 ]; // vertices {1,3,4};
    // e1 = {0, 1, 2, 1, 2, 3};
    // e2 = {1, 3, 4, 3, 4, 2};
    // // output
    // V = [ 1, 1, 1, 1, 1 ];
    // e = {0, 1, 2, 0, 3, 1, 1, 2, 3, 3, 4, 2};
    // merge might be easier than expected, just need to run Kruskal
    // if intersection (V1, V2)
    //   is empty : add lowest weight u, v edge with u in V1 and v in V2 if | intersection(V1, V2) > 1 | : remove common edges
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
