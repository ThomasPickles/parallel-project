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

int get_leader(int *parent_array, int elem)
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
    return get_leader(parent_array, parent);
  }
}
void make_set(int *parents_array, int t1, int t2)
{
  int l1 = get_leader(parents_array, t1);
  int l2 = get_leader(parents_array, t2);
  printf(" Setting leaders of %d and %d to be %d\n", t1, t2, l2);
  parents_array[l1] = l2;
}

int deref_pointer(const void *v1, const void *v2)
{
  const int i1 = **(const int **)v1;
  const int i2 = **(const int **)v2;
  return i1 < i2 ? -1 : (i1 > i2);
}

void sort_array(int **buf, int *arr, int count)
{

  int i;
  printf("Array of values:\n ");
  for (i = 0; i < count; i++)
    printf("%d ", arr[i]);

  printf("Array of pointers:\n ");
  for (i = 0; i < count; i++)
  {
    printf("%p ", buf[i]);
    buf[i] = &arr[i];
  }
  // qsort(buf, count, sizeof *buf, deref_pointer);

  // for (i = 0; i < count; i++)
  //   printf("%d ", arr[i]);
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

  int lowest_weight[5] = {10, 10, 10, 10, 10};
  int vertices[5] = {4, 0, 3, 1, 2}; // Start with 4
  int v, min_key, weight;
  int abs_lowest = 100;
  int *last_added = vertices, *ptr = NULL, *back = NULL;
  int mst = 0;

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

    for (int vertex_count = 1; vertex_count < 5; vertex_count++)
    {

      printf("Starting array [");
      for (int i = 0; i < 5; i++)
      {
        printf("%d,", vertices[i]);
      }
      printf("]\n");
      printf("Considering edges from vertex %d\n", *last_added);

      ptr = last_added + 1;
      back = vertices + 5; // beyond last elem
      while (ptr != back)
      {
        v = *ptr;
        weight = adj[*last_added + v * N];

        printf(" Edge weight for vertex (%d,%d) is %d\n", *last_added, v, weight);

        if (previous_candidate(lowest_weight[v]))
        {
          printf("  Vertex %d has been a candidate before with weight %d\n", v, lowest_weight[v]);
          if (has_edge(weight) && weight < lowest_weight[v])
          {
            printf("   Found an edge for %d with weight %d\n", v, weight);
            // we've found a lower weight edge, update it
            // update weight
            // 0: inf -> 1; 1: inf -> 5; 2: inf -> 1
            printf("   Update weight for vertex %d: %d -> %d\n", v, lowest_weight[v], weight);
            lowest_weight[v] = weight;
          }
          ptr++;
        }
        else
        {

          printf(" Vertex %d not previously a candidate for inclusion...\n", v);
          if (has_edge(weight))
          {
            printf("  Vertex %d now a candidate with edge weight %d\n", v, weight);
            // it's now a candidate, set its value
            lowest_weight[v] = weight;
            ptr++;
          }
          else
          {
            printf("  Vertex %d still not a candidate\n", v);
            // it's still not a candidate, swap with unseen elem at back
            // swap elems, but don't disturb pointers
            back--;
            printf("  Elem %d to back and %d to front\n", v, *back);
            *ptr = *back;
            *back = v;
          }
        }
      }

      printf("Candidate edges from all added %d\n", *last_added);

      ptr = last_added + 1;
      min_key = *ptr;
      abs_lowest = 10;
      while (ptr != back)
      {
        v = *ptr;
        printf(" Candidate %d, weight %d\n", v, lowest_weight[v]);

        if (lowest_weight[v] < abs_lowest)
        {
          min_key = v;
          abs_lowest = lowest_weight[v];
        }
        ptr++;
      }

      printf("Lowest weight candidate %d with weight %d, adding to set\n", min_key, abs_lowest);

      last_added++;
      *last_added = min_key;
      mst += abs_lowest;
      printf("Value of mst is %d\n\n", mst);
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
    int leaders[5];
    int *v1 = calloc(M, sizeof(int));
    int *v2 = calloc(M, sizeof(int));
    int *weights = calloc(M, sizeof(int));
    int w, t1, t2;
    int count = 0; // NOT THE SAME AS M, SINCE WE IGNORE SELF-LOOPS

    if (v1 == NULL)
    {
      fprintf(stderr, "malloc failed\n");
    }

    for (int i = 0; i < 5; i++)
    {
      // initialize leaders to point to self
      leaders[i] = i;
    }

    ptr = weights;
    int *ptr1 = v1, *ptr2 = v2;
    // adjancency matrix is symmetric, so only consider upper
    // triangular part
    for (int i = 0; i < N; i++)
    {
      for (int j = i + 1; j < N; j++)
      {
        weight = adj[i * N + j];
        if (has_edge(weight))
        {
          count++;
          // update and increment
          *ptr++ = weight;
          *ptr1++ = i;
          *ptr2++ = j;
        }
      }
    }

    // edge order
    // sort_array(&order, weights, 10);
    int **loc = calloc(count, sizeof(int));

    if (loc == NULL)
    {
      fprintf(stderr, "malloc failed\n");
    }

    if (count > 10)
    {
      printf("Array too small, aborting...");
      return;
    }

    printf("Starting array [");
    for (int i = 0; i < count; i++)
    {
      printf("%d,", weights[i]);
    }
    printf("]\n");

    int idx;
    int i;

    puts("Before sorting:");
    for (i = 0; i < count; i++)
    {
      loc[i] = &(weights[i]); // edge address
      printf("value: %-6d adress: %p\n", *loc[i], loc[i]);
    }

    // read through code
    // pointer to pointer is useful to give location of a result,
    // since we can only pass by value

    // loc is pointer to first elem
    qsort(loc, count, sizeof loc, deref_pointer);

    puts("After sorting:");
    for (i = 0; i < count; i++)
    {
      printf("value: %-6d posn: %d\n", *loc[i], (int)(loc[i] - weights));
    }

    for (i = 0; i < count; i++)
    {
      idx = loc[i] - weights;
      w = weights[idx];
      t1 = v1[idx]; // t1 < t2
      t2 = v2[idx];

      printf("Considering vertex of weight %d from %d to %d\n", w, t1, t2);

      int l1, l2;
      l1 = get_leader(leaders, t1);
      l2 = get_leader(leaders, t2);

      if (l1 != l2)
      {
        // add edge to set
        mst += w;
        printf(" Disjoint sets, adding edge from %d to %d with leaders %d and %d\n", t1, t2, l1, l2);
        // should this be other way around too?
        // nope, need to point the leader of t2 to be l1
        make_set(leaders, t1, t2);
      }
      else
      {
        printf(" Edge from %d to %d in same set with leaders %d and %d.  Discarding...\n", t1, t2, l1, l2);
      }
    }
    printf("MST has weight %d\n", mst);
    free(v1);
    free(v2);
    free(weights);
    free(loc);

    //////
  }
  else if (strcmp(algo_name, "prim-par") == 0)
  { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE
  }
  else if (strcmp(algo_name, "kruskal-par") == 0)
  { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE
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
