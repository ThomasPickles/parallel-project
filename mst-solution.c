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
  parents_array[l1] = l2;
}
void print_array(int *arr, int count, int log_level)
{
  if (log_level > 0)
  {
    printf("[");
    for (int i = 0; i < count - 1; i++)
    {
      printf("%d,", arr[i]);
    }
    printf("%d]\n", arr[count - 1]);
  }
}
void add_edge(int e1, int e2)
{
  printf("%d %d\n", e1, e2);
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
int get_global_id(int local_id, int vertices_per_proc, int proc_rank)
{
  return local_id + vertices_per_proc * proc_rank;
}

int not_in_tree(int local_id, int vertices_per_proc, int proc_rank, int *added)
{
  int id = get_global_id(local_id, vertices_per_proc, proc_rank);
  return *(added + id) == 0;
}
void compute_mst(
    int N,
    int M,
    int *adj,
    char *algo_name)
{
  int VERBOSE = 0;
  int proc_rank = 0, nb_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);

  int lowest_weight[5] = {10, 10, 10, 10, 10};
  int vertices[5] = {4, 0, 3, 1, 2}; // Start with 4
  int v, min_key, weight;
  int abs_lowest = 100;
  int *last_added = vertices, *ptr = NULL, *back = NULL;
  int mst = 0;
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

    for (int vertex_count = 1; vertex_count < 5; vertex_count++)
    {

      printf("Starting array [");
      for (int i = 0; i < 5; i++)
      {
        printf("%d,", vertices[i]);
      }
      printf("]\n");
      printf("Considering edges from vertex %d\n", *last_added);

      // pointer arithmetic allows us only to look at elements in the
      // candidate set, no need to look at edge weights for any vertices
      // that are not connected
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

    // To check for set membership for vertices, just keep an array of length N
    //  initialised to 0, and set value to 1 if vertex is added to set.
    int *added = calloc(N, sizeof(int));
    int *leaders = calloc(N, sizeof(int));
    int *v1 = calloc(M, sizeof(int));
    int *v2 = calloc(M, sizeof(int));
    int *weights = calloc(M, sizeof(int));
    int **loc = calloc(M, sizeof(int *)); // initialise array of pointers
    int w, t1, t2;
    int count = 0; // NOT THE SAME AS M, SINCE WE IGNORE SELF-LOOPS

    if (v1 == NULL)
    {
      fprintf(stderr, "malloc failed\n");
    }

    for (int i = 0; i < N; i++)
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

    if (loc == NULL)
    {
      fprintf(stderr, "malloc failed\n");
    }

    if (count > 10)
    {
      printf("Array too small, aborting...");
      return;
    }

    print_array(weights, count, VERBOSE);

    int idx;
    int i;
    for (i = 0; i < count; i++)
    {
      loc[i] = &(weights[i]); // edge address
    }

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

    // loc is pointer to first elem
    qsort(loc, count, sizeof loc, deref_pointer);

    if (VERBOSE > 0)
    {
      puts("After sorting:");
      for (i = 0; i < count; i++)
      {
        printf("value: %-6d posn: %d\n", *loc[i], (int)(loc[i] - weights));
      }
    }

    for (i = 0; i < count; i++)
    {
      idx = loc[i] - weights;
      w = weights[idx];
      t1 = v1[idx]; // t1 < t2
      t2 = v2[idx];

      if (VERBOSE > 0)
      {
        printf("Considering vertex of weight %d from %d to %d\n", w, t1, t2);
      }

      int l1, l2;
      l1 = get_leader(leaders, t1);
      l2 = get_leader(leaders, t2);

      if (l1 != l2)
      {
        // TODO: if both t1 and t2 already in set, no need to add edge as it's already
        // spanned.  If either of t1/t2 in set, edd edge and do a set union with other set
        // If neither in set yet, we create a new set which includes the two vertices
        mst += w;
        added[t1] = 1;
        added[t2] = 1;
        make_set(leaders, t1, t2);
        add_edge(t1, t2);

        if (added[t1] == 1 && added[t2] == 1)
        {
          if (VERBOSE > 0)
          {
            printf(" Neither %d or %d added, adding edge", t1, t2);
          }
        }
        else
        {

          // add edge to set
          if (VERBOSE > 0)
          {
            printf(" Disjoint sets, adding edge from %d to %d with leaders %d and %d\n", t1, t2, l1, l2);
          }
          // should this be other way around too?
          // nope, need to point the leader of t2 to be l1
        }
      }
      else
      {
        if (VERBOSE > 0)
        {
          printf(" Edge from %d to %d in same set with leaders %d and %d.  Discarding...\n", t1, t2, l1, l2);
        }
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
    if (N % nb_procs != 0)
    {
      printf("Procs doesn't exactly divide vertices, exiting...");
      return;
    }
    int root = 0;
    int nb_loc = N / nb_procs;
    int start_vertex = 1; // less than nb_loc
    int *added = calloc(N, sizeof(int));
    int *min_weight = calloc(nb_loc, sizeof(int));
    int *key = calloc(nb_loc, sizeof(int));
    int *gather_mins = calloc(3 * nb_procs, sizeof(int));
    int proc_min[3] = {0};
    int i, edge, next, v_i; //

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

          edge = adj[v_j + i * N];
          if (edge == 0)
            continue; // only want edges

          if (edge < min_weight[i] || min_weight[i] == 0) // and edge is lower weight
          {
            // update its key
            min_weight[i] = edge;
            key[i] = v_j;
            // TODO: feels costly to check whether we've set the proc_min
            // to zero every time, since we just want to do it once
            // for the first edge we find
            // then, if it's the lowest weight yet seen by process
            if (edge < proc_min[0] || proc_min[0] == 0)
            {
              // printf("lower weight edge R%dP%d %d:%d-%d\n", vertex_count, proc_rank, edge, v_i, v_j);
              proc_min[0] = edge;
              // Put into lexicographic order
              proc_min[1] = v_i < v_j ? v_i : v_j;
              proc_min[2] = v_i < v_j ? v_j : v_i;
            }
          }
        }
      }
      // printf("\nP%d: %d:%d-%d\n", proc_rank, proc_min[0], proc_min[1], proc_min[2]);
      // for (i = 0; i < nb_loc; i++)
      // {
      //   printf("%d ", min_weight[i]);
      // }
      // send mins back to root
      // TODO: can just do a reduce here, as we can discard any
      // values that are less than the minimum
      MPI_Gather(&proc_min, 3, MPI_INT, gather_mins, 3, MPI_INT, root, MPI_COMM_WORLD);

      if (proc_rank == root)
      {
        // printf("Round 1, Expect \n3 0 1 1 1 3 5 1 5\n");
        // for (int i = 0; i < 3 * nb_procs; i++)
        // {
        //   printf("%d ", gather_mins[i]);
        // }
        // todo: take this decl out of loop
        int min[3] = {0, 0, 0};
        for (int i = 0; i < nb_procs; i++)
        {
          if (gather_mins[3 * i] > 0)
          {
            if (min[0] == 0 ||
                (gather_mins[3 * i] < min[0]) ||
                ((gather_mins[3 * i] == min[0]) && (gather_mins[3 * i + 1] < min[1])) ||
                ((gather_mins[3 * i] == min[0]) && (gather_mins[3 * i + 1] == min[1]) && (gather_mins[3 * i + 2] < min[2])))
            {
              // printf("Edge %d %d, w %d \n", next, gather_mins[1], gather_mins[2]);

              min[0] = gather_mins[3 * i];
              min[1] = gather_mins[3 * i + 1];
              min[2] = gather_mins[3 * i + 2];
            }
          }
        }
        if (added[min[1]])
        {
          next = min[2];
        }
        else
        {
          next = min[1];
        }
        printf("%d %d\n", min[1], min[2]);
      }
    }
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
