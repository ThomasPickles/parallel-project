/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algo_name the name of the algorithm to be executed
 */

int has_edge(int edge_weight)
{
  return edge_weight > 0;
}
int previous_candidate(int min_weight)
{
  return min_weight != 10;
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

  if (strcmp(algo_name, "prim-seq") == 0)
  { // Sequential Prim's algorithm
    if (proc_rank == 0 && nb_procs != 1)
    {
      printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", nb_procs);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // BEGIN IMPLEMENTATION HERE

    int lowest_weight[5] = {10, 10, 10, 10, 10};
    int vertices[5] = {4, 0, 3, 1, 2}; // Start with 4
    int v, min_key, edge_weight;
    int abs_lowest = 100;
    int *last_added = vertices, *ptr = NULL, *back = NULL;
    int mst = 0;

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
        edge_weight = adj[*last_added + v * N];

        printf(" Edge weight for vertex (%d,%d) is %d\n", *last_added, v, edge_weight);

        if (previous_candidate(lowest_weight[v]))
        {
          printf("  Vertex %d has been a candidate before with weight %d\n", v, lowest_weight[v]);
          if (has_edge(edge_weight) && edge_weight < lowest_weight[v])
          {
            printf("   Found an edge for %d with weight %d\n", v, edge_weight);
            // we've found a lower weight edge, update it
            // update weight
            // 0: inf -> 1; 1: inf -> 5; 2: inf -> 1
            printf("   Update weight for vertex %d: %d -> %d\n", v, lowest_weight[v], edge_weight);
            lowest_weight[v] = edge_weight;
          }
          ptr++;
        }
        else
        {

          printf(" Vertex %d not previously a candidate for inclusion...\n", v);
          if (has_edge(edge_weight))
          {
            printf("  Vertex %d now a candidate with edge weight %d\n", v, edge_weight);
            // it's now a candidate, set its value
            lowest_weight[v] = edge_weight;
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

    // adjacency matrix size N**2
    // accordingly might not be feasible to convert to binary heap as this
    // requires adjacency list representation
  }
  else if (strcmp(algo_name, "kruskal-seq") == 0)
  { // Sequential Kruskal's algorithm
    if (proc_rank == 0 && nb_procs != 1)
    {
      printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", nb_procs);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // BEGIN IMPLEMENTATION HERE
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
