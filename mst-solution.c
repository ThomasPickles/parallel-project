/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algo_name the name of the algorithm to be executed
 */
void compute_mst(
    int N,
    int M,
    int *adj,
    char *algo_name)
{
  int proc_rank = 0, nb_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);

  if (strcmp(algo_name, "prim-seq") == 0) { // Sequential Prim's algorithm
    if (proc_rank == 0 && nb_procs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", nb_procs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algo_name, "kruskal-seq") == 0) { // Sequential Kruskal's algorithm
    if (proc_rank == 0 && nb_procs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", nb_procs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algo_name, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algo_name, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else { // Invalid algorithm name
    if (proc_rank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algo_name);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
