#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#include "mst-solution.c"

void read_graph(
        char *file_name,
        int *ptr_nb_vertex,
        int *ptr_nb_edge,
        int **ptr_adj){

    FILE *file = fopen(file_name, "r");
    int rank = 0,size = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (file == NULL && rank == 0) {
        printf("ERROR: Unable to open the file %s.\n", file_name);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int nb_vertex = 0, nb_edge = 0;
    int *adj = 0;
    int *temptr_adj = 0;

    if(!fscanf(file, " %d %d", ptr_nb_vertex, ptr_nb_edge)){
        printf("error while reading file!\n");
        exit(-1);
    }

    nb_vertex = *ptr_nb_vertex;
    nb_edge = *ptr_nb_edge;
    int *edge_left = NULL, *edge_right = NULL, *weight_l = NULL;

    if (rank == 0){
        edge_left = (int *) malloc(nb_edge * sizeof(edge_left[0]));
        edge_right = (int *) malloc(nb_edge * sizeof(edge_right[0]));
        weight_l = (int *) malloc(nb_edge * sizeof(weight_l[0]));
    }
    temptr_adj = (int *) malloc(nb_vertex * nb_vertex * sizeof(temptr_adj[0]));

    //Allocate the matrix and fill it with input data
    int nb_elements = nb_vertex*(int)ceil((float)nb_vertex/(float)size);
    if (rank == size-1){
        *ptr_adj = adj = (int *) malloc(
                (nb_vertex*nb_vertex-(size-1)*nb_elements)*sizeof(adj[0]));
    }else{
        *ptr_adj = adj = (int *) malloc(
                nb_elements*sizeof(adj[0]));
    }

    if (rank == 0){
        int i,j;
        for (i = 0; i < nb_vertex; i++){
            for (j = 0; j < nb_vertex; j++){
                temptr_adj[i*nb_vertex+j] = 0;
            }
        }

        for (i = 0; i < nb_edge; i++) {
            if(!fscanf(file, " %d %d %d", edge_left+i, edge_right+i, weight_l+i)){
                printf("error while reading file!\n");
                break;
            }
            temptr_adj[edge_left[i]*nb_vertex + edge_right[i]] = weight_l[i];
            temptr_adj[edge_right[i]*nb_vertex + edge_left[i]] = weight_l[i];
        }

        free(edge_left);
        free(edge_right);
        free(weight_l);

    }

    if (size > 1){
        //  MPI_Bcast(adj, nb_vertex*nb_vertex, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Comm comm_first;
        MPI_Comm_split(MPI_COMM_WORLD, rank/(size-1), rank, &comm_first);

        //SCATTER THE MATRIX FOR P-1 PROCESS
        if (rank != size-1)
            MPI_Scatter(temptr_adj, nb_elements, MPI_INT, adj, nb_elements,
                    MPI_INT, 0, comm_first);

        //SEND THE LAST ELEMENTS TO THE LAST PROCESS
        if (rank == 0){
            MPI_Send(temptr_adj+(size-1)*nb_elements,
                    nb_vertex*nb_vertex-(size-1)*nb_elements,
                    MPI_INT, size-1, 0, MPI_COMM_WORLD);
        }else if (rank == size-1){
            MPI_Recv(adj, nb_vertex*nb_vertex - (size-1)*nb_elements,
                    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }else{
        memcpy(adj, temptr_adj, nb_vertex*nb_vertex*sizeof(temptr_adj[0]));
    }

    free(temptr_adj);

    fclose(file);
}

void printUsage() {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf("Usage: mpirun -np [num-procs] ./mst [graph-file-name] [algo-name]\n"
                "Arguments:\n"
                "\t[num-procs]: The number of MPI ranks to be used. For the sequential algorithm it has to be 1.\n"
                "\t[graph-file-name]: Name of the graph file.\n"
                "\t[bisect-algo-name]: Name of the graph bisection algorithm. There are three possibilities:\n"
                "\t\tprim-seq: Sequential Prim's algorithm.\n"
                "\t\tkruskal-seq: Sequential Kruskal's algorithm.\n"
                "\t\tprim-par: Parallel Prim's algorithm.\n"
                "\t\tkruskal-par: Parallel Kruskal's algorithm.\n"
              );
    }
}


int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    if (argc < 3) {
        printUsage();
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int proc_rank=0, size=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char *algo_name = argv[2];
    int nb_vertex = 0, nb_edge = 0;
    int *adj = NULL;


    read_graph(argv[1], &nb_vertex, &nb_edge, &adj);
    MPI_Barrier(MPI_COMM_WORLD);

    double start_time = MPI_Wtime();
    compute_mst(nb_vertex, nb_edge, adj, algo_name);

    if (proc_rank == 0) {
        printf("compute_mst took %e seconds.\n", MPI_Wtime() - start_time);
    }

    MPI_Finalize();

    free(adj);
    return 0;
}
