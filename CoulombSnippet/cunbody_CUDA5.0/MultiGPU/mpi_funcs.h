
#ifndef MPI_FUNCS_H_INCLUDED
#define MPI_FUNCS_H_INCLUDED

#include "mpi.h"



struct _mpi_vars {
  int nb_procs; // Number of processes
  int rank;     // Rank of the process in the comm2d communicator
  int nb_gpu;   // number of GPUs;
  int nb_node_per_gpu; // number of node per GPU;
  int * coords;
  int * ID_node_GPU ; // MPI rank of node using GPU
  MPI_Comm comm2d; // Cartesian 2D-topology communicator
  MPI_Comm comm1d_h; // communicator_row
  MPI_Comm comm_gpu;
 // MPI_Group server_grp; // all the nodes in the server (intra-nodes)
 // MPI_Group gpu_grp; // all the nodes handling gpu
    int *counts, *disps;
    int *counts_col,*disps_col;
    int *counts_row,*disps_row,*disps_row_return;
};




void init_mpi(int argc, char *argv[]);

void clean_mpi();

#endif // MPI_FUNCS_H_INCLUDED
