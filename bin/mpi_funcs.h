
#ifndef MPI_FUNCS_H_INCLUDED
#define MPI_FUNCS_H_INCLUDED


#ifdef __MPI_ON__
#include <vector>
#include "ion.h"
#include "mpi.h"
#include "SLogger.h"

class _mpi_vars {
public:
    // structors
    _mpi_vars();
    ~_mpi_vars();
    //functions
    void ComputeCounts_Disps(int nrparticles);
#ifdef __CPUNBODY_ON__
    void Init_GlobalPos(int nrparticles);
#endif // CPU
    //variables
    
    int nb_procs; // Number of processes
    int rank;     // Rank of the process in the comm2d communicator
    MPI_Comm comm2d; // Cartesian 2D-topology communicator
    int *counts, *disps;
    int *counts_3, *disps_3;
    int nb_gpu;   // number of GPUs;
    int nb_node_per_gpu; // number of node per GPU;
    int n_part_tot;
    int * coords;
#if defined(__CUNBODY_ON__) || defined(__NBODY_ON__) || defined(__NBODY2_ON__)

    int * ID_node_GPU ; // MPI rank of node using GPU
    MPI_Comm comm1d_h; // communicator_row
    MPI_Comm comm_gpu;
    double *charge; // charge of ions
    int *counts_F4, *disps_F4;
    int *counts_F3, *disps_F3;
    int *counts_col,*disps_col;
    int *counts_row,*disps_row,*disps_row_return; 
    int *counts_Arow,*disps_Arow; // counts and disps total number of part in each row dim = nb_gpu
    MPI_Datatype MPI_FLOAT3;
    MPI_Datatype MPI_FLOAT4;
    void LoadCharge( vector<int > ion_charge);
#endif // GPU
    
#ifdef __CPUNBODY_ON__
    double (*pos_tot)[3]; // whole cloud
#endif // CPU
    
private:
    
};





#endif //__MPI_ON__
#endif // MPI_FUNCS_H_INCLUDED
