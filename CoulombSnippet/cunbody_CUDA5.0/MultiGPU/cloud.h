
#ifndef CLOUD_H_INCLUDED
#define CLOUD_H_INCLUDED

#include "mpi.h"



struct _cloud {
    int elements; // total number of particles
    int n_sub; // Number of particles in each mpi process
    int n_sub_gpu;// Number of bodies in gpu sum 
    double (*xi)[3]; // position array (total cloud)
    double (*sub_xi)[3]; // position array (MPI sub cloud)
    double (*xj)[3]; // position array (gpu sub cloud)
    double *mj; // mass
    //double sub_ai[numprocs_row*n_sub][3]; // acc array (gpu sub cloud)
    double (*sub_aj)[3]; // acc array (MPI sub cloud)
    double (*ai)[3]; // acc array (total cloud)
    double (*ai_temp)[3]; // acc array (total cloud)
};




void Initialization(int N_);

#endif // MPI_FUNCS_H_INCLUDED
