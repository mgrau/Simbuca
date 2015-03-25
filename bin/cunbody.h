#ifdef __CUNBODY_ON__
#ifndef CUNBODY_H_INCLUDED
#define CUNBODY_H_INCLUDED

#include "ioncloud.h"
#ifdef __MPI_ON__
#include "mpi_funcs.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
#endif // __MPI_ON__
extern "C" {
    void cunbody1_force(double xj[][3], // position of j-th particles
            double mj[],    // mass of j-th particles
            double xi[][3], // position of i-th particles
            double eps2,    // softening parameter
            double ai[][3], // force of i-th particles
            int ni,         // number of i-th particles
            int nj);        // number of j-th particles
    void cunbody1_force_end(double xj[][3], // position of j-th particles
            double mj[],    // mass of j-th particles
            double xi[][3], // position of i-th particles
            double eps2,    // softening parameter
            double ai[][3], // force of i-th particles
            int ni,         // number of i-th particles
            int nj);        // number of j-th particles
    void cunbody1_force_begin(double xj[][3], // position of j-th particles
            double mj[],    // mass of j-th particles
            double xi[][3], // position of i-th particles
            double eps2,    // softening parameter
            double ai[][3], // force of i-th particles
            int ni,         // number of i-th particles
            int nj);        // number of j-th particles
}

class _cunbody{

    public:
        //structors
        _cunbody();
        ~_cunbody();
        //functions
        void begin(const IonCloud & _cloud);
        void end(const IonCloud & _cloud);
        void AllocArrays(int nrparticles);
        //variables
        double (*ai)[3];
        double (*pos_tmp)[3];
        double (*m);
        float eps;
        int coords0;
        // MPI

#ifdef __MPI_ON__
        _mpi_vars * mpiv;
        void InitMPI(_mpi_vars * _mpiv);
        double (*pos_tot)[3];
        double (*ai_temp)[3];
#endif // __MPI_ON__
    private:
};

#ifdef after
void MultiGPUs_CUNBODY_begin();
void MultiGPUs_CUNBODY_end();
void CompaCPU();
#endif

#endif //  CUNBODY_H_INCLUDED
#endif // __CUNBODY_ON__
