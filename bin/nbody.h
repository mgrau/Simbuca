#ifdef __NBODY_ON__
#ifndef NBODY_H_INCLUDED
#define NBODY_H_INCLUDED
#include "ioncloud.h"

#ifdef __MPI_ON__
#include "mpi_funcs.h"
#endif // __MPI_ON__

#include "cuda.h"
#include "cuda_runtime_api.h"
#include <iostream>
#include <fstream>
#include <string> 
#include <vector>


using namespace std;
class _nbody{
public:
    //structor
    _nbody();
    ~_nbody();
    //function
    void Initialization(const IonCloud & _cloud, bool b_);
    void Finalize();
    int FindSupN(int n_);
    void LoadFile();    
    
    void LoadPosition(const IonCloud &_cloud);
    void SetParameters();
    void SetParameterWithoutFile();
    void SetParameterWithFile();
    void force_begin(const IonCloud & _cloud);
    void force_end(const IonCloud &_cloud);
    

    //variable
    float4 *host_pos;
    float3 *host_acc;

    int n_max_alloc;
    int n_max_file;
    int n_gpu;
    int n_cpu;
    int p_gpu;
    int q_gpu;
    int t_gpu;
    int sm;
    int cc;
    int gpu_name;
    int step;
    ifstream GPUfile;
    bool file;
    string name;
    std::vector< int > Vgpu_n;
    std::vector< int > Vgpu_p;
    std::vector< int > Vgpu_q;
    std::vector< int > Vgpu_ntiles;
    bool withCharge;
    cudaEvent_t end_gpu2;

#ifdef __MPI_ON__
    _mpi_vars * mpiv;
    void InitMPI(_mpi_vars * _mpiv);
    float4 *host_pos_sub;
    float3 *host_acc_sub;
#endif // __MPI_ON__
private:
    
};






#endif // NBODY_H_INCLUDED
#endif // __NBODY_ON__
