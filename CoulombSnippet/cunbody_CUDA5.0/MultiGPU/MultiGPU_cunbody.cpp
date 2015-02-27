#include "main.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
#include <iostream>
#include <math.h>
extern _cloud cloud;
extern _mpi_vars mpi_vars;

const double eps2 = 1.0e-18;

using namespace std;


void MultiGPUs_CUNBODY_begin()
{
    
    MPI_Allgather(&cloud.n_sub,1,MPI_INT,&mpi_vars.counts[0],1,MPI_INT,mpi_vars.comm2d);
       
    int ii;
    cloud.n_sub_gpu = 0;
    for(int i_=0;i_<mpi_vars.nb_node_per_gpu;i_++)
    {
        ii = i_ + mpi_vars.coords[0]*mpi_vars.nb_node_per_gpu;
        cloud.n_sub_gpu += mpi_vars.counts[ii];        
    }
    
    mpi_vars.counts[0] *= 3;
    mpi_vars.disps[0] =0;
    for(int i_=1;i_<mpi_vars.nb_procs;i_++)
    {
        mpi_vars.counts[i_] *= 3;
        mpi_vars.disps[i_] =  mpi_vars.disps[i_-1]+mpi_vars.counts[i_-1];
    }
    
    for(int i_=0;i_<mpi_vars.nb_node_per_gpu;i_++)
    {
        ii = i_ + mpi_vars.coords[0]*mpi_vars.nb_node_per_gpu;
        mpi_vars.counts_row[i_] = mpi_vars.counts[ii];
        if(i_==0)
        {
            mpi_vars.disps_row[0] = 0 ;
            mpi_vars.disps_row_return[0] = mpi_vars.disps[mpi_vars.ID_node_GPU[mpi_vars.coords[0]]];
       }
        else
        {
            mpi_vars.disps_row[i_] =  mpi_vars.disps_row[i_-1]+mpi_vars.counts_row[i_-1];
            mpi_vars.disps_row_return[i_] = mpi_vars.disps_row_return[i_-1]+mpi_vars.counts_row[i_-1];
        }
    }
   

    for(int dev=0;dev<mpi_vars.nb_gpu;dev++)
    {
        MPI_Gatherv(cloud.sub_xi[0], cloud.n_sub*3, MPI_DOUBLE, &cloud.xi[0][0], mpi_vars.counts, mpi_vars.disps, MPI_DOUBLE,mpi_vars.ID_node_GPU[dev],mpi_vars.comm2d);
    }
    MPI_Gatherv(cloud.sub_xi[0], cloud.n_sub*3, MPI_DOUBLE, &cloud.xj[0][0], mpi_vars.counts_row, mpi_vars.disps_row, MPI_DOUBLE,0,mpi_vars.comm1d_h);
    
    
    if(mpi_vars.coords[1]==0) // node using GPU
    {
        cudaSetDevice(mpi_vars.coords[0]);
        cunbody1_force_begin(cloud.xj, cloud.mj, cloud.xi, eps2, cloud.ai_temp, cloud.elements, cloud.n_sub_gpu );
    }
}

void MultiGPUs_CUNBODY_end()
{
    if(mpi_vars.coords[1]==0) // node using GPU
    {
        cunbody1_force_end(cloud.xj, cloud.mj, cloud.xi, eps2, cloud.ai_temp, cloud.elements, cloud.n_sub_gpu );
    }
    MPI_Allreduce(cloud.ai_temp[0],cloud.ai[0],cloud.elements*3,MPI_DOUBLE,MPI_SUM,mpi_vars.comm_gpu);
    MPI_Scatterv(&cloud.ai[0][0],mpi_vars.counts_row, mpi_vars.disps_row_return, MPI_DOUBLE,cloud.sub_aj[0],cloud.n_sub*3,MPI_DOUBLE,0,mpi_vars.comm1d_h);
    
}




void CompaCPU()
{
    double r_ij;
    double (*ai_cpu)[3] = new double [cloud.elements][3];
    
    for(int i=0; i<cloud.elements; i++)
    {
        for(int dim=0; dim<3; dim++)
            ai_cpu[i][dim]=0.0;
        
    }
    if(mpi_vars.rank==0)
    {
        for(int j_=0;j_ <cloud.elements;j_++)
        {
            for(int i_=0;i_ <cloud.elements;i_++)
            {
                if(i_ != j_)
                {
                    r_ij=sqrt((cloud.xi[i_][0]-cloud.xi[j_][0])*(cloud.xi[i_][0]-cloud.xi[j_][0])
                              +(cloud.xi[i_][1]-cloud.xi[j_][1])*(cloud.xi[i_][1]-cloud.xi[j_][1])
                              +(cloud.xi[i_][2]-cloud.xi[j_][2])*(cloud.xi[i_][2]-cloud.xi[j_][2]));
                    r_ij=fabs(1.0/(r_ij*r_ij*r_ij));
                   
                    ai_cpu[j_][0]-=(cloud.xi[j_][0]-cloud.xi[i_][0])*(r_ij);
                    ai_cpu[j_][1]-=(cloud.xi[j_][1]-cloud.xi[i_][1])*(r_ij);
                    ai_cpu[j_][2]-=(cloud.xi[j_][2]-cloud.xi[i_][2])*(r_ij);
                }
            }
        }
    }
    
    MPI_Bcast(&ai_cpu[0][0],cloud.elements*3,MPI_DOUBLE,0,mpi_vars.comm2d);
        int ii;

     
    double err_max = 0;
    double err_min = 1e8;
    double err_mean =0;
    double res;
    
    for(int i=0;i<cloud.n_sub;i++)
    {
        ii = i + mpi_vars.disps[mpi_vars.rank]/3;
        // compute cunbody error estimation
        res = (pow(ai_cpu[ii][0] - cloud.sub_aj[i][0],2)+
               pow(ai_cpu[ii][1] - cloud.sub_aj[i][1],2)+
               pow(ai_cpu[ii][2] - cloud.sub_aj[i][2],2));
        
        res = sqrt(res);
        
        res /= sqrt(ai_cpu[ii][0]*ai_cpu[ii][0]+
                    ai_cpu[ii][1]*ai_cpu[ii][1]+
                    ai_cpu[ii][2]*ai_cpu[ii][2]);
        
        if(res>err_max) err_max = res;
        if(res<err_min) err_min = res;
        err_mean += res;
        
    }
    err_mean/=cloud.n_sub;
    printf("rank %d : err %e\n",mpi_vars.rank,err_mean);
    delete [] ai_cpu;
}


