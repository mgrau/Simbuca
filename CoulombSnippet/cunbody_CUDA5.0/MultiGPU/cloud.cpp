#include "main.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
#include <iostream>
_cloud cloud;
extern _mpi_vars mpi_vars;

void Initialization(int N)
{
 
    cloud.elements = N;
    cloud.n_sub = N/mpi_vars.nb_procs;
    cloud.xi = new double [cloud.elements][3];
    cloud.sub_xi = new double [cloud.n_sub][3];
    cloud.xj= new double [mpi_vars.nb_node_per_gpu *   cloud.n_sub][3]; // position array (gpu sub cloud)
    cloud.mj= new double [cloud.elements]; // mass
    cloud.sub_aj= new double [cloud.n_sub][3]; // acc array (MPI sub cloud)
    cloud.ai= new double [cloud.elements][3]; // acc array (total cloud)
    cloud.ai_temp= new double [cloud.elements][3]; // part of acc array (total cloud)
    
    
    // INITIALIZATION ON NODE 0
    srand(0x19740526);
    if(mpi_vars.rank==0)
    {
        for(int i=0;i<cloud.elements;i++)
		{
			for(int j=0;j<3;j++)
			{
                cloud.xi[i][j]= rand()/1e8;
                
			}
           // std::cout << cloud.xi[i][0] << " " << cloud.xi[i][1] << " " << cloud.xi[i][2] << " " << std::endl;
		}
	}
    for(int i=0;i<cloud.elements;i++)
    {
        cloud.mj[i] =1.;
    }
    for(int i=0;i<mpi_vars.nb_procs;i++)
	{
        mpi_vars.counts[i] = 3*cloud.n_sub;
        mpi_vars.disps[i] = 3*cloud.n_sub*i;
	}
    for(int i=0;i<mpi_vars.nb_gpu;i++)
	{
        mpi_vars.counts_col[i] = 3*cloud.n_sub*mpi_vars.nb_node_per_gpu;
        mpi_vars.disps_col[i] = 3*cloud.n_sub*i*mpi_vars.nb_node_per_gpu;
	}
    
    MPI_Scatterv(&cloud.xi[0][0], mpi_vars.counts, mpi_vars.disps, MPI_DOUBLE,cloud.sub_xi[0],cloud.n_sub*3,MPI_DOUBLE,0,mpi_vars.comm2d);
}

