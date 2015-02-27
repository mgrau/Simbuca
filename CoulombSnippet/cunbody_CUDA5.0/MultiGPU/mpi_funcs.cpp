#include "main.h"
#include "cuda.h"
#include "cuda_runtime_api.h"

_mpi_vars mpi_vars;
void init_mpi(int argc, char *argv[])
{
    const int ndims = 2;
    const int reorder = 0;
    int periods[2] = { 0, 0 };
    
    int ierr;
    mpi_vars.coords = new int [2];
    int dims[2];

    
    MPI_Comm *comm2d = &(mpi_vars.comm2d);
    MPI_Comm *comm1d_h = &(mpi_vars.comm1d_h);
    MPI_Comm *comm1d_v = &(mpi_vars.comm_gpu);
    int remain_dims_h[2] = {0,1}; // communication only btw  row
    int remain_dims_v[2] = {1,0}; //communication only btw  column
    ierr = MPI_Init(&argc,&argv);
    
    MPI_Comm_size(MPI_COMM_WORLD,&(mpi_vars.nb_procs));
    
    
    cudaGetDeviceCount(&(mpi_vars.nb_gpu));
    if(mpi_vars.nb_procs%mpi_vars.nb_gpu!=0)
    {
        if (mpi_vars.rank==0) {
            printf("\n");
            printf("ERROR: Number of nodes is not modulo the number of GPU");
            printf("\n");
            exit(EXIT_FAILURE);
        }
 
    }
    mpi_vars.nb_node_per_gpu =   mpi_vars.nb_procs /mpi_vars.nb_gpu;
    mpi_vars.ID_node_GPU = new int [mpi_vars.nb_gpu];
    
    
    // Create 2D cartesian topology
    dims[0]=mpi_vars.nb_gpu;
    dims[1]=mpi_vars.nb_node_per_gpu;
    MPI_Dims_create(mpi_vars.nb_procs,ndims,dims);
    
       
    
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,comm2d);
    MPI_Comm_rank(*comm2d,&(mpi_vars.rank));
   
    
    MPI_Cart_coords(*comm2d,mpi_vars.rank,ndims,(mpi_vars.coords));
    
    if (mpi_vars.rank==0) {
        printf("\n");
        printf(" MPI execution with  %d processes (%d x %d)\n",mpi_vars.nb_procs,dims[0],dims[1]);        
        printf("\n");
    }
    MPI_Cart_sub(*comm2d,remain_dims_h,comm1d_h);
    MPI_Cart_sub(*comm2d,remain_dims_v,comm1d_v);
    for(int dev=0;dev<mpi_vars.nb_gpu;dev++)
        mpi_vars.ID_node_GPU[dev] = dev*mpi_vars.nb_node_per_gpu;
    
    
    
    mpi_vars.counts = new int[mpi_vars.nb_procs];
    mpi_vars.disps  = new int[mpi_vars.nb_procs];
    mpi_vars.counts_col = new int[mpi_vars.nb_node_per_gpu];
    mpi_vars.disps_col  = new int[mpi_vars.nb_node_per_gpu];
    mpi_vars.counts_row = new int[mpi_vars.nb_node_per_gpu];
    mpi_vars.disps_row  = new int[mpi_vars.nb_node_per_gpu];
    mpi_vars.disps_row_return = new int[mpi_vars.nb_node_per_gpu];    
}

void
clean_mpi()
{
    int ierr;
    //MPI_Comm_free(&mpi_vars.comm_gpu);
    MPI_Comm_free(&mpi_vars.comm1d_h);
    MPI_Comm_free(&mpi_vars.comm2d);
    delete [] mpi_vars.coords;
    delete [] mpi_vars.ID_node_GPU;
    delete [] mpi_vars.disps  ;
    delete []mpi_vars.counts ;
    delete [] mpi_vars.counts_row ;
    delete [] mpi_vars.disps_row  ;
    delete [] mpi_vars.disps_row_return  ;
    delete [] mpi_vars.counts_col ;
    delete [] mpi_vars.disps_col  ;
    //MPI_Group_free(&mpi_vars.server_grp);
    //MPI_Group_free(&mpi_vars.gpu_grp);
    ierr = MPI_Finalize();
}
