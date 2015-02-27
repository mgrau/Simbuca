


#ifdef __MPI_ON__
#include "mpi_funcs.h"
#include <iostream>
#include "stdio.h"
#if defined(__CUNBODY_ON__) || defined(__NBODY_ON__)
#include "cuda.h"
#include "cuda_runtime_api.h"

#endif //GPU


using namespace std;

_mpi_vars::_mpi_vars()
{

    MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
   
    

    const int ndims = 2;
    const int reorder = 0;
    int periods[2] = { 0, 0 };
    int dims[2];

    coords = new int [2];
#if defined(__CUNBODY_ON__) || defined(__NBODY_ON__)    
    cudaGetDeviceCount(&(nb_gpu));
    if(nb_procs%nb_gpu!=0)
    {
        if (rank==0) {
            SLogger slogger("trapparameters");
            slogger << ERROR << "Number of nodes is not modulo the number of GPU" << SLogger::endmsg;
            exit(EXIT_FAILURE);
        }
 
    }
    ID_node_GPU = new int [nb_gpu];
    int remain_dims_h[2] = {0,1}; // communication only btw  row
    int remain_dims_v[2] = {1,0}; //communication only btw  column
#else
    nb_gpu =1 ;
#endif
    //nb_gpu =1 ; // short cut REMOVE AT THE END
    nb_node_per_gpu =   nb_procs /nb_gpu;
    
  
    // Create 2D cartesian topology
    dims[0]=nb_gpu;
    dims[1]=nb_node_per_gpu;
    MPI_Dims_create(nb_procs,ndims,dims);
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm2d);
    MPI_Comm_rank(comm2d,&(rank)); 
    MPI_Cart_coords(comm2d,rank,ndims,(coords));
    
    if (rank==0) {
        
        cout << "\n MPI execution with " <<  nb_procs <<" processes ("<<dims[0]<<" x "<<dims[1]<<")\n"<<endl;
       
    }
    
#if defined(__CUNBODY_ON__) || defined(__NBODY_ON__)     
    MPI_Cart_sub(comm2d,remain_dims_h,&comm1d_h);
    MPI_Cart_sub(comm2d,remain_dims_v,&comm_gpu);
    for(int dev=0;dev<nb_gpu;dev++)
        ID_node_GPU[dev] = dev*nb_node_per_gpu;

    
    
    counts_F4 = new int[nb_procs];
    disps_F4  = new int[nb_procs];
    counts_F3 = new int[nb_procs];
    disps_F3  = new int[nb_procs];
    counts_col = new int[nb_node_per_gpu];
    disps_col  = new int[nb_node_per_gpu];
    counts_row = new int[nb_node_per_gpu];
    disps_row  = new int[nb_node_per_gpu];
    disps_row_return = new int[nb_node_per_gpu];

    counts_Arow = new int[nb_gpu];
    disps_Arow  = new int[nb_gpu];
    
    // Create Datatype MPI_FlOAT3 for broadcasting calculated acc from GPU
	MPI_Type_contiguous(3, MPI_FLOAT, &MPI_FLOAT3);
    MPI_Type_commit(&MPI_FLOAT3);
	// Create Datatype MPI_FlOAT4 for gathering positions for GPU
	MPI_Type_contiguous(4, MPI_FLOAT, &MPI_FLOAT4);
    MPI_Type_commit(&MPI_FLOAT4);
    
    
#endif    
    // counts and disps
    counts = new int[nb_procs];
    disps  = new int[nb_procs];
    counts_3 = new int[nb_procs];
    disps_3  = new int[nb_procs];

    
    MPI_Barrier(comm2d);
}

_mpi_vars::~_mpi_vars()
{
    
    MPI_Comm_free(&comm2d);
    
    delete [] counts;
    delete [] disps;
    delete [] counts_3;
    delete [] disps_3;
#ifdef __CPUNBODY_ON__
    //delete [] pos_tot;
#endif
    
#if defined(__CUNBODY_ON__) || defined(__NBODY_ON__)   
    //MPI_Comm_free(&comm_gpu);
    MPI_Comm_free(&comm1d_h);
    
    delete [] coords;
    delete [] ID_node_GPU;
    delete [] disps_F4  ;
    delete []counts_F4 ;
    delete [] disps_F3  ;
    delete []counts_F3 ;
    delete [] counts_row ;
    delete [] disps_row  ;
    delete [] disps_row_return  ;
    delete [] counts_col ;
    delete [] disps_col  ;
    delete [] counts_Arow;
    delete [] disps_Arow ;
    //MPI_Group_free(&server_grp);
    //MPI_Group_free(&gpu_grp);
#endif // GPU
}


void _mpi_vars::ComputeCounts_Disps(int nrparticles)
{
    MPI_Allgather(&nrparticles,1,MPI_INT,counts,1,MPI_INT,comm2d);
    disps[0] = 0;
    disps_3[0] = 0;
    n_part_tot = 0;
    for(int i=0;i<nb_procs;i++)
    {
        counts_3[i] = counts[i]*3;
        n_part_tot += counts[i];
    }
    for(int i=1;i<nb_procs;i++)
    {
        disps[i] = counts[i-1] + disps[i-1];
        disps_3[i] = counts_3[i-1] + disps_3[i-1];
    }
   
#if defined(__CUNBODY_ON__) || defined(__NBODY_ON__)      
    int ii;
    if(nb_gpu>1)
    {
        //ROW
        for(int i_=0;i_<nb_node_per_gpu;i_++)
        {
            ii = i_ + coords[0]*nb_node_per_gpu;
            counts_row[i_] = counts_3[ii];
            if(i_==0)
            {
                disps_row[0] = 0 ;
                disps_row_return[0] = disps_3[ID_node_GPU[coords[0]]];
            }
            else
            {
                disps_row[i_] =  disps_row[i_-1]+counts_row[i_-1];
                disps_row_return[i_] = disps_row_return[i_-1]+counts_row[i_-1];
            }
        }
        //AROW
        for(int j_=0;j_<nb_gpu;j_++)
        {
            counts_Arow[j_] =0;
            
            for(int i_=0;i_<nb_node_per_gpu;i_++)
            {
                ii = i_ + coords[0]*nb_node_per_gpu;
                counts_Arow[j_]+= counts[ii];
            }
           
            if(j_==0)
            {
                disps_Arow[j_]=0;
  
            }
            else
            {

                disps_Arow[j_] = disps_Arow[j_-1] + counts_Arow[j_-1];
            }
 
            //if(rank==0)
            //   printf("%d %d %d \n",j_,counts_Arow[j_],disps_Arow[j_]);

        }
    }
    //printf("%d %d %d %d\n",rank,counts[rank], disps[rank],n_part_tot);
    //printf("%d %d %d %d\n",rank,counts_3[rank], disps_3[rank],n_part_tot);
    charge = new double [n_part_tot];
#endif
}
#ifdef __CPUNBODY_ON__
void _mpi_vars::Init_GlobalPos(int nrparticles)
{
    if(n_part_tot==0)
        ComputeCounts_Disps(nrparticles);
        
    if(n_part_tot!=0)
        pos_tot = new double[n_part_tot][3];
    else
    {
        if(rank==0)
            cout << "n_part_tot = 0 in MPI VARS" << endl;
        
    }
}
#else
void _mpi_vars::LoadCharge( vector<int > ion_charge)
{
    MPI_Allgatherv(&ion_charge[0],ion_charge.size(),MPI_INT,charge,counts,disps,MPI_DOUBLE,comm2d);
}
#endif // __CPUNBODY_ON__




#endif // __MPI_ON__



