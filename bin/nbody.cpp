#ifdef __NBODY_ON__
#include "nbody.h"


#define is_power_of_two(x) ( ((x&(x-1)) == 0) && (x!=0) )



using namespace std;

//gpu functions
extern "C" 	{
    void initialize(int n_,bool withCharge);
    float finalize(void);
    void simulation_Async(float4* h_p,float3 * res,int nbody_,int p_, int q_, int ntiles_);
    void simulation_Async_begin(float4* h_p,float3 * res,int nbody_,int p_, int q_, int ntiles_);
    void simulation_Async_end(float4* h_p,float3 * res,int nbody_,int p_, int q_, int ntiles_);
}

_nbody::_nbody(){}
_nbody::~_nbody(){}


void _nbody::Initialization(const IonCloud & _cloud, bool _withCharge) {
#ifdef __MPI_ON__
    n_cpu = mpiv->n_part_tot;
    if(mpiv->coords[1]==0) // node using GPU
    {
#else
        n_cpu = _cloud.nrparticles;
#endif
        cout << "Begin NBODY init" << endl;
        int devID;
        cudaDeviceProp props;
        cudaGetDevice(&devID);
        cout << "cudaGetDeviceProperties begin" << endl;
        cudaGetDeviceProperties(&props, devID);
        cout << "cudaGetDeviceProperties end" << endl;
        name = props.name;
        sm = props.multiProcessorCount;
        step = sm*16;
        cc = props.maxThreadsPerBlock;
        cout << "devID: " << devID << endl;
        cout << "Name: " << name << endl;
        cout << "Processors: " << sm << endl;
        n_max_alloc = FindSupN(n_cpu);

        cout << "Nmax: " << n_max_alloc << endl;
        LoadFile();
        cout << "Loaded File..." << endl;
        cout << "Will cudaHostAlloc " <<n_max_alloc * sizeof (float4) << " plus " << n_max_alloc * sizeof (float3) << "\n";

        cudaHostAlloc( (void**)&host_pos,n_max_alloc*sizeof(float4),cudaHostAllocDefault);
        cudaHostAlloc( (void**)&host_acc,n_max_alloc*sizeof(float3),cudaHostAllocDefault);

        cout << "Memory allocated..." << endl;

        //host_pos  = (float4 *) malloc (n_max_alloc * sizeof (float4));
        //host_acc  = (float3 *) malloc (n_max_alloc * sizeof (float3));
        initialize(n_max_alloc,_withCharge);
#ifdef __MPI_ON__
    }
    MPI_Bcast(&n_max_alloc,1, MPI_INT,0,mpiv->comm1d_h);
    if(mpiv->coords[1]!=0) // node using GPU
    {
#endif
        cout << "Will malloc " <<n_max_alloc * sizeof (float4) << " plus " << n_max_alloc * sizeof (float3) << "\n";
        host_pos  = (float4 *) malloc (n_max_alloc * sizeof (float4));
        host_acc  = (float3 *) malloc (n_max_alloc * sizeof (float3));
#ifdef __MPI_ON__

    }
    host_pos_sub = (float4 *) malloc (_cloud.nrparticles * sizeof (float4));
    host_acc_sub = (float3 *) malloc (_cloud.nrparticles * sizeof (float3));
#endif
}

int _nbody::FindSupN(int n_) {
    int n_temp = step;
    if(n_ < step)
    {
        cout << " n too small" << endl;
        exit(EXIT_FAILURE);
    }
    while(n_>n_temp)
    {
        n_temp +=step;

    }
    return n_temp;
}

void _nbody::Finalize() {
#ifdef __MPI_ON__
    if(mpiv->coords[1]==0) // node using GPU
    {
#endif
        finalize();
        cudaFreeHost(host_pos);
        cudaFreeHost(host_acc);
#ifdef  __MPI_ON__
    }
    else
    {
        free(host_pos);
        free(host_acc);
    }
    free(host_pos_sub);
    free(host_acc_sub);
#endif
}

void _nbody::LoadPosition(const IonCloud &_cloud) {
#ifndef  __MPI_ON__

    for(int i=0;i<n_cpu;i++)
    {
        host_pos[i].x = _cloud.pos2[i][0];
        host_pos[i].y = _cloud.pos2[i][1];
        host_pos[i].z = _cloud.pos2[i][2];
        host_pos[i].w = 1;// CHARGE
    }

    // ghost particles
    if(n_cpu<n_gpu)
        for(int i=n_cpu;i<n_gpu;i++)
        {
            host_pos[i].x = 1e8;
            host_pos[i].y = 0;
            host_pos[i].z = 0;
            host_pos[i].w = 0;
        }
#else
    MPI_Barrier(mpiv->comm2d);
    for(int i=0;i<_cloud.nrparticles;i++)
    {
        host_pos_sub[i].x = _cloud.pos2[i][0];
        host_pos_sub[i].y = _cloud.pos2[i][1];
        host_pos_sub[i].z = _cloud.pos2[i][2];
        host_pos_sub[i].w = 1; // CHARGE

    }

    MPI_Gatherv(&host_pos_sub[0],_cloud.nrparticles,mpiv->MPI_FLOAT4,&host_pos[0],mpiv->counts,mpiv->disps,mpiv->MPI_FLOAT4,0,mpiv->comm2d); // here

    if(mpiv->coords[1]==0) // node using GPU
    {
        if(n_cpu<n_gpu)
            for(int i=n_cpu;i<n_gpu;i++)
            {
                host_pos[i].x = 1e8;
                host_pos[i].y = 0;
                host_pos[i].z = 0;
                host_pos[i].w = 0;
            }
    }

#endif
}

void _nbody::SetParameterWithoutFile() {
    int n,p,q,t;
    n_gpu =  FindSupN(n_cpu);

    n = n_gpu;
    p = min(cc,n);
    q = 1;
    t = sm;
    while( n/t>cc)
    {
        t*=2;
    }
    p = n/ t;
    q=1;
    while((n%(q*p))==0&&(q*p<=cc)&&(q<p))
    {
        q*=2;
    }
    q/=2;
    if(!is_power_of_two(n)&&n<6400)
        if(p%64==0)
        {
            p/=2;
            q*=2;
            t*=2;
        }
    p_gpu =p ;
    q_gpu =q ;
    t_gpu =t ;
}
void _nbody::SetParameterWithFile() {
    int gpu_i = 0;

    while(n_cpu>Vgpu_n[gpu_i])
    {
        gpu_i++;
    }

    n_gpu = Vgpu_n[gpu_i];
    p_gpu = Vgpu_p[gpu_i];
    q_gpu = Vgpu_q[gpu_i];
    t_gpu = Vgpu_ntiles[gpu_i];
}

void _nbody::SetParameters() {
    if(!file||(n_cpu>n_max_file))
    {
        SetParameterWithoutFile();

    }
    else
    {
        SetParameterWithFile();
    }
}
void _nbody::force_begin(const IonCloud &_cloud) {
#ifdef  __MPI_ON__
    n_cpu = mpiv->n_part_tot;
#else
    n_cpu = _cloud.nrparticles;
#endif

#ifdef  __MPI_ON__
    if(mpiv->coords[1]==0)
        SetParameters();
#else
    SetParameters();
#endif

    LoadPosition(_cloud);

#ifdef  __MPI_ON__

    if(mpiv->coords[1]==0)
    {
        cudaEventCreate (&end_gpu2);
        simulation_Async(host_pos,host_acc,n_gpu, p_gpu, q_gpu, t_gpu);
        cudaEventRecord(end_gpu2,0);

    }

#else
    cudaEventCreate (&end_gpu2);
    simulation_Async(host_pos,host_acc,n_gpu, p_gpu, q_gpu, t_gpu);
    cudaEventRecord(end_gpu2,0);
#endif
}


void _nbody::force_end(const IonCloud &_cloud) {
#ifdef  __MPI_ON__
    int iii=0;
    if(mpiv->rank==0)
    {
        while( cudaEventQuery(end_gpu2) == cudaErrorNotReady )
        {

        }
        cudaEventDestroy(end_gpu2);

    }
    MPI_Scatterv(host_acc,mpiv->counts,mpiv->disps,mpiv->MPI_FLOAT3,host_acc_sub,_cloud.nrparticles,mpiv->MPI_FLOAT3,0,mpiv->comm2d);

#else
    while( cudaEventQuery(end_gpu2) == cudaErrorNotReady )
    {

    }
    cudaEventDestroy(end_gpu2);
#endif
}



void _nbody::LoadFile() {
    int tmp_n,tmp_p,tmp_q,tmp_ntiles;
    double tmp1;


    unsigned found = name.find(' ');
    while(found>name.size())
    {
        name.replace(found,1,"_");
        found = name.find(' ');
    }
    stringstream filename;
    //filename << "../../libraries/NBODY/" << name << "_NB0DY.dat";
    filename << "" << name << "_NB0DY.dat";


    GPUfile.open(filename.str().c_str(),ios::in);

    if(!GPUfile)
    {
        file = false;
        return;
    }
    else
    {
        file = true;
        while(!GPUfile.eof())
        {
            GPUfile >> tmp_n >> tmp_p >> tmp_q  >> tmp_ntiles >> tmp1 ;//>> tmp2;
            //  cout << tmp_n << " " << tmp_p << " " << tmp_q << " " << tmp_ntiles << " " << tmp1 << " " << endl;
            Vgpu_n.push_back(tmp_n);
            Vgpu_p.push_back(tmp_p);
            Vgpu_q.push_back(tmp_q);
            Vgpu_ntiles.push_back(tmp_ntiles);
        }
        ;
        Vgpu_n.pop_back();
        Vgpu_p.pop_back();
        Vgpu_q.pop_back();
        Vgpu_ntiles.pop_back();

        n_max_file = tmp_n;
        n_max_alloc = max(n_max_alloc,n_max_file);
        GPUfile.close();
        return;
    }
}

#ifdef __MPI_ON__
void _nbody::InitMPI(_mpi_vars * _mpiv)
{
    mpiv = _mpiv;
}
#endif //MPI

#endif // __NBODY_ON__

