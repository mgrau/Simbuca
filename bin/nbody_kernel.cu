/*
NBODY algorithm from NVIDIA SDK. 
Changes:
1) bodyBodyInteraction function computes Coulomb interaction without the factor qq/(4 pi eps0)
2) gravitation function, solving a problem with multithreadBodies ON 
*/


#include <math.h>
#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>

// Macros to simplify shared memory addressing
#define SX(i) sharedPos[i+blockDim.x*threadIdx.y]
// This macro is only used when multithreadBodies is true (below)
#define SX_SUM(i,j) sharedPos[i+blockDim.x*j]

__constant__ float softeningSquared = 0.0000001f;

struct SharedMemory
{
    __device__ inline operator       float4*()
    {
        extern __shared__ int __smem[];
        return (float4*)__smem;
    }

    __device__ inline operator const float4*() const
    {
        extern __shared__ int __smem[];
        return (float4*)__smem;
    }
};


static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) 
{
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


#define HANDLE_NULL( a ) {if (a == NULL) { \
                            printf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
                            exit( EXIT_FAILURE );}}




int nbody;
// body value arrays (pos + mass, velocity)
//float4 *p_m;
//float3 *acc;


// body value arrays for device
float4 *p_m_Dev;
float3 *acc_Dev;

//CUDA timing
//cudaEvent_t start, stop;

//CUDA PROPS
int p_nb = 512;
int q = 1;
bool bMT =false;
dim3 grid, threads;
int sharedMemSize;


void copyDataToDevice_Async(float4* h_p)
{
	HANDLE_ERROR (cudaMemcpyAsync(p_m_Dev, h_p, nbody * sizeof (float4), cudaMemcpyHostToDevice));
}

void copyDataToDevice_Sync(float4* h_p)
{
	HANDLE_ERROR (cudaMemcpy(p_m_Dev, h_p, nbody * sizeof (float4), cudaMemcpyHostToDevice));
}



void copyDeviceToHost_Async(float3 * res)
{	
	cudaMemcpyAsync (res, acc_Dev, nbody * sizeof (float3), cudaMemcpyDeviceToHost);
}

void copyDeviceToHost_Sync(float3 * res)
{	
	cudaMemcpy (res, acc_Dev, nbody * sizeof (float3), cudaMemcpyDeviceToHost);		
}

void freeDevice(void)
{
	cudaFree (p_m_Dev);
	cudaFree (acc_Dev);
	
}

void initDataToDevice(void)
{
	HANDLE_ERROR (cudaMalloc ((void **) &p_m_Dev, nbody * sizeof (float4)));
	HANDLE_ERROR (cudaMalloc ((void **) &acc_Dev, nbody * sizeof (float3)));
}



extern "C" void initialize(int n_,int p_, int q_ ,int gridx_)
{		
 	nbody= n_;
	//p_m = (float4 *) malloc (nbody * sizeof (float4));
	//acc = (float3 *) malloc (nbody * sizeof (float3));
	
	
 	initDataToDevice();
	
	cudaDeviceProp props;
	
	if(p_*q_*gridx_==0)
	{
		cudaGetDeviceProperties(&props, 0);
	
		p_nb=min(nbody,p_nb);
		
		while ((nbody > 0) && p_nb > 1 && (nbody / p_nb < int((unsigned)props.multiProcessorCount)))
       		{
            		p_nb /= 2;
            		q *= 2;
        	}

       		grid.x = (int)(nbody + (p_nb-1))/p_nb;
	}
	else
	{
		p_nb = p_;
		q = q_;
		grid.x = gridx_;
	}
	
	threads.x = p_nb;
	threads.y = q;
	threads.z = 1;
	grid.y = 1;
	grid.z = 1;
	printf("######################################\n");
	printf("            NBODY ALGO USED         \n");
	printf(" tiles of %d x %d bodies              \n",p_nb,p_nb);
	printf(" %d thread(s) per body               \n",q);
	printf(" %d blocks used                       \n",grid.x);
	printf("######################################\n");
	
	if (grid.x > 0 && threads.y == 1)
        {
          	bMT = false;
	  	
        }
       	else if (grid.x > 0)
       	{
          	bMT = true;
	}
	
	sharedMemSize = p_nb * q  * 4 * sizeof(float);	
}


extern "C" float finalize(void)
{
	float elapsed=0;
	freeDevice();
	return elapsed;
}


__device__  float3 bodyBodyInteraction(float3  ai, float4 bi, float4 bj) 
{
   float3 r;

    // r_ij  [3 FLOPS]
    r.x = bj.x - bi.x;
    r.y = bj.y - bi.y;
    r.z = bj.z - bi.z;

    // distSqr = dot(r_ij, r_ij) + EPS^2  [6 FLOPS]
    float distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
    distSqr += softeningSquared;

    // invDistCube =1/distSqr^(3/2)  [4 FLOPS (2 mul, 1 sqrt, 1 inv)]
    float invDist = rsqrt(distSqr);
    float invDistCube =  (invDist * invDist * invDist);

    // s = m_j * invDistCube [1 FLOP]
    // float s = bj.w * invDistCube;

    // a_i =  a_i + s * r_ij [6 FLOPS]
    /*
    ai.x += r.x*s;// * s;
    ai.y += r.y*s;// * s;
    ai.z += r.z*s;// * s;
    */
    
    ai.x += r.x*invDistCube;// * s;
    ai.y += r.y*invDistCube;// * s;
    ai.z += r.z*invDistCube;// * s;
    
    return ai;
}


__device__ float3 gravitation(float4 iPos, float3  accel)
{
    float4 * sharedPos = SharedMemory();
    //extern __shared__ float4 sharedPos[];
    // The CUDA 1.1 compiler cannot determine that i is not going to 
    // overflow in the loop below.  Therefore if int is used on 64-bit linux 
    // or windows (or long instead of long long on win64), the compiler
    // generates suboptimal code.  Therefore we use long long on win64 and
    // long on everything else. (Workaround for Bug ID 347697)

    unsigned long j = 0;


    // Here we unroll the loop to reduce bookkeeping instruction overhead
    // 32x unrolling seems to provide best performance

    // Note that having an unsigned int loop counter and an unsigned
    // long index helps the compiler generate efficient code on 64-bit
    // OSes.  The compiler can't assume the 64-bit index won't overflow
    // so it incurs extra integer operations.  This is a standard issue
    // in porting 32-bit code to 64-bit OSes.

#pragma unroll 32
    for (unsigned int counter = 0; counter < blockDim.x; counter++ ) 
    {
        accel = bodyBodyInteraction(accel, iPos, SX(j++)); 
    }

    return accel;
}

// WRAP is used to force each block to start working on a different 
// chunk (and wrap around back to the beginning of the array) so that
// not all multiprocessors try to read the same memory locations at 
// once.
#define WRAP(x,m) (((x)<m)?(x):(x-m))  // Mod without divide, works on values from 0 up to 2m

template <bool multithreadBodies>
__global__ void computeBodyAccel(float4 * positions, float3 * DevA, int numBodies)
{
    float4* sharedPos = SharedMemory();
    //extern __shared__ float4 sharedPos[];
    float3 acc = {0.0f, 0.0f, 0.0f};

    int p_nb = blockDim.x;
    int q = blockDim.y;
    int n = numBodies;
    int numTiles = n / (p_nb * q);
   // float3 * globalA = (float3*) DevA;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    float4 myposition = positions[index];
 

    for (int tile = blockIdx.y; tile < numTiles + blockIdx.y; tile++) 
    {
        sharedPos[threadIdx.x+blockDim.x*threadIdx.y] = 
            multithreadBodies ? 
            positions[WRAP(blockIdx.x + q * tile + threadIdx.y, gridDim.x) * p_nb + threadIdx.x] :
        positions[WRAP(blockIdx.x + tile,                   gridDim.x) * p_nb + threadIdx.x];

        __syncthreads();

        // This is the "tile_calculation" function from the GPUG3 article.
        acc = gravitation(myposition, acc);

        __syncthreads();
    }
	
    // When the numBodies / thread block size is < # multiprocessors (16 on G80), the GPU is 
    // underutilized.  For example, with a 256 threads per block and 1024 bodies, there will only 
    // be 4 thread blocks, so the GPU will only be 25% utilized. To improve this, we use multiple 
    // threads per body.  We still can use blocks of 256 threads, but they are arranged in q rows 
    // of p threads each.  Each thread processes 1/q of the forces that affect each body, and then 
    // 1/q of the threads (those with threadIdx.y==0) add up the partial sums from the other 
    // threads for that body.  To enable this, use the "--p=" and "--q=" command line options to 
    // this example. e.g.: "nbody.exe --n=1024 --p=64 --q=4" will use 4 threads per body and 256 
    // threads per block. There will be n/p = 16 blocks, so a G80 GPU will be 100% utilized.

    // We use a bool template parameter to specify when the number of threads per body is greater 
    // than one, so that when it is not we don't have to execute the more complex code required!
    if (multithreadBodies)
    {
        SX_SUM(threadIdx.x, threadIdx.y).x = acc.x;
        SX_SUM(threadIdx.x, threadIdx.y).y = acc.y;
        SX_SUM(threadIdx.x, threadIdx.y).z = acc.z;
        __syncthreads();
        // Save the result in global memory for the integration step
       // if (threadIdx.y == 0) 
        {
	
	acc.z =0.0f;
	acc.x =0.0f;
	acc.y =0.0f;
	
            for (int i = 0; i < blockDim.y; i++) 
            {
                acc.x += SX_SUM(threadIdx.x,i).x;
                acc.y += SX_SUM(threadIdx.x,i).y;
                acc.z += SX_SUM(threadIdx.x,i).z;
            }
        }
	//__syncthreads();
    } 
    //float3 acc3 = {acc.x,acc.y,acc.z};
    // globalA[index] = acc3;
    DevA[index] =acc;
}

	
extern "C"  void simulation_Async(float4 * h_p,float3 * res)
	{
	copyDataToDevice_Async(h_p);
	
	if(bMT)
	computeBodyAccel<true><<< grid,threads,sharedMemSize,0 >>>(p_m_Dev,acc_Dev,nbody);
	else
	computeBodyAccel<false><<< grid,threads,sharedMemSize,0 >>>(p_m_Dev,acc_Dev,nbody);
	
	copyDeviceToHost_Async(res);
	return;
	}
	
extern "C"  void simulation_Sync(float4 * h_p,float3 * res)
	{
	copyDataToDevice_Sync(h_p);
	
	if(bMT)
	computeBodyAccel<true><<< grid,threads,sharedMemSize,0 >>>(p_m_Dev,acc_Dev,nbody);
	else
	computeBodyAccel<false><<< grid,threads,sharedMemSize,0 >>>(p_m_Dev,acc_Dev,nbody);
	
	copyDeviceToHost_Sync(res);
	return;
	}
