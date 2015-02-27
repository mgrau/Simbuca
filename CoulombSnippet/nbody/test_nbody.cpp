/******************************************

Test of the NBODY algorithm for SIMBUCA 3.0
Based on the NBODY algorithm of NVIDIA SDK

author: Pierre Dupré
email: pierre.dupre@csnsm.in2p3.fr
*****************************************/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector_types.h>
#include <vector>
#include <time.h>
#include "mtrand.h"
// CUDA includes
#include "cuda.h"
#include "cuda_runtime_api.h"

using namespace std;



#define EPSILON 0.0000001f
#define G 1.0f


int n = 1024;
int t = 1;



// random function
MTRand algrand;

//gpu functions
extern "C" 	{
			void initialize(int n_, int p_,int q_ , int gridx_);
			float finalize(void);	
			//float3 * simulation(double * px_, double * py_, double * pz_,double * coulombfactor_);
			void simulation_Async(float4* h_p,float3 * res);
			void simulation_Sync(float4* h_p,float3* res);
		}

// cpu
double *px, *py, *pz, *cf;
double *ax, *ay, *az;
float4 *h_p;
double rtemp;
double temp;
int n_verbose=100;
inline double dist2_cpu(double x1, double y1, double z1, double x2, double y2 , double z2);
bool compare_to_cpu;
bool synchro;
bool benchmark;
bool verbose_cpu = true;
bool verbose_gpu = true;
//gpu
float3 *acc_res;

int p_gpu = 0;
int q_gpu = 0;
int gridx_gpu = 0;;

//CUDA timing
cudaEvent_t start, stop;
double start_cpu;
double stop_cpu;



/******* COMMAND OPTION *****/
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}



int main(int argc, char *argv[])
{
	start_cpu = time(0);
	
	synchro = cmdOptionExists(argv,argv+argc,"sync");
	benchmark = cmdOptionExists(argv,argv+argc,"benchmark");
	compare_to_cpu = cmdOptionExists(argv,argv+argc,"cpu");
	
	if(argc>2)
	{
		n = atoi(argv[1]);
		t = atoi(argv[2]);
	}
	/*
	if(argc==4)
	{
		n = atoi(argv[1]);
		t = atoi(argv[2]);
		p_gpu = atoi(argv[3]);
	}
	if(argc==5)
	{
		n = atoi(argv[1]);
		t = atoi(argv[2]);
		p_gpu = atoi(argv[3]);
		q_gpu = atoi(argv[4]);
	}
	if(argc==6)
	{
		n = atoi(argv[1]);
		t = atoi(argv[2]);
		p_gpu = atoi(argv[3]);
		q_gpu = atoi(argv[4]);
		gridx_gpu = atoi(argv[5]);
	}		
	*/	
	printf("\n\tNBODY test for SIMBUCA V3.0\n");
	printf("\tNumber of particles : %d\n",n);
	printf("\tNumber of iteration : %d\n",t);
	if(synchro)
	{
		printf("\tSynchro mode enable\n");
	}
	else
	{
		printf("\tAsynchro mode enable\n");
	}	
	if(benchmark)
	printf("\tBenchmark mode enable\n");
	if(compare_to_cpu)
	printf("\tComparison with CPU calculation\n");
// cpu
	px = new double[n];
	py = new double[n];
	pz = new double[n];
	cf = new double[n];
	ax = new double[n];
	ay = new double[n];
	az = new double[n];	
	memset(ax,0,n*sizeof(double));
	memset(ay,0,n*sizeof(double));
	memset(az,0,n*sizeof(double));
	
// gpu
	
	if(synchro)
	{
		h_p = (float4*) malloc (n*sizeof(float4));
		acc_res  = (float3 *) malloc (n * sizeof (float3));
	}	
	else
	{
		cudaHostAlloc((void**)&h_p,n*sizeof(float4),cudaHostAllocDefault);
		cudaHostAlloc( (void**)&acc_res,n * sizeof (float3),cudaHostAllocDefault);
	}	
	

// init		
	double seed=time(0);
	algrand.seed(seed);
	
	printf("\n*** CUDA SIMULATION STARTS ***\n");
	float elapsed = 0.0f;
	double dcounter=0.; 
	unsigned long int counter=0;
	if(benchmark)
	{
		for(int i=0;i<n;i++)
		{
			px[i] = algrand()*2-1;
			py[i] = algrand()*2-1;
			pz[i] = algrand()*2-1;
			cf[i] = 1.0;
			h_p[i].x = px[i];
			h_p[i].y = py[i];
			h_p[i].z = pz[i];
			h_p[i].w = cf[i];
			//cout << px[i] << ", " << py[i] << ", " <<pz[i] << endl;
			//cout << ax[i] << ", " << ay[i] << ", " <<az[i] << endl;
		}
	
		initialize(n,p_gpu,q_gpu,gridx_gpu);	
		cudaEventCreate (&start);
		cudaEventCreate (&stop);
		cudaEventRecord (start, 0);	

		if(synchro)
		{
	
			for (int steps = 0; steps < t; steps++)
			{
				simulation_Sync(h_p,acc_res);
			}
		}
		else
		{
			for (int steps = 0; steps < t; steps++)
			{
				simulation_Async(h_p,acc_res);
			}
		}	
	
		cudaEventRecord(stop, 0);	
		
    		while( cudaEventQuery(stop) == cudaErrorNotReady )
    		{
        		dcounter += 1.;
    		}
		
  		
		
		
		finalize();
		cudaEventElapsedTime(&elapsed, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		
		printf("*** CUDA SIMULATION ENDS ***\n\n");
		//if (DEBUG) printBodies(20); // print results
		printf("CUDA time elapsed: %f ms\n", elapsed);
		double tmp =(double) t*n*n*1e-6;
		
		tmp*=20.;
		
		tmp/=elapsed;
		
		printf("PERF = %f GFLOPS/s\n",tmp);
		
		if(!synchro)
			printf("Number of double addition on the CPU during GPU calculation : %d\n\n\n",(int)dcounter);
		
		
	}
	
	
	double res = 0;
	double err_max = 0;
	double err_mean = 0;
	double err_min = 1e8;
	if(compare_to_cpu)
	{
		
		initialize(n,p_gpu,q_gpu,gridx_gpu);
			
		for(int k=0;k<t;k++)
		{
			for(int i=0;i<n;i++)
			{
				px[i] = algrand()*2-1;
				py[i] = algrand()*2-1;
				pz[i] = algrand()*2-1;
				cf[i] = 1.0;
				h_p[i].x = px[i];
				h_p[i].y = py[i];
				h_p[i].z = pz[i];
				h_p[i].w = cf[i];
				ax[i] = 0;
				ay[i] = 0;
				az[i] = 0;
				//cout << px[i] << ", " << py[i] << ", " <<pz[i] << endl;
				//cout << ax[i] << ", " << ay[i] << ", " <<az[i] << endl;
			}
			
			if(synchro)
			{
				simulation_Sync(h_p,acc_res);
			}
			else
			{
				simulation_Async(h_p,acc_res);
			}
		
	
			for(int j_=0;j_ <n;j_++)
			{
          			for(int i_=0;i_ <n;i_++)
	  			{
        				rtemp = EPSILON + dist2_cpu(px[i_],py[i_],pz[i_],px[j_],py[j_],pz[j_]);
                			temp = 1./ sqrtf(rtemp*rtemp*rtemp);
					ax[j_]-=(px[j_]-px[i_])*temp;
					ay[j_]-=(py[j_]-py[i_])*temp;
					az[j_]-=(pz[j_]-pz[i_])*temp;     
                		}
			}
			
			
			//err_max = 0;
			for(int i=0;i<n;i++)
			{
			// compute cunbody error estimation
				res = (pow(acc_res[i].x - ax[i],2) + pow( acc_res[i].y - ay[i],2) + pow(acc_res[i].z - az[i],2));
				res = sqrt(res);
				res /= sqrt(ax[i]*ax[i] + ay[i]*ay[i] + ay[i]*ay[i]);
				if(res>err_max) err_max = res;
				if(res<err_min) err_min = res;
				err_mean += res;
	
			}
			
		} // end iteration loop
		err_mean/=(n*t);
		printf("*** CUDA SIMULATION ENDS ***\n\n");
		finalize();
		printf("Maximal error  between CPU and GPU : %e\n",err_max);
		printf("Minimal error  between CPU and GPU : %e\n",err_min);
		printf("Mean error between CPU and GPU : %e\n",err_mean);
	}	// end if cpu 	
	

	free(px); free(py); free(pz); free(ax); free(ay); free(az);
	if(synchro)
	{
	free(h_p);
	free(acc_res);
	}
	else
	{
	cudaFreeHost(h_p);
	cudaFreeHost(acc_res);
	}
	
	return 0;
}

inline double dist2_cpu(double x1, double y1, double z1, double x2, double y2 , double z2)
{
return ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}
