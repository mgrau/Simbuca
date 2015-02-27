#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPSILON 0.0000001f
#define DEBUG false
inline float dist2(float x1, float y1, float z1, float x2, float y2 , float z2);
int n; // number of part

float *px, *py, *pz;
float *ax, *ay, *az;

float rtemp;
float temp;
using namespace std;
int main(int argc, char *argv[])
{
	if(argc >1)
	{
		n = atoi(argv[1]);
	}
	else
	{
		n = 10;
	}
	
	px = new float[n];
	py = new float[n];
	pz = new float[n];
	ax = new float[n];
	ay = new float[n];
	az = new float[n];	
	memset(ax,0,n*sizeof(float));
	memset(ay,0,n*sizeof(float));
	memset(az,0,n*sizeof(float));
	
	//srand(0);
	
	for(int i=0;i<n;i++)
	{
		px[i] = rand()*1e-10;
		py[i] = rand()*1e-10;
		pz[i] = rand()*1e-10;
		//cout << px[i] << ", " << py[i] << ", " <<pz[i] << endl;
		//cout << ax[i] << ", " << ay[i] << ", " <<az[i] << endl;
	}	
	
	

 	for(int j_=0;j_ <n;j_++)
	{
          	for(int i_=0;i_ <n;i_++)
	  	{
        	rtemp = EPSILON + dist2(px[i_],py[i_],pz[i_],px[j_],py[j_],pz[j_]);
                temp = 1/ sqrtf(rtemp*rtemp*rtemp);
		ax[j_]+=(px[j_]-px[i_])*temp;
		ay[j_]+=(py[j_]-py[i_])*temp;
		az[j_]+=(pz[j_]-pz[i_])*temp;     
                }
	}
	
	for(int i=0;i<n;i++)
	{
		//cout << px[i] << ", " << py[i] << ", " <<pz[i] << endl;
		cout << ax[i] << ", " << ay[i] << ", " <<az[i] << endl;
	}	
		
return 1;	
}


inline float dist2(float x1, float y1, float z1, float x2, float y2 , float z2)
{
return ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}
