/* sample.cpp
 * by Tsuyoshi Hamada
 * adapted by Simon Van Gorp 
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
using namespace std;

#define NMAX (25)
double xj[NMAX][3];
double mj[NMAX];
double xi[NMAX][3];
double ai[NMAX][3];

/* 
 * cunbody1_force is the only subroutine of libcunbody1.a.
 * cunbody1_force calculates gravitational interactions between
 * i-th and j-th particles.
 */
extern "C" {
   void cunbody1_force(double xj[][3], // position of j-th particles
                    double mj[],    // mass of j-th particles
                    double xi[][3], // position of i-th particles
                    double eps2,    // softening parameter
                    double ai[][3], // force of i-th particles
                    int ni,         // number of i-th particles
                    int nj);        // number of j-th particles
}

/*  main routine */

int main()
{
  struct timeval start, end;
  long mtime, seconds, useconds;    
  gettimeofday(&start, NULL);

  int i, dim, ni, nj;
  double eps2 = 1e-20;   
  
  double mass=1.0;  //35*1.660e-27; //mass of p.e. 35Argon in kg

  srand(0x19740526);

  double sim_start=time(0);
  time_t rawtime2;time ( &rawtime2 );
  /* SETUP I,J-PARTICLES */
  ni = nj = NMAX;

  for(i=0; i<ni; i++){
    for(dim=0; dim<3; dim++)
      xi[i][dim] = xj[i][dim] = rand();
      ai[i][dim]=0.0;
    mj[i] = mass;
  }
  //example 1
//   xi[0][0]=xj[0][0]=0.01;
//   xi[1][0]=xj[1][0]=0.02;
  //both particles on x-axis. one at 1cm other at 2cm
  // -> 2 forces in the same direction.  
  
 

   //example 3
   //particle 1
   xi[0][0]=xj[0][0]=-4.0;
   xi[0][1]=xj[0][1]=0.0;
   
   //particle 2
   xi[1][0]=xj[1][0]=4.0;
   xi[1][1]=xj[1][1]=0.0;
   
   //particle 3
   xi[2][0]=xj[2][0]=0.0;
   xi[2][1]=xj[2][1]=3.0;
   mj[0]=1.0;mj[1]=1.0;mj[0]=1.0;


  //example from Physics for scientist and Engineers Volume II, Fishbane, Gasiorowicz and Thornton pg 621
  //the magnitude of force should be (without ke)

  //in both examples F1 should be negative and F2 positive.
  
  /* CALL GeForce8800GTX */
     /*      for(int i=0; i< NMAX;i++){
                 for(int j=0; j<3;j++){
                 xj[i][j]=0.0;xi[i][j]=0.0;ai[i][j]=0.0;mj[i]=1.0;eps2=1e-10;
                 }
         }*/
cout<<"nrparticles= "<<NMAX<<endl;  
  cunbody1_force(xj, mj, xi, eps2, ai, NMAX, NMAX); 
  
 /*CPU*/
/*	 double r_ij;
 	for(int j_=0;j_ <NMAX;j_++){
          for(int i_=0;i_ <NMAX;i_++){
                 if(i_ != j_){
                      r_ij=sqrt((xi[i_][0]-xi[j_][0])*(xi[i_][0]-xi[j_][0])
		      	       +(xi[i_][1]-xi[j_][1])*(xi[i_][1]-xi[j_][1])
			       +(xi[i_][2]-xi[j_][2])*(xi[i_][2]-xi[j_][2]));
                      r_ij=fabs(1.0/(r_ij*r_ij*r_ij));
                      //i.e. Force particle i exerts(uitoefend) on particle j   
                      ai[j_][0]+=(xi[j_][0]-xi[i_][0])*(r_ij);
                      ai[j_][1]+=(xi[j_][1]-xi[i_][1])*(r_ij);
                      ai[j_][2]+=(xi[j_][2]-xi[i_][2])*(r_ij);  
                 }
          }
      }
  
  */
  /* OUTPUT RESULTS */

for(i=0; i<2; i++)
   for(dim=0; dim<3; dim++)
     printf("a[%d][%d] = %e\n",i, dim, -ai[i][dim]);


 //for(i=0; i<2; i++)
//  i = 0;
//    for(dim=0; dim<3; dim++)
//      printf("FCoulomb[%d][%d] = %e\n",i, dim, -ai[i][dim]*(8.987551e9)*mj[i]);
 
/*
   Gravity -> Coulomb:
           multiplication factor for accelerations. Supossing that all masses and charges are the same.
                          =q*q*ke/(G*m*m) 
 */



    gettimeofday(&end, NULL);
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    printf("Simulation time: %ld milliseconds\n", mtime);


  return (0);
}
