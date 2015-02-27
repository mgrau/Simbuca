/* sample.cpp
 * by Tsuyoshi Hamada
 */
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "builtin_types.h"
#include "octgrav.h"
#include <cstring>
#include <vector>
using std::vector;

#include <iostream>
using namespace std;

/* 
The tree code by Gaburov et all. from Aug. 2010.
float4: .x ,.y, .z the positions, .w can be used as mass or charge.
The function "tree.evaluate_gravity(pos, acc);"
calculates the interaction between the particles via this equation:

acc.x = [pos.w / (r*r)]_x
acc.y = [pos.w / (r*r)]_y
acc.z = [pos.w / (r*r)]_z

OPEN QUESTIONS:
-why does it only work from >= 9 particles ?
-why is acc a float4 and not float3, so what does acc.w mean ?

 */

/*  main routine */

const int nrparticles = 25000;

int main()
{ 

     struct timeval start, end;
     long mtime, seconds, useconds;    
     gettimeofday(&start, NULL);


    const double amu = 1.6605402E-27;//massa van een ion in kg
    const double ke=8.987551787E9; //=1/(4*pi*eps0) [N*(m/C)^2]
    const double el_charge  = 1.60217733e-19;	//el lading [C]
    
      srand(0x19740526);


    std::vector< float4 > pos(nrparticles);
    std::vector< float4 > acc(nrparticles);
    std::vector< float3 > vel(nrparticles);
    std::vector< float3 > vel1(nrparticles);
    //supossing each particles has the same mass and charge
    //otherwise one has to input the charge or mass in the .w value of the pos array.
    std::vector<float> mass(nrparticles);



    static octgrav tree;
    float eps=1e-20;
    float theta = 0.5;
    
  //initialize the tree
  tree.set_softening(eps);
  tree.set_opening_angle(theta);
  
    
  int i;
  double eps2 = 1e-20;

  srand(0x19740526);
  double sim_start=time(0);
  time_t rawtime2;time ( &rawtime2 );
  /* SETUP I,J-PARTICLES */

  //example 1
//   xi[0][0]=xj[0][0]=0.01;
//   xi[1][0]=xj[1][0]=0.02;
  //both particles on x-axis. one at 1cm other at 2cm
  // -> 2 forces in the same direction.  
  
  for(int i=0; i<nrparticles;i++){ 
      mass[i]=1.0;
      pos[i].x=rand();
      pos[i].y=rand();
      pos[i].z=rand();
      vel[i].x=0.0;
      vel[i].y=0.0;
      vel[i].z=0.0;
      acc[i].x=0.0;
      acc[i].y=0.0;
      acc[i].z=0.0;   
               
  }
  
  
  //example 2
   //particle 1
   pos[0].x=-3.0;
   pos[0].y=0.0;
   pos[0].z=0;
   
   pos[0].w=pos[1].w=pos[2].w=1.0;
   acc[0].w=acc[1].w=acc[2].w=1.0;
   //particle 2
   pos[1].x=3.0;
   pos[1].y=0.0;
   pos[1].z=0;
   //particle 3
   pos[2].x=0.0;
   pos[2].y=4.0;
   pos[2].z=0;
   
    

  //example from Physics for scientist and Engineers Volume II, Fishbane, Gasiorowicz and Thornton pg 621
  //the magnitude of force should be (without ke)
  //
  //


  //in both examples F1 should be negative and F2 positive.
  
  // creating data
 /* for(i=0;i<nrparticles;i++)
    {
      pos[i].w = 35*amu;
    }
  for(i=0;i<nrparticles;i++)
    {
      pos[i].x=-0.001+algrand()*0.002;
      pos[i].y=-0.001+algrand()*0.002;
      pos[i].z=-0.001+algrand()*0.002;
    }
 */ 
  

  
  
  
  /* THE FORCE!*/
  tree.evaluate_gravity(pos, acc);

 
  
  /* OUTPUT RESULTS */
 
   //print out the output of the program, (calculates just F=1/(r*r)
 /*  for(i=0; i<4; i++){
	   cout<<endl;
      printf("Fi[%d].x = %e\n",i, acc[i].x);
      printf("Fi[%d].y = %e\n",i, acc[i].y);
      printf("Fi[%d].z = %e\n",i, acc[i].z);
      printf("Fi[%d].w = %e\n",i, acc[i].w);
    }
  */
  
  //print out Coulomb Force
/*  for(i=0; i<2; i++){
	  cout<<endl;
      printf("F_Coulomb[%d].x = %e\n",i, -acc[i].x*mass[i]*mass[i]*el_charge*el_charge*ke);
      printf("F_Coulomb[%d].y = %e\n",i, -acc[i].y*mass[i]*mass[i]*el_charge*el_charge*ke);
      printf("F_Coulomb[%d].z = %e\n",i, -acc[i].z*mass[i]*mass[i]*el_charge*el_charge*ke);
      printf("F_Coulomb[%d].w = %e\n",i, -acc[i].w);
  }
 */ 
  
    /* OUTPUT RESULTS */

/*
 for(i=0; i<2; i++){
      printf("acc[%d].x = %e\n",i,-acc[i].x);
      printf("acc[%d].y = %e\n",i,-acc[i].y);
      printf("acc[%d].z = %e\n",i,-acc[i].z);
      }
*/
cout<<"nrparticles= "<<nrparticles<<endl;  

 for(i=0; i<2; i++){
      printf("FCoulomb[%d].x = %e\n",i,-acc[i].x*(8.987551e9)*pos[i].w);
//      printf("FCoulomb[%d].y = %e\n",i,-acc[i].y*(8.987551e9)*pos[i].w);
//      printf("FCoulomb[%d].z = %e\n",i,-acc[i].z*(8.987551e9)*pos[i].w);
      }
      
    gettimeofday(&end, NULL);
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    printf("Simulation time: %ld milliseconds\n", mtime);     
      
      
      
 /*
   Gravity -> Coulomb:
           multiplication factor for accelerations. Supossing that all masses and charges are the same.
                          =q*q*ke/(G*m*m) 
 */
 
  return (0);
}
