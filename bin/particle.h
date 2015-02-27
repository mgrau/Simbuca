//Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H
#include <iostream>
#include <vector>
using namespace std;

//about Particle
class Particle
{   
public:
      //global vars 
      double px,py,pz,pvx,pvy,pvz;
   // vector<  double  > pos;
   // vector<  double  > vel;
      //accessable functions
      Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz); 
      Particle(); //default constructor he. Is made for Global Pool Particle Variables.

};

ostream& operator<<(ostream& os,Particle _p);


#endif

