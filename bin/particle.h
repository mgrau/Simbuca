#ifndef PARTICLE_H
#define PARTICLE_H
#include <iostream>
using namespace std;

struct Particle
{   
      double px,py,pz,pvx,pvy,pvz;
      Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz); 
      Particle(); //default constructor he. Is made for Global Pool Particle Variables.
};

ostream& operator<<(ostream& os,Particle _p);

#endif

