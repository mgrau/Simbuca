//Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H
#include <iostream>
using namespace std;

class Particle
{   
    public:
        double px,py,pz,pvx,pvy,pvz;
        Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz); 
        Particle();
};

ostream& operator<<(ostream& os,Particle _p);

#endif

