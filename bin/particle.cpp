#include "particle.h"

Particle::Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz) {
    px = _x;
    py = _y;
    pz = _z;
    pvx = _vx; 
    pvy = _vy;
    pvz = _vz;
}

Particle::Particle() {
    px=py=pz=pvx=pvy=pvz=0.0;
}

ostream& operator<<(ostream& os,Particle& _p) {
    os<< _p.px<<"\t"<<_p.py<<"\t"<<_p.pz<<"\t"<<_p.pvx<<"\t"<<_p.pvy<<"\t"<<_p.pvz;
    return os;                 
}
