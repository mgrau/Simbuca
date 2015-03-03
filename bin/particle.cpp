//Particle.cpp
#include "particle.h"


Particle::Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz){
	px = _x;
	py = _y;
	pz = _z;
	pvx = _vx; 
	pvy = _vy;
	pvz = _vz;
//    pos.resize(3);
//    vel.resize(3);
//    pos[0]=_x;
//    pos[1]=_y;
//    pos[2]=_z;
//    vel[0]=_vx;
//    vel[1]=_vy;
//    vel[2]=_vz;
    
}

Particle::Particle(){
	px=py=pz=pvx=pvy=pvz=0.0;
   // pos.resize(3);
   // vel.resize(3);
   // pos[0]=pos[1]=pos[2]=0;
   // vel[0]=vel[1]=vel[2]=0;
}


ostream& operator<<(ostream& os,Particle& _p){
     os<< _p.px<<"\t"<<_p.py<<"\t"<<_p.pz<<"\t"<<_p.pvx<<"\t"<<_p.pvy<<"\t"<<_p.pvz;
   // os<< _p.pos[0]<<"\t"<<_p.pos[1]<<"\t"<<_p.pos[2]<<"\t"<<_p.vel[0]<<"\t"<<_p.vel[1]<<"\t"<<_p.vel[2];
     return os;                 
}
