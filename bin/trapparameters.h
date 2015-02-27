
//
//  trapparameters.h
//
//
//  
//
//

////////////////////////////////////////////////////////////////////////
//
//  Trap parameters manager class
//
//  correction factor of the electric potential:
//  phi = U0*c0
//      + U0*c2*(z^2 - 0.5 r^2)
//      + U0*c4*(z^4 - 3*r^2*z^2 =3/8 r^4)
//      + U0*c6*(z^6 - 15/2*r^2*z^4 + 45/8*r^4*z^2 - 5/16*r^6)
//
//  correction factor of the magnetic field
//  Bz = B0*(1 - B2*(z^2 - r^2/2) )
////////////////////////////////////////////////////////////////////////
#ifndef TRAPPARAMETERS_H // header guards
#define TRAPPARAMETERS_H
#include <sstream>
#include <string>
#include <algorithm>    // std::find_if
#include <fstream>
#include <iostream>
#include "math.h"
#include "SLogger.h"

using namespace std;
class _trap_param{
public:
    
    _trap_param();
    ~_trap_param();
    
    int ReadFile(string filename); // return 1 if error
    void Print();
    //magnetic field
    double B0;
    double B2;
    // electrode radius
    double r_0;
    //electric potential
    double Ud2; // U/d^2 defines w_z2 = qU/(d^2 m)
    double U0;

    double c2; //   w_z2 = 2 q U0 c2/(d^2 m) -> c2 = 1/(2 d^2)
    double c4;
    double c6;
    double a;
    double d2; // = ((z_0^2 +r_0^2/2)/2)^(1/2)
    // diaphragm
    double PUMPING_DIAPRHAGM_RADIUS;

    // trap config number
    //0 IDEAL TRAP
    //1 NONIDEAL TRAP
    //2 REAL TRAP i.e. EM field map
    int trap_config;
    
};
#endif
