//ode.h
/*        Calculate Dormand Prince for all particles together and the error for all the particles together 
          inside the Dormand_Prince_rk function then calculate the hnext and assign it to the next h.
          Make sure beta=0.08 otherwise you`re error becomes >1 and you have to calculate 
          Dormand Prince all over again! this is solving the equations for all the particles, AGAIN!
          see note: 

          The beta=0.08 (because 5th order rk) is the control mechanism to hold h fluently. 
          If beta is 0.04 or even 0 then your timestep(h) is altered by letting your error
          become bigger and just then (approx around every 10 cycles your err>1) h is set 
          gradually smaller.
          The control method with b=0.08 is preferable, since h changes rather smoothly this way.
          In general, when putting beta=0.08 your step stays around 10-9 so you should have no problems
          at all! Great :-)
          */
#ifndef ODE_H
#define ODE_H

#include "globals.h"
#include "particle.h"
#include "ion.h"
#include "ioncloud.h"
#include "mpi_funcs.h"
#include "trapparameters.h"
#include "SLogger.h"

#ifdef __NBODY_ON__
#include "nbody.h"
#endif //__NBODY_ON__

#ifdef __CUNBODY_ON__
#include "cunbody.h"
#endif //__CUNBODY_ON__

#include "logfile.h"
#include <vector>
#include <math.h>
#include <limits>
#include <iostream>
using namespace std;

struct _force_vars {
    _force_vars();
    ~_force_vars();
    double (*derivs)[3];

    double coulomb_scale;
    bool coulomb_interaction;

    bool trap;
    bool tof;

    void reset_ops();
    _trap_param trap_param;
};

struct _ode_vars {
    _ode_vars();
    ~_ode_vars();
    // ode function
    void Init_ode(int _ode_order, double _timestep, bool _adaptive_stepsize);
    void Reset_ode();
    void Initpoolvectors(int nrParticles);
    // ode vars
    double httry;
    double h;

    bool RK4;
    bool DP5;
    bool change_stepsize;
    bool Gear_initialized;
    double hnext; //pay attenion for hnext (well and for h, where it iss changed...)    
    /////////////////
    //relative and absolute error definition
    double rk_atol;
    double rk_rtol;
    //Gear stepsize (1e-9 is precise enough. No difference with 1e-10)
    double Gdt;
    /////////////////
    double EPS;        
    double errold; //standard also defined
    bool reject; //standard false!
    bool poolvectorsInitialized;      //standard false he    
    double reset_timestep;        
    //Pool Particles 
    //vector< vector<double> > k1,k2,k3,k4,k5,k6,k7;
    double (*k1)[6],(*k2)[6],(*k3)[6],(*k4)[6],(*k5)[6],(*k6)[6],(*k7)[6];
    vector< double > yerr;
    vector<vector<double> >  accel_new_suggestion, a_diff;
    unsigned i,j;
    double sk,err;
    //from RK4...		
    //just the same as k1,2,...
    vector< vector<double> > bn,cn,dn,en,fn,b,c,d,e,f;
    double h2f,h3f,h4f,h5f,h6f;	 
    int countstep;
    // force variable
    _force_vars forcev;
    double time_ini_ope;
    // charge
    bool with_charge;
    // PRINT PARTICLES

    bool print_after_operation;
    bool print_at_x;
    double print_x;

    // time of simu
    double initial_time; // initial time of the simu
    double final_time; // final time of the operation
    double total_time; // final time of the simulations

    // mpi funcs
#ifdef __MPI_ON__
    //pointeur -> can be create and delete when I want
    _mpi_vars  * mpiv;
    //used in Operations::Launch();
    inline void create_mpiv(){mpiv = new _mpi_vars;};
    inline void delete_mpiv(){delete mpiv;};
#endif // __MPI_ON__

    // NBODY ALGO
#ifdef __CUNBODY_ON__
    _cunbody cunbody;
#endif //__CUNBODY_ON__

#ifdef __NBODY_ON__
    _nbody nbody;
#endif //__NBODY_ON__
};

void step(IonCloud &_cloud,_ode_vars &odev);
double GetTimeStep(_ode_vars &odev);
void SetMinStepsize(double _h_min);
void ResetOde();

#endif
