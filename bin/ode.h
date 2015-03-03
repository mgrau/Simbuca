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
#include "Potmap.h"
#include "fieldmap.h"
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


class _force_vars{
public:
    _force_vars();
    ~_force_vars();
    //vector<vector<double> > derivs;
    double (*derivs)[3];
    bool excitation_type[15];
    /*
    0 DIP
    1 QUAD
    2 OCT
    3 AXIAL (check the equation in force.cpp)
    4 RW
    5 Anti RW
    6 SIMCO
    7 DUMM
    8 AR
    9 FB
    10 EXC_EMAP 
     */
    double coswTtimeTU;
    double sinwTtimeTU;
    double coswTtimeTU2; // for SIMCO
    double coswTtimeTU3; // for SIMCO
    double coswTtimeTU4; // for SIMCO
    double cos2wTtimeTU; // for RW/AW
    double sin2wTtimeTU; // for RW/AW
    double scaledCoulombFactor ;
    bool coulombinteraction;
    // trapparamaeters
    _trap_param trap_param;
    PotMap Potential_map;
    fieldmap * exc_emap;
    
    void Reset_excitation_type();
    inline void Set_excitation_type(int _exc){excitation_type[_exc]=true;};
    void Load_EXC_EMAP(string _file_name,double _factor);
    void Reset_EXC_EMAP();
private:
};



class _ode_vars {
public:
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
    bool VV;
    bool change_stepsize;
    bool Gear_initialized;
    bool VV_initialized;
    double hnext; //pay attenion for hnext (well and for h, where it iss changed...)    
    /////////////////
    /////////////////
    //relative and absolute error definition
    double rk_atol;
    double rk_rtol;
    //Gear stepsize (1e-9 is precise enough. No difference with 1e-10)
    double Gdt;
    /////////////////
    /////////////////
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
    // force function
    void InitExcitationVars(int _order, bool _rotatingwall, bool _antirotatingwall, double _U_exc, double _w_exc);
    // force variable
    // excitation variable
    _force_vars forcev;
    // excitation
    double time_ini_ope;
    double U_exc, U_exc2, U_exc3, U_exc4;
    double w_exc, w_exc2,w_exc3, w_exc4;
    // charge
    bool withCharge;
    inline void SetwithCharge(bool b_){withCharge = b_;};
    inline bool GetwithCharge(){return withCharge;};
    // PRINT PARTICLES
    inline void SetPrintAfterOperation(bool b_){PrintAfterOperation=true;}
    bool PrintAfterOperation;
    bool PrintatZpos_bool;
    double PrintZpos;
    inline void SetPrintatZpos_bool(bool b_){PrintatZpos_bool = b_;};
    inline void SetPrintZpos(double z_){PrintZpos = z_;};

    // sweep
    bool sweep_flag; // flag for linear sweep
    
    double sweep_wi; // initial pulsation for sweep
    double sweep_wf; // final pulsation for sweep
    double sweep_par;
    // time of simu
    double initial_time; // initial time of the simu
    double final_time; // final time of the operation
    double total_time; // final time of the simulations
    //swift
    string swift_file_name;
    vector< double > swift_time;
    vector< double > swift_amp;
    vector< double > swift_amp_cos;
    vector< double > swift_amp_sin;
    bool swift_flag;
    bool swift_RW;
    inline void SetSWIFT(bool b_){swift_flag = b_;};
    inline void SetSWIFT_RW(bool b_){swift_RW = b_;};
    double swift_dt;
    int swift_nsample;
    double SetSWIFT_function(string nf); // return t_end
    double SetSWIFT_function_RW(string nf); // return t_end
    double Interpolate_SWIFT(double t_);
    double Interpolate_SWIFT(double t_,int sc);
    void UnsetSWIFT();
    
    
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
    
private:
    
};



void step(IonCloud &_cloud,_ode_vars &odev);
//tis should be followed by
//            particles_time +=h;
//            h= Gethnext(); //hnext is calculated in error stuff  

//void InitOde(int _ode_order, double _timestep, bool _adaptive_stepsize);
double GetTimeStep(_ode_vars &odev);
//int GetCountStep();
//make sure the stepsize h stays always minimum this size... can be removed
void SetMinStepsize(double _h_min);
void ResetOde();


#endif
