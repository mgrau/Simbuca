#ifndef IONFLY_H
#define IONFLY_H


using namespace std;

#include "globals.h"
#include "particle.h"
#include "ioncloud.h"
#include "force.h"
#include "ode.h"
#include "ImageCharges.h"
#include "logfile.h"
#include "trapparameters.h"
#include "SLogger.h"

#include <math.h>
#include <limits>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>

#ifdef __GUI_ON__
#include <QObject>
#include "../GUI/counter.h"
#endif

#ifdef __GUI_ON__
void SetPercentagePointer(Counter *c);
#endif


void MoveParticles(double _time_movement,IonCloud &_cloud,_ode_vars &odev);

//Structors
void InitIonFly(const char * _filenamebegin,int _ode_order, double _timestep, bool _adaptive_stepsize,IonCloud &_cloud,_ode_vars & odev); //InitstheFiles its uses.
void ExitIonFly(IonCloud &_cloud);


void AddParticle(double x, double y, double z, double vx, double vy, double vz,Ion &_ion, IonCloud &_cloud);
void DelParticle(int _index, IonCloud &_cloud);

void DoNoExcitation(double _time_movement, bool _buffergas, double _p_buffergas,IonCloud &_cloud,_ode_vars & odev);
void ChangeEfieldmap(char * _trapErz);
void DoTransfer(double _time_movement, bool _buffergas, double _p_buffergas,char * _transferErz, char * _trapErz,IonCloud &_cloud,_ode_vars & odev);

// ACCESSOR
void SetCoulomb(bool _coulombinteraction);
void UseScaledCoulomb(double _ScaledCoulombFactor);
//Dipole Excitations
void DoDipoleExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_d,IonCloud &_cloud,_ode_vars & odev);
void DoDipoleExcitationWithBuffergas(double _time_movement, long double _exc_w, double _U_d, double p_buffergas,IonCloud &_cloud,_ode_vars & odev);

//Quadrupole Excitations
void DoQuadrupoleExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_q,IonCloud &_cloud,_ode_vars & odev);
void DoQuadrupoleExcitationWithBuffergas(double _time_movement, long double _exc_w, double _U_q, double p_buffergas,IonCloud &_cloud,_ode_vars & odev);

//Octupole Excitations
void DoOctupoleExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_q,IonCloud &_cloud,_ode_vars & odev);
void DoOctupoleExcitationWithBuffergas(double _time_movement, long double _exc_w, double _U_q, double p_buffergas,IonCloud &_cloud,_ode_vars & odev);

//Rotating wall Excitations
void DoRotatingWall(int _order,double _time_movement, long double _exc_w, double _U_exc, bool _buffergas, double _p_buffergas,IonCloud &_cloud,_ode_vars & odev);
void DoAntiRW(int _order,double _time_movement, long double _exc_w, double _U_exc, bool _buffergas, double _p_buffergas,IonCloud &_cloud,_ode_vars & odev);

// Axial coupling excitations
void DoAxialCouplingExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc,IonCloud &_cloud,_ode_vars & odev);

// SIMCO
void DoSIMCOWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, long double _exc_w2, double _U_exc2,IonCloud &_cloud,_ode_vars & odev);
// NOT CODED YET
void DoAxialQuadCoulpingExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, long double _exc_w2, double _U_exc2,
        long double _exc_w3, double _U_exc3,long double _exc_w4, double _U_exc4,IonCloud &_cloud,_ode_vars & odev);

void DoARexcitation(double _time_movement, double _U_exc,double _exc_w, double _exc_w2, IonCloud & _cloud,_ode_vars & odev);

void DoFBexcitation(double _time_movement, double _U_exc,double _exc_w, double _exc_w2, IonCloud & _cloud,_ode_vars & odev);

//Scans
void DoFrequencyScan(double _time_movement, double _U_exc, int _order, bool _rotatingwall, double _p_buffergas, double _central_freq, double _freq_dev, int _nr_steps,IonCloud &_cloud,_ode_vars & odev);

void IncludeElectrodeBoundaries(bool _bool);


void SetPrintInterval(double _print_interval);


void SetTotalTime_of_Simu(double t_,_ode_vars & odev);
void Print_particles(IonCloud &_cloud);
//eventualle, we can change this later together with setPrintInterval in one function he!
void use_particle_file(bool _bool,IonCloud &_cloud); //standard negative

void SetLifetime(double _lifetime, IonCloud &_cloud); // Set the lifetime of the ioncloud


void SetBufferGas(bool b_);

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w);

#endif
