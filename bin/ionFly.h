#ifndef IONFLY_H
#define IONFLY_H
#include <math.h>
#include <limits>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include "globals.h"
#include "particle.h"
#include "ioncloud.h"
#include "force.h"
#include "ode.h"
#include "logfile.h"
#include "trapparameters.h"
#include "SLogger.h"
using namespace std;

#ifdef __GUI_ON__
#include <QObject>
#include "../GUI/counter.h"
#endif

#ifdef __GUI_ON__
void SetPercentagePointer(Counter *c);
#endif

void MoveParticles(double _time_movement,IonCloud &_cloud,_ode_vars &odev);

void InitIonFly(const char * _filenamebegin,int _ode_order, double _timestep, bool _adaptive_stepsize,IonCloud &_cloud,_ode_vars & odev); //InitstheFiles its uses.
void ExitIonFly(IonCloud &_cloud);


void AddParticle(double x, double y, double z, double vx, double vy, double vz,Ion &_ion, IonCloud &_cloud);
void DelParticle(int _index, IonCloud &_cloud);

void DoNoExcitation(double _time_movement, bool _buffergas, double _p_buffergas,IonCloud &_cloud,_ode_vars & odev);
void ChangeEfieldmap(char * _trapErz);

void SetCoulomb(bool _coulombinteraction);
void UseScaledCoulomb(double _ScaledCoulombFactor);

void IncludeElectrodeBoundaries(bool _bool);

void SetPrintInterval(double _print_interval);

void SetTotalTime_of_Simu(double t_,_ode_vars & odev);
void Print_particles(IonCloud &_cloud);

void use_particle_file(bool _bool,IonCloud &_cloud); //standard negative

void SetLifetime(double _lifetime, IonCloud &_cloud); // Set the lifetime of the ioncloud

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w);

#endif
