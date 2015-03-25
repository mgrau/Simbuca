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

void move_particles(double _time_movement,IonCloud &_cloud,_ode_vars &odev);

void InitIonFly(const char * _filenamebegin,int _ode_order, double _timestep, bool _adaptive_stepsize,IonCloud &_cloud,_ode_vars & odev); //InitstheFiles its uses.

void ExitIonFly(IonCloud &_cloud);

void AddParticle(double x, double y, double z, double vx, double vy, double vz,Ion &_ion, IonCloud &_cloud);

// void DelParticle(int _index, IonCloud &_cloud);

void normal_operation(double time, IonCloud &_cloud,_ode_vars & odev);

void SetCoulomb(bool _coulombinteraction);

void UseScaledCoulomb(double _ScaledCoulombFactor);

void SetPrintInterval(double _print_interval);

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w);

#endif
