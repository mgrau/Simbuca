//coll.h.

#ifndef COLL_H
#define COLL_H

#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include "globals.h"
#include "mtrand.h"
#include <time.h>

bool gas_interact(double &vx, double &vy, double &vz, const double mass, const double timestep);             
//returns true if there was a collision
//timestep = h
int& GetNrCollisions();

void InitColl();

void SetPressure(double _pressure_bar);

void ChangeBuffergasType(double mass_amu, double radii_m, double polarization_factor);

void CheckBuffergasCreation();
#endif
