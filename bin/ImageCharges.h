//coll.h.

#ifndef IMAGECHARGES_H
#define IMAGECHARGES_H

#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include "globals.h"
#include "mtrand.h"
#include "numrec.h"


void InitImageCharges(int & nrparticles);

double EEDz(double el_meas);
double GetImageChargeApprox(int & charge, double & z);
double GetImageChargeBessel(int,double,double);

void TankCircuit(double t, double IQ, double &dV);
#endif
