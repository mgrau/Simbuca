#ifndef FORCE_H
#define FORCE_H




#include "globals.h"
#include "trapparameters.h"
#include "particle.h"
#include "ion.h"
#include "ioncloud.h"
#include "logfile.h"
#include "ode.h"
#include <iostream>


using namespace std;
#include <vector>

void force(const IonCloud &_cloud,_ode_vars &odev);
//returns the acceleration, not the force Pay attention!!

void TmpSetVars(double _wcm,double _wz2m,double _kq,double _kd);


#endif // FORCE_H
