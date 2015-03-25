#ifndef FORCE_H
#define FORCE_H

#include "ode.h"

void force(const IonCloud &_cloud,_ode_vars &odev);
//returns the acceleration, not the force Pay attention!!

#endif // FORCE_H
