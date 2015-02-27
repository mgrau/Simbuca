//Initialization.h
#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "ion.h"
#include "ioncloud.h"
#include "ode.h"
#include "SimParser.h"
#include <vector>
#include <iostream>
void Run(int argc, char * argv[],SimParser & sparser);
void Run(SimParser & sparser);
void InitCloud(int nparticles, double eV_max_boltz,int seed,double semiaxis[3], double offset[3], std::vector<Ion> Ions, IonCloud &_cloud);

#endif
