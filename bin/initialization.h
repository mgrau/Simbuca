#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include <vector>
#include <iostream>
#include "ion.h"
#include "ioncloud.h"
#include "ode.h"
#include "SimParser.h"
#include "inireader.h"

void Run(INIReader & parser);
void InitCloud(int nparticles, double eV_max_boltz,int seed,double semiaxis[3], double offset[3], std::vector<Ion> Ions, IonCloud &_cloud);
#endif
