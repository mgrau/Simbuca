//fieldcalc.h
#ifndef FIELDCALC_H
#define FIELDCALC_H


#include <iostream>
using namespace std;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <cstdlib> // for exit function
#include <sstream>
#include "math.h"
#include <stdexcept>
#include "matrix.h"
#include <map>
#include <vector>
#include "globals.h"
#include "numrec.h"
#include "mtrand.h"
#include <time.h>

void InitFieldCalc(double parameter1, double parameter2);

void ChangePotential(double _t);

void CalcElectricField();
void CalcMagneticField();

vector<double> getElectrCalc(double _x,double _y,double _z);
vector<double> getMagnCalc(const double& _x,const double& _y,const double& _z);

#endif
