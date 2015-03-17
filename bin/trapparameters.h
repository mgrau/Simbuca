#ifndef TRAPPARAMETERS_H // header guards
#define TRAPPARAMETERS_H
#include <sstream>
#include <string>
#include <algorithm>    // std::find_if
#include <fstream>
#include <iostream>
#include "math.h"
#include "SLogger.h"

using namespace std;
struct _trap_param {
    _trap_param();
    ~_trap_param();
    
    int ReadFile(string filename); // return 1 if error
    void Print();
    // Voltages
    double Vrf;
    double Vdc;
    // geometry
    double kappa;
    double r0;
    double z0;
    // trap RF frequency, in Hz
    double freq_rf;
    // type, ideal=0
    int trap_type;
};
#endif
