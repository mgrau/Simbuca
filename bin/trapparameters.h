#ifndef TRAPPARAMETERS_H // header guards
#define TRAPPARAMETERS_H
#include <string>

using namespace std;
struct _trap_param {
    _trap_param();
    ~_trap_param();
    
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

    // other operational parameters
    double E_kick;

    double newVrf;
    double newVdc;
};
#endif
