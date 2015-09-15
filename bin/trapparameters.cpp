#include "trapparameters.h"
_trap_param::_trap_param() {
    Vrf = 50;
    Vdc = 5;
    kappa = 0.001;
    r0 = 0.04;
    z0 = 0.04;
    freq_rf = 50e5;
    trap_type = 0;
}

_trap_param::~_trap_param() {}
