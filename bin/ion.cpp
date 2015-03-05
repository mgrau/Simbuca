#include "ion.h"

Ion::Ion(long double _mass_amu) {
    mass = _mass_amu*amu; 
    charge = 1;
}

Ion::Ion(long double _mass_amu, int _charge) {
    mass = _mass_amu*amu; 
    charge = _charge;
}

Ion::Ion(){
    mass = 0.0*amu; 
    charge = 1;           
}

void Ion::SetParameters(double _mass_amu) {
    mass = _mass_amu*amu; 
    charge = 1;
}

void Ion::SetParameters(double _mass_amu, double _charge) {
    mass = _mass_amu*amu; 
    charge = _charge;
}

void Ion::Printall() {
    cout<<mass<<" "<<charge<<endl;
}

ostream& operator<<(ostream& os,Ion _ion) {
    os << _ion.GetName();
    return os;                 
}

void operator<<(LogFile& lf, Ion _ion) {
    lf<<(_ion.Getmass() / amu);         
}

long double& Ion::Getmass() {
    return mass;       
}
int& Ion::Getcharge() {
    return charge;  
}
