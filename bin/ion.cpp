//ion.cpp
#include "ion.h"

Ion::Ion(long double _mass_amu){
    mass = _mass_amu*amu; 
    charge=1;
}

Ion::Ion(long double _mass_amu, int _charge){
    mass = _mass_amu*amu; 
    charge=_charge;
}

Ion::Ion(){
    mass = 0.0*amu; 
    wc	= 0.0;
    wz2 = 0.0;
    wmin	= 0.0;
    wplus	= 0.0;
    kd_div_U = 0.0; 
    kq_div_U = 0.0; 
    ko_div_U = 0.0;  
    charge=1;           
}

void Ion::SetParameters(double _mass_amu, double _frac,_trap_param & trap_param){
    charge=1;
    double r_0 = trap_param.r_0;
    mass = _mass_amu*amu; 
    wc	= (charge*el_charge*trap_param.B0)/mass;
    wz2 = (charge*el_charge*trap_param.Ud2)/mass;					//de 2 duid op het kwadraad van w_z
    wmin	= (wc-(sqrt(wc*wc-2.0*wz2)))*0.5;		//benaderende waarde: Ud2/(2*B)
    wplus	= (wc+(sqrt(wc*wc-2.0*wz2)))*0.5;		//benaderende waarde: Ud2/(2*B)
    kd_div_U = charge*trap_param.a*(el_charge/mass)*1.0/(r_0); 
    kq_div_U = charge*2.0*trap_param.a*(el_charge/mass)*1.0/(r_0*r_0);      
    ko_div_U = charge*4.0*trap_param.a*(el_charge/mass)*1.0/(r_0*r_0*r_0*r_0);

    //fraction=_frac;
}

void Ion::SetParameters(double _mass_amu, double _frac, double _charge,_trap_param & trap_param){
    charge=_charge;
    mass = _mass_amu*amu; 
    double r_0 = trap_param.r_0;
    mass = _mass_amu*amu;
    wc	= (charge*el_charge*trap_param.B0)/mass;
    wz2 = (charge*el_charge*trap_param.Ud2)/mass;					//de 2 duid op het kwadraad van w_z
    wmin	= (wc-(sqrt(wc*wc-2.0*wz2)))*0.5;		//benaderende waarde: Ud2/(2*B)
    wplus	= (wc+(sqrt(wc*wc-2.0*wz2)))*0.5;		//benaderende waarde: Ud2/(2*B)
    kd_div_U = charge*trap_param.a*(el_charge/mass)*1.0/(r_0);
    kq_div_U = charge*2.0*trap_param.a*(el_charge/mass)*1.0/(r_0*r_0);
    ko_div_U = charge*4.0*trap_param.a*(el_charge/mass)*1.0/(r_0*r_0*r_0*r_0);

    //fraction=_frac;
}

void Ion::SetParameters(_trap_param & trap_param){
    double r_0 = trap_param.r_0;

    wc	= (charge*el_charge*trap_param.B0)/mass;
    wz2 = (charge*el_charge*trap_param.Ud2)/mass;					//de 2 duid op het kwadraad van w_z
    wmin	= (wc-(sqrt(wc*wc-2.0*wz2)))*0.5;		//benaderende waarde: Ud2/(2*B)
    wplus	= (wc+(sqrt(wc*wc-2.0*wz2)))*0.5;		//benaderende waarde: Ud2/(2*B)
    kd_div_U = charge*trap_param.a*(el_charge/mass)*1.0/(r_0);
    kq_div_U = charge*2.0*trap_param.a*(el_charge/mass)*1.0/(r_0*r_0);
    ko_div_U = charge*4.0*trap_param.a*(el_charge/mass)*1.0/(r_0*r_0*r_0*r_0);
}

void Ion::Printall(){
     cout<<mass<<" "<<wc<<" "<<wz2<<" "<<wmin<<" "<<wplus<<" "<<kd_div_U<<" "<<kq_div_U<<endl;
}


ostream& operator<<(ostream& os,Ion _ion){
     //os<< (_ion.Getmass() / amu);
     os << _ion.GetName();
     return os;                 
}

void operator<<(LogFile& lf, Ion _ion){
     lf<<(_ion.Getmass() / amu);         
}

long double& Ion::Getmass(){
       return mass;       
}

long double& Ion::Getwc(){
       return wc;
}
long double& Ion::Getwplus(){
       return wplus;
}
long double& Ion::Getwmin(){
       return wmin;
}
long double& Ion::Getwz2(){
       return wz2;
}
long double& Ion::Getkd_div_u(){
       return kd_div_U;
}
long double& Ion::Getkq_div_u(){
       return kq_div_U;
}
long double& Ion::Getko_div_u(){
       return ko_div_U;  
}

int& Ion::Getcharge(){
       return charge;  
}
