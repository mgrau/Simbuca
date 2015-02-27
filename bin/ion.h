//Ion.h
#ifndef ION_H
#define ION_H
#include "globals.h"
#include "trapparameters.h"
#include "math.h"
#include "logfile.h"
#include <iostream>
using namespace std;



class Ion
{
  long double mass; // *amu so in kg  XXX
  long double wc, wplus, wmin, wz2; //as angular frequency
  long double kd_div_U, kq_div_U, ko_div_U; //amplitude parameters divided by the amplitude(is choosen later on he)
  int charge;
  //double fraction;
  string name;
public:
  //constructor
  Ion(long double _mass_amu);  //by default the charge is positive
  Ion(long double _mass_amu, int _charge);
  Ion(); //default constructor, default euuh...
  //maybe I should also make a copy constructor... default copys all parameters off the thing he!
  
  void Printall();
  void SetParameters(double _mass_amu, double _frac,_trap_param & trap_param);
  void SetParameters(double _mass_amu, double _frac, double _charge, _trap_param & trap_param);
  void SetParameters(_trap_param & trap_param);

  long double& Getmass();
  long double& Getwc();
  long double& Getwplus();
  long double& Getwmin();
  long double& Getwz2();
  long double& Getkd_div_u();
  long double& Getkq_div_u();
  long double& Getko_div_u();
  int& Getcharge();
  //inline double GetFraction(){return fraction;};
  //inline void SetFraction(double f_){fraction = f_;};
  inline string GetName(){return name;};
  inline void SetName(string name_){name= name_;};
};

ostream& operator<<(ostream& os, Ion _ion);
void operator<<(LogFile& lf, Ion _ion);

#endif

