//Ion.h
#ifndef ION_H
#define ION_H
#include "globals.h"
#include <math.h>
#include "logfile.h"
#include <iostream>
using namespace std;



struct Ion
{
  long double mass; // reads in amu, stored in kg
  int charge; // read and stored as units of charge
  //double fraction;
  string name;

  Ion();
  Ion(long double _mass_amu); //by default the charge is positive
  Ion(long double _mass_amu, int _charge);
  
  void Printall();
  void SetParameters(double _mass_amu);
  void SetParameters(double _mass_amu, double _charge);

  long double& Getmass();
  int& Getcharge();
  inline string GetName(){return name;};
  inline void SetName(string name_){name= name_;};
};

ostream& operator<<(ostream& os, Ion _ion);
void operator<<(LogFile& lf, Ion _ion);

#endif

