//force.h
//still to make a newer version off the cosinus!

#ifndef FORCE_H
#define FORCE_H




#include "globals.h"
#include "trapparameters.h"
#include "particle.h"
#include "ion.h"
#include "ioncloud.h"
#include "fieldmap.h"
#include "fieldcalc.h"
#include "logfile.h"
#include "ode.h"
#include <iostream>


using namespace std;
#include <vector>
//#include <math.h>




void force(const IonCloud &_cloud,_ode_vars &odev);
//returns the acceleration, not the force Pay attention!!

//void InitExcitationVars(int _order, bool _rotatingwall, bool _antirotatingwall, double _U_exc,double _w_exc);

void VV_BField_off(); //required for the Velocity Verlet integrator


void NoIdealTrap();
void NoIdealTrap(char *_Erz_filename, char *_B_filename);
void NoIdealTrap(char *_Er_filename, char *_Ez_filename, char *_B_filename);
void ChangeEfield(char *_Erz_filename);

void InitARexcitation(char *Erz_filename, const char * filename_begin);
void InitFBexcitation(char *Erz_filename, const char * filename_begin);
void TmpSetVars(double _wcm,double _wz2m,double _kq,double _kd);

//void SetCoulomb(bool _coulombinteraction);

//void UseScaledCoulomb(int _ScaledCoulombFactor);


#endif // FORCE_H
