//
//  Potmap.h
//
//
//
//
//




#ifndef POTMAP_h
#define POTMAP_h
#include <cstdlib> // for exit function
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include "matrix.h"
#include "SLogger.h"


using namespace std;
class PotMap{
    
public:	
    PotMap();
    ~PotMap();
   void ReadEPotential(string Potmap_name);

   int GetIndex(const double & r,const double & z);
    double GetPotential(const double &x,const double &y, const double & z);
    double GetEr(const double &x,const double & y, const double & z);
    double GetEz(const double &x,const double & y, const double & z);
    pair<double,double > GetE(const double &x,const double & y, const double & z); //Er,Ez
    vector<double > GetEField_from_Pot(const double &x,const double & y, const double & z); //Ex,Ey,Ez
private:
    vector<double > potential;
    vector<double > rp;
    vector<double > zp;
	double dr;
    double dz;
    double rmin;
    double rmax;
    double zmin;
    double zmax;
    int numr;
    int numz;
    int maxIndex;
    
    bool verbose;
    
    // interpolator
    vector<double > ap;
    vector<double > fp;
};

#endif // POTMAP
