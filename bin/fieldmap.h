//fieldmap.h
#ifndef FIELDMAP_H
#define FIELDMAP_H


#include "matrix.h"
#include "SLogger.h"
#include <map>
#include <vector>

#include <iostream>
using namespace std;


class fieldmap{
private:
	double r_min,r_max,r_step,z_min,z_max,z_step;
    double factor;
	double nrows;
	double ncols;
	Matrix<std::pair<double,double> > * matrix;     //(Ez,Er)


	double r,dzr;
	double z_index,r_index;
	double z1,z2,r1,r2;
	int z1_index,z2_index,r1_index,r2_index;

	pair<double, double> m12,m22,m11,m21, mzr;
	
	bool IsElement(double _value);
	inline double round(double d);
	bool interpolatepoints;

  public:
	fieldmap(int _nrows, int _ncols);
	~fieldmap();
	
	void ReadField(char * Frz_filename);
	void ReadField(char *_Fr_filename,char *_Fz_filename);
	void Print();
	void PrintOnAxis();
    inline void SetFactor(double _f){factor = _f;};
    inline void SetRmin(double _rmin){r_min = _rmin;};
    inline void SetZmin(double _zmin){z_min = _zmin;};
	vector<double> getField(const double& _x,const double& _y,const double& _z);
    vector<double> getFieldXY(const double& _x,const double& _y);
	void SetInterpolate(double _interpolate); // (true) interpolate the fielmap points yes or no...
};					     

#endif
