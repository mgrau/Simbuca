//fieldcalc.cpp

#include "fieldcalc.h"

const double a=1.8/1000.; // [mm]
const double d=0.14/1000.; // [mm]

const double l1=20./1000.; // [mm]
const double l2=1.33094/1000.; // [mm]
const double l3=0.440187/1000.; // [mm]


const double z[10]={0. , l1 , l1+d , l1+d+l2 , l1+d+l2+d , l1+d+l2+d+l3 , l1+d+l2+d+l3+d , l1+d+l2+d+l3+d+l2 , l1+d+l2+d+l3+d+l2+d , l1+d+l2+d+l3+d+l2+d+l1};
const double lambda = z[9]; //=43.6621 [mm]
const double lambdainv=1./lambda;
double V[6]={0.,0.,0.880428,1,0.880428,0.}; // [V] note V[0] is not used.
//double V[6]={0.,0.,0.0,0.0,0.0,0.}; // [V] note V[0] is not used.

const double B0=1.20156; // [T]
double B2=301171; // [T/m^2]

double noise_sigma;

vector<double> E_field(3);
vector<double> B_field(3);


//random number
MTRand calcrand;

void InitFieldCalc(double parameter1, double parameter2){
     calcrand.seed(time(0)+23); //collsion random number generator
     B2=parameter1;
     noise_sigma=parameter2;
}


double RandGauss(double mean,double stddev){// create randoms with gaussian
// distribution mean = mean, standard deviation = stddev
         double bla1, bla2, woe, boe1;
          do {
                 bla1 = 2.0 * calcrand() - 1.0;
                 bla2 = 2.0 * calcrand() - 1.0;
                 woe = bla1 * bla1 + bla2 * bla2;
         } while ( woe >= 1.0 || woe == 0.);

         woe = sqrt( (-2.0 * log( woe ) ) / woe );
         boe1 = bla1 * woe;
         return (boe1*stddev + mean);
}


void ChangePotential(double _t){
	for(signed int i=0; i<6;i++){
		V[i] = V[i]+RandGauss(0.,noise_sigma);
	}	
}



void CalcElectricField(){
//initialize calculation for Efield

	//calculate axial field in z and r direction as a check
	ofstream tmp;
	tmp.open("tmpEz.txt");
	for(double z=-5.0/1000. ; z < 5.0/1000.; z += 0.1/1000.){
		E_field = getElectrCalc(0,0,z);
		tmp<<z<<"\t"<<E_field[2]<<endl;	
	}
	tmp.close();
	tmp.open("tmpEr.txt");
	for(double r=0.0/1000. ; r < 1.8/1000.; r += 0.05/1000.){
		E_field = getElectrCalc(r,0,0);
		tmp<<r<<"\t"<<E_field[0]<<endl;	
	}
	tmp.close();
}


void CalcMagneticField(){
//initialize calculation for Bfield
	//calculate axial Bfield in z and r direction as a check
	ofstream tmp;
	tmp.open("tmpBz.txt");
	for(double z=-5.0/1000. ; z < 5.0/1000.; z += 0.1/1000.){
		B_field = getMagnCalc(0,0,z);
		tmp<<z<<"\t"<<B_field[2]<<endl;	
	}
	tmp.close();
	tmp.open("tmpBr.txt");
	for(double r=0.0/1000. ; r < 1.8/1000.; r += 0.05/1000.){
		B_field = getMagnCalc(r,0,0);
		tmp<<r<<"\t"<<B_field[0]<<endl;	
	}
	tmp.close();
}


vector<double> getElectrCalc(double _x,double _y,double _z){
	double gamma,kn,tmp;
	double Er = 0.;
	double Ez = 0.;	
	double r=sqrt(_x*_x+_y*_y);
	//double potential = 0.;
	for(signed int n=1; n< 100; n++){
		kn = n*pi/lambda;
		tmp = 0;
		for(signed int i =1; i<5; i++){
			tmp += (V[i+1]-V[i])/(kn*kn*d)*(sin(kn*z[2*i])-sin(kn*z[2*i-1]));
		}
		gamma = 2./(lambda * BESSI0(kn*a)) * ( (V[1]*cos(kn*z[0])-V[5]*cos(kn*lambda))/kn + tmp ); 
		//potential += gamma*BESSI0(kn*r)*sin(kn*(_z+lambda*0.5));
		Er += -gamma*kn*BESSI1(kn*r)*sin(kn*(_z+lambda*0.5));
		Ez += -gamma*kn*BESSI0(kn*r)*cos(kn*(_z+lambda*0.5));	
	}//for(n)
	//cout<<" potential at point "<<_x<<" "<<_y<<" "<<_z<<" = "<<potential<<" V\n";
	E_field[0]=Er * _x/r;	 
	E_field[1]=Er * _y/r; 
	E_field[2]=Ez;
	return E_field;
}


vector<double> getMagnCalc(const double& _x,const double& _y,const double& _z){
	double br,bz;	
	double r = sqrt(_x*_x+_y*_y);
	bz = B0 + B2*(_z*_z-r*r*0.5);
	br = -B2* r*_z;

	B_field[0]=br*(_x/r);
        B_field[1]=br*(_y/r);
        B_field[2]=bz;

	return B_field;
}




