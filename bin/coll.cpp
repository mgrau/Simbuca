//coll.cpp
#include "coll.h"

//evironment vars
int T_n = 293;				
int T = 293;                //actually I also use T_buff now, quite confusing...
//float K0 = 21.5e-4;			//voor K in He
float p_n=1.013e5 ;  //1 atm pressure in pascal
float p = 1e-7*1e5;			//Buffergas Pressure [Pascal] 1bar=1e5 Pascal

//conserning the Buffergas
//for ion-properties http://www.chemicool.com/elements/helium.html
double m_gas = amu*massa_He4_amu;
double diameter_gas = 2.0*atomic_radii_He4; //neutral buffergas


const double diameter_ion = 2.0*1.88E-10;//0.93E-10 diameter of the molecule 2.20 A for K, 0.31 A for He

double alfa_e = 2.433E-41;//by 1atm and 140degrees C  has to be in [C*m*m/V]  
double correction_factor = 1.0;
 
// 2.433E-41 ->Liquid helium      http://www.physics.uq.edu.au/people/jones/phys6040/tut04/node4.html
// 7.4383e-39

/*atomic radii He: 0.31

  atomic radii K39: 2.27A
  ionic radii K39: 1.52A
  
  atomic radii C35: 0.97A
  ionic radii C35:  xxx
        
  atomic radii Ar35:0.88A
  ionic radii AR35:  xxx
*/

//gas variables
//allround global vars, be sure these are initiated before calling em.
double randcoll, probcoll, theta_gas, phi_gas, vmb, randmb, v_coll, gas_n, gas_sigma, 
gas_s, gas_lambda, vrel_bar, mum, toni, n;
double gas_v[3];



//Simion collision models related variables.

//mean free path and s
double mfp,s = 0.0;

// Collision-cross section (m^2)
// (The diameter of the cross-sectional area is roughly
//  the sum of the diameters of the colliding ion and buffer gas particles.)
// (2.1E-19 is roughly for two Helium atoms--Atkins1998-Table 1.3)
// (Note: the Van der Waals radius for He is 140 pm = 1.40 angstrom.
//   -- http://www.webelements.com and http://en.wikipedia.org/wiki/Helium --
//   i.e. 2.46e-19 collision cross-section)
// (2.27E-18 is for collision between He and some 200 amu ion with combined
//  collision diameter of 2 + 15 angstroms.  It is used in some benchmarks.)
//so sigma_m = pi*(radii_buffergas+radii_atom)^2
double _sigma_m2=2.1e-19;

//-- mean gas speed (m/s)
double c_bar_gas;
//-- median gas speed (m/s)
double c_star_gas;   
double c_bar_rel;

int number_coll=0;
//pool
double pvx,pvy,pvz,pmass,pdiameter; //still have to construct pdiameter XXX

double spr=1;
double v1,v2;
//random number
MTRand collrand;


/////////////////////Buffergasstuff/////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


void InitColl(){
     number_coll=0;
     //-- Compute mean gas speed (m/s)
     c_bar_gas = sqrt(8.0*kb*T/pi/(m_gas));
     //-- Compute median gas speed (m/s)
     c_star_gas = sqrt(2.0*kb*T/(m_gas));     
     //cross section it pi*r_buffergas^2+pi*r_ion^2
     _sigma_m2 = pi*(atomic_radii_He4+diameter_ion*0.5)*(atomic_radii_He4+diameter_ion*0.5); 
     //cout<<"sigma changed in "<<_sigma_m2<<endl;
     collrand.seed(time(0)+23); //collsion random number generator
}

//double erf(double z){
/*-- Error function (erf).
--   erf(z) = (2/sqrt(pi)) * integral[0..z] exp(-t^2) dt
-- This algorithm is quite accurate.  It is based on
-- "Q-Function Handout" by Keith Chugg:
--   http://tesla.csl.uiuc.edu/~koetter/ece361/Q-function.pdf
-- See also http://www.theorie.physik.uni-muenchen.de/~serge/erf-approx.pdf
-- I also find that the following makes a reasonable approximation:
--   1 - exp(-(2/sqrt(pi))x - (2/pi)x^2)

    double z2 = fabs(z);
    double t = 1.0 / (1.0 +0.32759109962 * z2);
    double res = (    -1.061405429) * t;
    res = (res +1.453152027 ) * t;
    res = (res -1.421413741 ) * t;
    res = (res +0.2844966736) * t;
    res =((res -0.254829592 ) * t) * exp(-z2*z2);
    res = res + 1.0;
    if (z < 0.0){ res = -res;}
    return res;
}
*/

double gaussian_random(){
//--using the box muller algorithm
	s=1;
     do{
       v1 = 2*collrand()-1;   
       v2 = 2*collrand()-1;

       spr=v1*v1+v2*v2; 
     }
     while(spr >= 1);
     //cout<<v1<<" kkk "<<v2<<" ppp "<<s<<endl;
     double test = v1*sqrt(-2.0*log(spr)/spr);
     return test;
}

void create_random_gas(){
     
     //see Petersson how I achieved this (from IsCool)
     vmb = sqrt(-2.0*kb*T/m_gas*log(1-sqrt(collrand()))); 
     //for velocity
     theta_gas = collrand()*2.0*pi;
     phi_gas = acos(collrand()*2.0-1.0);
	 
     gas_v[0] = vmb*cos(theta_gas)*sin(phi_gas);
     gas_v[1] = vmb*sin(theta_gas)*sin(phi_gas);
     gas_v[2] = vmb*cos(phi_gas);
   //cout<<"normal: "<<gas_v[0]<<"  "<<gas_v[1]<<"  "<<gas_v[2]<<endl;   
       
     // for location
     theta_gas = collrand()*2.0*pi;
     phi_gas = acos(collrand()*2.0-1.0);
}    

/*
void create_random_gas(){
     double randmb = 0;
     double scale = 
     do{
     vmb = collrand()*10000;//maxwell boltzmann speed
     randmb = collrand();
     }
     while(randmb < len/scale);
     
     //for velocity
     theta_gas = collrand()*2.0*pi;
     phi_gas = 2.0*pi*collrand() - pi;
     gas_v[0] = vmb*cos(theta_gas)*sin(phi_gas);
     gas_v[1] = vmb*sin(theta_gas)*sin(phi_gas);
     gas_v[2] = vmb*cos(phi_gas);
  	// for location
     theta_gas = collrand()*2.0*pi;
     phi_gas = 2.0*pi*collrand() - pi; 
}
*/

void create_random_gas_HS1(){
     double vr_stdev_gas = sqrt((kb*T)/(m_gas));

     double v_ion = sqrt(pvx*pvx+pvy*pvy+pvz*pvz);
     double scale = v_ion +vr_stdev_gas *sqrt(3.0)*3;
     double len = 0.0;
     do{
        gas_v[0] = vr_stdev_gas*gaussian_random();
        gas_v[1] = vr_stdev_gas*gaussian_random();
        gas_v[2] = vr_stdev_gas*gaussian_random();
	//cout<<"HS1: "<<gas_v[0]<<"  "<<gas_v[1]<<"  "<<gas_v[2]<<endl;cin.get();
	len = sqrt((gas_v[0]-pvx)*(gas_v[0]-pvx)+(gas_v[1]-pvy)*(gas_v[1]-pvy)+(gas_v[2]-pvz)*(gas_v[2]-pvz));
     }
     while(collrand() < len /scale);
     
     
     // for location
     theta_gas = collrand()*2.0*pi;
     phi_gas = 2.0*pi*collrand() - pi; 
}

void do_collision(){
     
     // particle number 'particle_index' experienced a collision
     double vrel;
     double vx1, vy1, vz1, vx2, vy2, vz2;
     double vx1r, vy1r, vz1r;
     double st, ct, sp, cp;
     
     do{
       
       do{
    	//normal way to do it
        create_random_gas();
       //(useless) SIMION HS1 way to do it
       //create_random_gas_HS1(); 

       
	                               vx2 = gas_v[0];
                                       vy2 = gas_v[1];
                                       vz2 = gas_v[2];
                                       vrel = sqrt((vx2-vx1)*(vx2-vx1)
                                       + (vy2-vy1)*(vy2-vy1)
                                       + (vz2-vz1)*(vz2-vz1));
          }while (vrel==0);

     
     vx1 = pvx;
     vy1 = pvy;
     vz1 = pvz;   
     
     //boost reference frame so that gas particle ar rest
     vx1 = vx1 - vx2;
     vy1 = vy1 - vy2;
     vz1 = vz1 - vz2;
     
     //rotate reference frame so that gas particle on the z-axis
	 //theta_gas angle is from Z-axis toward xy-plane
	 //phi_gas angle is from X-axis towards YZ plane
     st = sin(theta_gas);
     ct = cos(theta_gas);
     sp = sin(phi_gas);
     cp = cos(phi_gas);
     
     vx1r=ct*cp*vx1 + ct*sp*vy1 - st*vz1;
     vy1r=cp*vy1 - sp*vx1;
     vz1r=st*cp*vx1 + st*sp*vy1 + ct*vz1;
     
     }while(vz1r<=0); //no collision (ion has to move towards gas particle which is on the z axis)     
     
     //collision
     vz1r = vz1r*((pmass-m_gas)/(pmass+m_gas));
     
     //rotate refernce frame back + unboost
     
     vx1=ct*cp*vx1r - sp*vy1r + st*cp*vz1r + vx2;
     vy1=ct*sp*vx1r + cp*vy1r + st*sp*vz1r + vy2;
     vz1=ct*vz1r - st*vx1r                 + vz2;
     
     //assign new velocity
     
     pvx = vx1;
     pvy = vy1;
     pvz = vz1;
     
}


bool collision(double timestep){
     // k = particle index
     mum = pmass*m_gas/(pmass+m_gas);
     n = (p/(kb*T)); //volume density
     toni= sqrt(alfa_e/mum);
     randcoll = collrand();
     probcoll = 1-exp(-timestep*toni*2.0*ke*el_charge*n*pi);
     //cout<<"sams method: "<<probcoll<<endl;
     if(probcoll > 0.1){ cout<<"BANG"<<endl;}
     if (randcoll<probcoll){return true;}
     else {return false;}
}


double K_n(const double &vd){//return 0.002;}
       //for heavy masses
       //return ((11590+vd)/(0.19008*(11590+vd)*(11590+vd)-3420*(11590+vd)+2.084e7)); //is the construction for heavy masses he!
       
       /*for light masses
       a=10720
       b0=2.0833e7
       b1=-3067.13
       b2=0.151211
       */
       return ((10720.0+vd)/(0.151211*(10720.0+vd)*(10720.0+vd)-3067.13*(10720.0+vd)+2.0833e7)); //is the construction for heavy masses he!       
       }    


bool collision_via_KO(const double &timestep, const double &v){
     double mum = pmass*m_gas/(pmass+m_gas);
     double vd2 = (pmass*v*v-3*kb*T)/(pmass+m_gas);
     if (vd2 <= 0 ){ vd2 = 0;}
     else{ vd2 = sqrt(vd2);}
     double K = K_n(vd2)*T*p_n/(T_n*p);
     //average collision rate:
     double avg_coll_rate = el_charge /(mum*K);
     randcoll = collrand();
     probcoll = 1-exp(-timestep*avg_coll_rate);
     //cout<<"K0 method: "<<probcoll<<endl;cin.get();
     if (randcoll<probcoll){return true;}
     else {return false;}
   }     
     
   
bool simion_collision_simple(double timestep, double velocity){
   //this collision is like described in the 2001 paper from A.D. Appelhans and D.A. Dahl: "Measurement of external ion infjection adn trapping efficiency in the ion trap mass
   //spectrometer and comparison with a predictive model"
   mfp = kb*T*0.707106781/(pi*2.3e-10*2.3e-10*p); //mfp = kt /(pi*d*d*P*sqrt(2) with d=2.3 angstrum P is the pressure

   randcoll = collrand();
   probcoll =1-exp(-timestep*velocity/mfp);
   if(randcoll<probcoll){return true;}
   else {return false;}
}   


   
bool simion_collision_hs1(double timestep, double velocity){
   //this collision is the HS1 model revision 4 from 2006 as delivered with SIMION V8.x.x (corect version?)

   //-- Compute mean relative speed (m/s) between gas and ion.
   s = velocity / c_star_gas;
   c_bar_rel = c_bar_gas * ((s + 1.0/(2.0*s)) * 0.5 * sqrt(pi) * erf(s) + 0.5 * exp(-s*s));
   //-- Compute mean-free-path (m)
   mfp = kb * T * (velocity / c_bar_rel) / (p * _sigma_m2);
   
   
   //-- Store data about this calculation.
   //last_speed_ion = speed_ion
   //last_ion_number = ion_number
   //--print("DEBUG:ion[c],gas[c_bar],c_bar_rel,MFP=",
   //--      speed_ion, c_bar_gas, c_bar_rel, effective_mean_free_path_mm)
   //-- Note: The following is a simpler and almost as suitable
   //-- approximation for c_bar_rel, which you may used instead:
   //-- c_bar_rel = sqrt(speed_ion^2 + c_bar_gas^2)
       
    /*-- Compute probability of collision in current time-step.
    -- > For an infinitesimal distance (dx) traveled, the increase in the
    --   fraction (f) of collided particles relative to the number
    --   of yet uncollided particles (1-f) is equal to the distance
    --   traveled (dx) over the mean-free-path (lambda):
    --     df/(1-f) = dx / lambda
    --   Solving this differential equation gives
    --     f = 1 - exp(- dx / lambda) = 1 - exp(- v dt / lambda)
    --   This f can be interpreted as the probability that a single
    --   particle collides in the distance traveled.
   */
   randcoll = collrand();
   probcoll =1.0-exp(-timestep*velocity/mfp);
   if(randcoll<probcoll){return true;}
   else {return false;}

}    
     
bool gas_interact(double &vx, double &vy, double &vz, const double mass, const double timestep){
        //assign to poolvalues
        pvx=vx;
        pvy=vy;
        pvz=vz;
        pmass=mass;
	 //if(simion_collision_hs1(timestep,sqrt(vx*vx+vy*vy+vz*vz))){
     //if(collision_via_KO(timestep,sqrt(vx*vx+vy*vy+vz*vz))){
     if(collision(timestep)){
		//cout<<"collision nr: "<<number_coll<<endl;
        //cout<<"initial velocity's: \n";
	    //cout<<"\t pvx: "<<pvx<<"\n \t pvy: "<<pvy<<"\n \t pvz: "<<pvz<<endl;
	    do_collision();
        number_coll++;
        //reassign poolvalues!
        vx=pvx;
        vy=pvy;
        vz=pvz;
        //cout<<"end velocity's: \n";
        //cout<<"\t pvx: "<<pvx<<"\n \t pvy: "<<pvy<<"\n \t pvz: "<<pvz<<endl;
		return true;
     }else{
        return false;
     }
}                  
                  
                  
int& GetNrCollisions(){
     return number_coll;
}


void SetPressure(double _pressure_bar){
     p = _pressure_bar*1e5;//so to pascal he!     
     //alfa_e=epsilon_0*kb*T*6.02553e-9*correction_factor;
     //cout<<"alfa_e changed in "<<alfa_e<<endl;
}

void ChangeBuffergasType(double mass_amu, double radii_m, double polarization_factor){
     m_gas = mass_amu*amu;
     diameter_gas = 2.0*radii_m; 
     correction_factor = correction_factor*polarization_factor; 
}

void CheckBuffergasCreation(){
     ofstream bgasfile;
     bgasfile.open("maxwellboltz_test.txt");
     bgasfile<<"vmb vx\t vy\t vz \t theta_gas \t phi_gas \n";
     for(int k=0; k <100000; k++){
        create_random_gas();
        bgasfile<<vmb<<"\t"<<gas_v[0]<<" "<<gas_v[1]<<" "<<gas_v[2]<<" "<<theta_gas<<"\t"<<phi_gas<<"\n";
     }
     bgasfile.close();        
}
