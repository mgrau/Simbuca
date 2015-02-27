//CreateRecoilDistr.cpp
#include "CreateRecoilDistr.h"


#define SQRT2 1.414213562373095048

double M,F,R,smalla,smallb,me,E0;
double Rmax,Pmax;

double Temp, charge_state, halflife, mass,tob;

char beta_decay = '+'; //beta + or beta - 


MTRand mrand;



double Pr(double R){
       
       double Prob = 0.0;
       R = R*charge_state;
       Prob = 1.0/6.0*M*F*(-1.0-smalla)*

       (
        (
         ( ((E0-R+SQRT2*sqrt(M*R))*(E0-R+SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R+2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        *( ((E0-R+SQRT2*sqrt(M*R))*(E0-R+SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R+2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        *( ((E0-R+SQRT2*sqrt(M*R))*(E0-R+SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R+2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        )
       -(
         ( ((E0-R-SQRT2*sqrt(M*R))*(E0-R-SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R-2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        *( ((E0-R-SQRT2*sqrt(M*R))*(E0-R-SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R-2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        *( ((E0-R-SQRT2*sqrt(M*R))*(E0-R-SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R-2.0*SQRT2*sqrt(M*R)+2.0*me)+me )     	)
       )
       
       +0.25*M*F*(smalla*(E0+me)-smallb+E0+me)*
       (
        (
         ( ((E0-R+SQRT2*sqrt(M*R))*(E0-R+SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R+2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        *( ((E0-R+SQRT2*sqrt(M*R))*(E0-R+SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R+2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        )
        -(
          ( ((E0-R-SQRT2*sqrt(M*R))*(E0-R-SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R-2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
         *( ((E0-R-SQRT2*sqrt(M*R))*(E0-R-SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R-2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        )
       )
       +0.5*M*F*(0.5*smalla*(2.0*M*R+me*me-(E0+me)*(E0+me))+smallb*(E0+me))*
       (
        (
         ( ((E0-R+SQRT2*sqrt(M*R))*(E0-R+SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R+2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        -(
          ( ((E0-R-SQRT2*sqrt(M*R))*(E0-R-SQRT2*sqrt(M*R)))/(2.0*E0-2.0*R-2.0*SQRT2*sqrt(M*R)+2.0*me)+me )
        )
       )
       )
      ;       
       return Prob;
}

double getPmax(){
       
       double step = 1.0E-6; //1eV
       double Pmax = 0.0;
       double Prob = 0.0;
       double R = 0.0;
       while(R<=Rmax){
                      Prob=Pr(R);
                      if(Prob>Pmax){Pmax = Prob;}
                      R+=step;
                      }
       return Pmax;
       
}

double Erec(){
       double E;
       double P,y;
       do{
          E = mrand()*Rmax;
          P = Pr(E);
          y = mrand()*Pmax*1.5;
         }
       while(y>P);
       return E;
       }
       
double TOB(){
       double x;
       x = mrand();
       x = -(halflife/0.693)*log(x);
     return x;      
}

//bool cut_off(double energy){
//     bool cut_off = false;
//     double chance = mrand();
//	 double U = Ud2 *0.04*0.04; //this is the endcapvoltage
//     if(chance<sqrt(charge_state*U/(energy))) cut_off = true;     
//     
//     return cut_off;
//     }
     

void CreateRecoilSpectraDifferential(char * filename){
//writes a recoil spectrum to a file 
//and generates (If I feel like it) smalla gnuplot file 
  
	
   ofstream dif_file(filename); 
  double energyloop = 0.0;
  for(energyloop = 0.0;energyloop<=Rmax ; energyloop+=1e-6)
    dif_file<<energyloop*1000000<<"\t"<<Pr(energyloop)<<"\n";               
  cout<<"Differential Recoil Spectra created in "<<filename<<endl;
  dif_file.close();
}

void GenerateRandomEnergy(char * filename,int N){
//writes out a file with N random energies written to it.
 ofstream dif2_file(filename); 
  double loop = 0.0;
  for(loop = 0.0;loop<=N ; loop++)
    dif2_file<<Erec()*1e6<<"\n";               
  cout<<"Endpoint "<<Rmax*1000000<<" eV.\n";
  cout<<"Random Recoil Energys created in "<<filename<<endl;
  dif2_file.close();
}


double MonteCarloEnergy(){
//throws back a certain energy value in eV. Accordingly to the recoil energy spectra.
	return Erec()*1e6;
}


void InitRecoilDistr(){
    mrand.seed(time(0));
    me=0.510998;
    Temp=273;
    F = 1;
    smallb = 0;
    smalla = 0.9004;
    charge_state=1;
  
  //type the mass here...
  int species = 61; 
  //case loop to put the vars in position.
  switch ( species )
      {
         case 35: //35Ar
 		mass = 34.975257;
  		E0=5.9653;
  		halflife = 1.775;
  		smalla=0.9004;
  		beta_decay = '+';
	     break;

         case 61: //Mn61
                mass = 60.94465;
  		E0=7.178;
  		halflife = 0.67;
  		smalla=-1/3;
  		beta_decay = '-';
	     break;

         case 62: //Mn62
            	mass = 61.94843;
  		E0=10.697;
  		halflife = 0.625;
  		smalla=-1/3;
  		beta_decay = '-';
	     break;

         case 63: //Mn63
            	mass = 62.95024;
  		E0=8.750;
  		halflife = 0.275;
  		smalla=-1/3;
  		beta_decay = '-';
	      break;	    
         
	 default:
            cout<<"choose a proper ion species to create the recoil spectra\n";
      }
  

     //beta + or -
        M = mass*mass_Mev;
	if(beta_decay == '+'){
		Rmax = ((E0-me)*(E0-me)-me*me)/(2.0*charge_state*M);}
	else{
		Rmax = ((E0+me)*(E0+me)-me*me)/(2.0*charge_state*M);}
	//important line to initialize everything.

	Pmax = getPmax();      
   //print out file information
      switch ( species )
      {
        case 35: //Ar35
 		CreateRecoilSpectraDifferential("differential_spectrum_35Ar.txt"); 
            break;
	case 61: //Mn61
                CreateRecoilSpectraDifferential("differential_spectrum_Mn61.txt"); 
            break;
	case 62: //Mn62
            	CreateRecoilSpectraDifferential("differential_spectrum_Mn62.txt"); 
            break;
	case 63: //Mn63
            	CreateRecoilSpectraDifferential("differential_spectrum_Mn63.txt"); 
            break;
      }
     GenerateRandomEnergy("randomenergies.txt",20000);
     
}

