//ImageCharges.h
#include "ImageCharges.h"

    
    //calculate the image charge & current
    const double elec_z1 = -0.125; //U5: -0.125 -0.095        U8: -0.233 -0.193
    const double elec_z2 = -0.095; 
    const double elec_gamma = 1.3; //correction factor; roughly scaled such that Approx method is similar to the Bessel one

    //MRE electrode size
    const double mreLength = 0.526;
    const double electrodeRadius = 0.04;
    const double trapCenter = -0.032;

    const double electrodeLeft = elec_z1 - trapCenter;        //shifted position of the left edge of the pickup electrode
    const double electrodeRight = elec_z2 - trapCenter;        //shifted position of the right edge of the pickup electrode


    //to calulate the EED BASE dimensions
    const double lr = 0.00092 ; // ring electrode
    const double lc = 0.00285; // correction electrode
    const double le = 0.01 ;  // endcap electrode
    const double d = 0.0007; //electrode spacer
    const double a0 = 0.0035; // trap diameter

    //calculate the image charge & current BASE dimensions
/*    const double elec_z1 = lr+d/2; //-0.127;
    const double elec_z2 = lr+d/2+lc; //-0.127+0.03;
    const double elec_R =  a0;
    const double elec_gamma = 1.; //correction factor;
*/
    //variables for the tank circuit
    const double freq = 13.2e6; //resonance frequency
    const double Lx = 0.511e-6; // Henry
    const double Cx = 3e-12; // Farad
    const double Rx = 2*pi*freq*Lx; // effective resistance of the external circuit [Ohm]
    double IQold, dWold, iL, told, dT;

////////////////////////////////////////////////////////////
///////////////////// ImageCharges /////////////////////////
////////////////////////////////////////////////////////////


void InitImageCharges(int & nrparticles){
	IQold = 0;
        dT = 0;
}


// calculate the Effective Electrode Distance (EED) for axial detection in a cylindrical Penning trap
// argument is the length of the measurement electrode. The measurement electrode starts next to the ring electrode (lc/2+d) and stops at (lc/2+d+el_meas).

double EEDz(double el_meas){
      	double EED = 0.;
	double lambda = 2*(lc+le+2*d)+lr;
	//calculate sum
	for (int n=1; n < 100; n++){
		EED+=sin(n*pi/lambda*el_meas)*sin(n*pi/lambda*(2*d+lr+el_meas))/BESSI0(2*n*pi/lambda*a0);
	}
	EED = 1./EED;//^ -1
	EED = lambda/8.*EED;//*lambda/8
	return EED;
}



double GetImageChargeApprox(int & charge, double & z){
	double zdown = elec_z2-z;
	double zup = elec_z1-z;
  	//compute image charge information!
	return  (-1*elec_gamma*el_charge*charge*0.5*( (zdown/sqrt(electrodeRadius*electrodeRadius+zdown*zdown)) - (zup/sqrt(electrodeRadius*electrodeRadius+zup*zup)) ) ); //multiplication wtih -1 singe i_image = -1* i_ring;
}


//---------------------------------------------------------------------------------------------------------------------
//  precise Image Charge via Bessel functions
//  returns induced image charge on the pickup electrode for the considered particle with (charge, Zpos)
double GetImageChargeBessel(int charge, double r, double z)
{
    //debug: force r=0

    double shifted_z = z - trapCenter;               // shift coordinate system to the center of the trap

    // max number of spatial modes to be considered:
    // as we approximated the particles to be on axis of the trap we already implemented a source of error
    // therefore we need to decide which order fits best to the error we made (######## TO BE DONE ########)
    const unsigned int orderOfSum = 100;

    //storage variables
    //double sinL1, sinL2, cosL1, cosL2;
    double qN = 0.0; //initialization to zero bec we sum over it

/*
    for(int n = 1; n <= orderOfSum; n++)
    {
        sinL1 = sin(electrodeLeft*(2*n-1)*PI/mreLength);
        sinL2 = sin(electrodeRight*(2*n-1)*PI/mreLength);

        cosL1 = cos(electrodeLeft*2*n*PI/mreLength);
        cosL2 = cos(electrodeRight*2*n*PI/mreLength);

        //
        qN += cos(((2*n-1)*PI/mreLength)*shifted_z) *  BESSI0(r*(2*n-1)*PI/mreLength) /( BESSI0(electrodeRadius*(2*n-1)*PI/mreLength)*(2*n-1)*PI ) * ( sinL1 - sinL2 );
        qN -= sin((2*n*PI/mreLength)*shifted_z) * BESSI0(r*(2*n-1)*PI/mreLength) /( BESSI0(electrodeRadius*2*n*PI/mreLength)*2*n*PI ) * ( cosL1 - cosL2 );

    }
*/
    for(int n = 1; n <= orderOfSum; n++)
    {

        qN += sin(2*pi*n*shifted_z/mreLength)*sin(n*pi*(electrodeRight-electrodeLeft)/mreLength)*sin(n*pi*(electrodeLeft+electrodeRight)/mreLength)* BESSI0(r*2*n*pi/mreLength) /( BESSI0(electrodeRadius*2*n*pi/mreLength) * n);

    }


    //returns the image charge induced by the particle to the U5 pickup electrode
    //return (2*charge/electrodeRadius)*qN;
    return -4.0*charge*qN/pi;
}

//-----------------------------------------------------------------------------------------------------------------------



// MC gets image current and converts it into a voltage difference
// MC given by the passage of the image current through the external circuit
void TankCircuit(double t, double IQ, double &dV){
  //dW = 0;
  //ImI = 0;


  dT = t-told;
  //dV = Rx*IQ;

  double dW = (dWold*(Cx/dT-dT/2/Lx)+IQ-iL)/(1/Rx+Cx/dT+dT/2/Lx);
  IQold = IQ;
  iL = iL+1/Lx*(dW+dWold)*dT/2;
  dWold = dW;
  dV = 12900*dW;
  /*
  if (stepp%10==0)
  {
   FILE *dvfile;
   dvfile = fopen( "imageI.dat", "a");
   if ( dvfile)
   {
    fprintf(dvfile, "%.4e %.6e %.6e\n", t*1e6, IQ, dW);
   }
   fclose ( dvfile);
  }// */
}


