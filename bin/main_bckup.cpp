/*known errors:
        in parser.cpp: the radii taken is always the radii of Ar354 0.88E-10m
                       one can make a case loop for all mases it reads from the file.
*/


#include "MPI_simbuca.h"


#include <iostream>
//using namespace std;
#include <vector>
#include <fstream>
#include <iterator> 
#include <math.h>
#include "globals.h"
#include "trapparameters.h"
#include "ionFly.h"
#include "ion.h"
#include "fieldmap.h"
#include "mtrand.h"
#include "logfile.h"
#include "parser.h"
#include "IonTable.h"
#include "Operation.h"
#include "Initialization.h"
#include <sstream>

#ifdef __NBODY_ON__
#include "nbody.h"
#endif //__NBODY_ON__




////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//Generally, segmentation faults are because a pointer of yours is either
//NULL, or points to random memory (probably never initialized to anything)
//, or points to memory that was deleted.
LogFile logger;





/* to simulate positrons or protons. Default the charge is positive, but in the definition of ion the charge can be given as a 2nd parameter as well) */
Ion positron(massa_electron_amu);
Ion antiproton(massa_proton_amu,-1);
Ion proton(massa_proton_amu);





int main( int argc,  char* argv[] ){
	
	/*************** initialize MPI ****************/
#ifdef __MPI_ON__
	//cout << "MPI ON "<< endl;
	double startwtime = 0.0, endwtime;
	int    namelen;
	char   processor_name[MPI_MAX_PROCESSOR_NAME];	
	MPI::Status status;
	MPI::Init(argc,argv);
	int myid = MPI::COMM_WORLD.Get_rank();
	int numprocs = MPI::COMM_WORLD.Get_size();
	MPI::Get_processor_name(processor_name,namelen);
	
	//std::cout << "Process "<< myid<< " of "<< numprocs << " is on " << processor_name <<std::endl;

    MPI::COMM_WORLD.Barrier();

	if (myid == 0)
		startwtime = MPI_Wtime();
#endif // __MPI_ON__
	
	// SIMBUCA FILE INPUT
    	string s_arg;
    	string extension_simbuca = ".sim";
    	size_t found_ext;
    	// if argc==2 -> check if argv[1] is the simbuca file
  
    IonCloud Cloud;
    _ode_vars Ode_vars;
    
#ifdef __MPI_ON__
    Ode_vars.create_mpiv();
#endif // _MPI_ON__
    	if(argc==2)
    	{
           
        	s_arg = argv[1];
        	found_ext=s_arg.rfind(extension_simbuca);
        	if (found_ext!=string::npos)
        	{ 
            		DoSimulation_with_input_file(s_arg, Cloud, Ode_vars);
	    		       
        	}
    	}
	else // argc==2
	{
	
	
#ifdef todo
	/**************************************/
	
	//initialize the number of particles to be simulated (note: when running on a GPU also change this in force_gpu.cpp
	int NRPARTICLES = 1;
	//initialize buffer gas pressure (in mbar)
	double  p_buff_mbar=1e-4;
	//set the order of the integration method. Either use 1 or 5 here. 
	const int ODEORDER = 5;
	bool adaptative_stepsize = true;
	double timestep = 1e-9;
	double printinterval_timestep = 1e-8;
	// Coulomb interaction
	bool Coulomb_enable = false;
	double Coulomb_scale_factor = 1;
        // Print all particles in one file or separate
	bool separate_particle_file = true;	
	
	/******* initialize GPU *********************/ 
	
#ifdef __NBODY_ON__
#ifdef __MPI_ON__  
	if(Coulomb_enable)
	{
	
		if(myid==0)
		{ 
			initialize_gpu(NRPARTICLES,0,0,0);
		} 
	}	
#else   
	if(Coulomb_enable)
	{
		initialize_gpu(NRPARTICLES,0,0,0);
	}
#endif //__MPI_ON__
	
#endif //     __NBODY_ON__     
	
	/**************************************/    	
	
	//initialize the prefix of the filenames
#ifdef __MPI_ON__
	char filename_prefix[250];
	stringstream ssltemp2;ssltemp2<<"test" << "-" << myid;
	//stringstream ssltemp2;ssltemp2<<"test" << myid;
	ssltemp2>>filename_prefix;
#else // __MPI_ON__  
	const char * filename_prefix = "myfirstsim";
#endif // __MPI_ON__   
	
	//initialize the random seed
	unsigned long seed;	
	seed = 0;

	// open logger 
	stringstream ssltemp;ssltemp<<filename_prefix<<"_logfile.txt";
	char filename_logger[250];ssltemp>>filename_logger;
	logger.open(filename_logger);
	
	//Print out the particle information in
	//true: each particle has his own file (default)
	//false: print out all information in one .._pAll.txt
	use_particle_file(separate_particle_file);
	
	
	//Initialize the ODE. First the filename_prefix, then the Order of the ODE can be 1,4 or 5, next is the initial timestep (should be around 1e-9), next is a boulean
	//true= with adaptive stepsize, false is without adaptive stepsize.
	InitIonFly(filename_prefix,ODEORDER,timestep,adaptative_stepsize);
#ifdef __MPI_ON__
    MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__  

	//print out the ions informations to the file every certain timestep
	SetPrintInterval(printinterval_timestep);
	//with or without Coulomb Interaction
	SetCoulomb(Coulomb_enable);
	//scaled Coulomb Factor, if used
	UseScaledCoulomb(Coulomb_scale_factor);
	//with or without including the electrode Boundaries. I.e. an ion is lost if it hits the electrode walls
	// electrode Boundaries not tested with MPI + GPU !
	IncludeElectrodeBoundaries(false);
	//In case you want to include external fieldmaps
	//NoIdealTrap("Er_cooler_trapping.txt","Ez_cooler_trapping.txt","6T fieldmap z-35to35 r0to25.txt");
	NoIdealTrap("../fieldmaps/Erz_FI_well.txt","../fieldmaps/B_fieldmap500x500_80Amps.txt");

        //Import the data from previous simulations
	//Or create the ions according to a maxwell boltzmann distribution with the maxima around 1.5eV 
	//MPI_ImportData("Rb512-dip");	
	double semiaxis[3]={0.001,0.001,0.001};
	double offset[3]={0.,0.02,-0.15};
	vector<Ion > Ions;
	//Ion Rb87(massa_87_Rb_amu);
	//Rb87.SetFraction(25.);
   	//Ion Rb85(massa_85_Rb_amu);
	//Rb85.SetFraction(75.); 
	//Ions.push_back(Rb87);
	//Ions.push_back(Rb85);
	Ions.push_back(antiproton);
	
	
	InitCloud(NRPARTICLES,1.5,seed,semiaxis,offset, Ions);
	//ImportData("import",Table);

	

#ifdef __MPI_ON__ 
	MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__
	
	/* What simulations do you want to do? */

	//DoNoExcitation(1e-3,false,p_buff_mbar/1000.0); 
	//DoDipoleExcitationWithBuffergas(0.025, Rb87.Getwmin(), 0.05, p_buff_mbar/1000.0);
    	//DoQuadrupoleExcitationWithBuffergas(0.075, Rb87.Getwc(), 0.115, p_buff_mbar/1000.0);
	
        //If you want to simulate a transfer between 2 penning traps. For advanced users only :)
	//DoTransfer(0.000065, false, 0.0,"transfer_Er.txt","transfer_Ez.txt","traps_Er.txt","traps_Ez.txt"); 
	
	

	//Here you can put what you want to test...
	CalcElectricField();



	/* close all outputfiles */ 	
	ExitIonFly();   
	
	/* PostProces the data */
	//   BundleData(filename_prefix,NRPARTICLES);
	//email_it("G100_bundled.txt" , "blabla@blablab.com");
	//see http://cc.byexamples.com/20070120/print-color-string-without-ncurses/
#ifdef __MPI_ON__ 
	if(myid==0)
#endif // __MPI_ON__	
	printf("%c[%d;%dmProgam Ended%c[%dm\n",27,4,31,27,0);   
	
    /* clear GPU memory */
    
#ifdef __NBODY_ON__ 
	if(Coulomb_enable)
	{ 
#ifdef __MPI_ON__  
	if(myid==0) finalize_gpu();
#else //__MPI_ON__
	finalize_gpu();
#endif  //__MPI_ON__
	}
#endif //     __NBODY_ON__
	
	
	
    /***************  finalize  MPI ****************/
#ifdef __MPI_ON__
	MPI::COMM_WORLD.Barrier();
	if(!separate_particle_file)
 		MPI_Concatenate_Data(NRPARTICLES,filename_prefix);
	MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__
#endif // 0	
	}  // argc!=2 .sim file 
	
#ifdef __MPI_ON__	
	MPI::COMM_WORLD.Barrier();
    if (myid == 0) 
	{
		endwtime = MPI_Wtime();
		std::cout << "wall clock time = " <<  endwtime-startwtime << std::endl;
	}
    MPI::Finalize();
#endif // __MPI_ON__
    /******************************************/
    
	return 0;
}




