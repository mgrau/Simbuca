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
#include "fieldcalc.h"
#include "mtrand.h"
#include "logfile.h"
#include "parser.h"
#include "IonTable.h"
#include "PDGTable.h"
#include "SimParser.h"
#include "Operation.h"
#include <sstream>
#include "Initialization.h"

#ifdef __linux
#include <sys/stat.h>
#include "errno.h"
#endif

#define is_power_of_two(x) ( ((x&(x-1)) == 0) && (x!=0) )

LogFile tmplogger;

MTRand algrand;
IonTable Table("../libraries/ame/mass.mas12");
PDGTable pdgTable("../libraries/pdg/pdg_table.txt");


int closest_sup_power_of_two(int n);

//containers
	vector< double > x_;
	vector< double > y_;
	vector< double > z_;
	vector< double > vx_;
	vector< double > vy_;
	vector< double > vz_;

int countword(string line)
{
    int wordcount;
    if (line.at(0)==' ')
		wordcount=0;
	else
		wordcount=1;
	for(int i=1; i<line.length() ;i++)
		if (line[i] != ' ' && line[i-1] == ' ')
            wordcount++;
	//cout<<"The string has "<<wordcount<<" words"<<endl;
    return wordcount ;
}

double RandomGaussian(double mean,double stddev){// create randoms with gaussian
// distribution mean = mean, standard deviation = stddev
         double bla1, bla2, woe, boe1;
          do {
                 bla1 = 2.0 * algrand() - 1.0;
                 bla2 = 2.0 * algrand() - 1.0;
                 woe = bla1 * bla1 + bla2 * bla2;
         } while ( woe >= 1.0 || woe == 0);

         woe = sqrt( (-2.0 * log( woe ) ) / woe );
         boe1 = bla1 * woe;
         //y2 = x2 * w;
         
         return (boe1*stddev + mean);
}

inline bool IsInteger(std::string & s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;
   char * p ;
   strtol(s.c_str(), &p, 10) ;
   return (*p == 0) ;
}

#ifdef __linux__
void mkpath(std::string s,mode_t mode){

    size_t pre=0,pos;
    std::string dir;
    int mdret;

    if(s[s.size()-1]!='/'){
        // force trailing / so we can handle everything in loop
        s+='/';
    }

    while((pos=s.find_first_of('/',pre))!=std::string::npos){
        dir=s.substr(0,pos++);
        pre=pos;
        if(dir.size()==0) continue; // if leading / first time is 0 length
        if((mdret=mkdir(dir.c_str(),mode)) && errno!=EEXIST){
           cout<<"error?"; 
        }
    }
}
#endif


void CreateCloud(int nparticles, vector<Ion> Ions,IonCloud &_cloud, _trap_param & trap_param){
	int j = 0;
	int myid=0;
	double frac_ =0;
#ifdef __MPI_ON__	
	//MPI
	myid = MPI::COMM_WORLD.Get_rank();
	int numprocs = MPI::COMM_WORLD.Get_size();
	int ii;
	// Broadcoast vector 
	MPI::COMM_WORLD.Bcast(&x_[0],x_.size(),MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(&y_[0],y_.size(),MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(&z_[0],z_.size(),MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(&vx_[0],vx_.size(),MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(&vy_[0],vy_.size(),MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(&vz_[0],vz_.size(),MPI::DOUBLE,0);
    

	frac_ = 0;
    int last_part = nparticles % numprocs;
    int * counts = new int[numprocs];
    int * displ = new int[numprocs];
    displ[0] = 0;
    for(int i=0;i<numprocs;i++)
    {
       counts[i] = nparticles/numprocs;
        if(i<last_part)
            counts[i] += 1;
    }

    for(int i=1;i<numprocs;i++)
    {
        displ[i] = displ[i-1] +counts[i-1];
    }

	for(int i=0;i < counts[myid] ;i++)
    	{
            ii = i + displ[myid];
		AddParticle(x_[ii],y_[ii],z_[ii],vx_[ii],vy_[ii],vz_[ii],Ions[ii],_cloud);
		//cout<<"particle added\t"<<x_[i] <<" "<< y_[i] <<" "<< z_[i] <<" "<< vx_[i] <<" "<< vy_[i] <<" "<< vz_[i] <<" "<< Ions[j]<<endl;cin.get();
    	}
    
    delete [] counts;
    delete [] displ;
   
    _cloud.UpdateIDs(nparticles);
   
#else // __MPI_ON__
	frac_ = 0;
	j =0;
	for(int i=0;i<nparticles;i++)
	{
		AddParticle(x_[i],y_[i],z_[i],vx_[i],vy_[i],vz_[i],Ions[i],_cloud);

		 //cout<<"particle added\t"<<x_[i] <<" "<< y_[i] <<" "<< z_[i] <<" "<< vx_[i] <<" "<< vy_[i] <<" "<< vz_[i] <<" "<< Ions[i]<<endl;cin.get();
	}
#endif // __MPI_ON__
	_cloud.UpdateIonParameters( trap_param);
}


void InitCloud(int nparticles, double eV_max_boltz,int seed,double semiaxis[3], double offset[3], vector<Ion> Ions, IonCloud &_cloud, _trap_param & trap_param)
{       /*
	Creation of the ions with a Maxwell Boltzmann distribution with ion positions in an ellipsoid 
        eV_max_boltz is the Energy in eV where there is the maximum of the Maxwell Boltzmann distibution.
 	note: The total speed 'vmb' is calculated with a Maxwell Boltzmann distribution (more or less. Check 		Petterson thesis about ISCool RFQ. 
	*/
	
	/*
	 the number of particles must be divible by numprocs!
	 */
	int myid=0;
#ifdef __MPI_ON__	
	//MPI
	myid = MPI::COMM_WORLD.Get_rank();
	int numprocs = MPI::COMM_WORLD.Get_size();
	int ii;
#endif // __MPI_ON__

	//data
	algrand.seed(seed);
	double mass;// = Cs133.Getmass();
	double vmb;
	
	//see Petersson how I achieved this (from IsCool). I checked it and it works fine :)
	// vmb = sqrt(-2.0*(eV_max_boltz*eV_to_Joule)/(mass)*log(1-sqrt(algrand())));
	
	double theta_gas;
	double phi_gas;
	
   	// containers
	x_.resize(nparticles);
	y_.resize(nparticles);
	z_.resize(nparticles);
	vx_.resize(nparticles);
	vy_.resize(nparticles);
	vz_.resize(nparticles);
   	int j = 0;
	double frac_ =0;
	
	if(myid==0)
	{		
		//printf("Initialization with %d particles\n",nparticles);	
		for(int i=0;i < nparticles;i++)
		{
			// positions
			while (true)   
			{    // three random values for one random point in 3D space   
				x_[i] = algrand()*2.-1.;   
				y_[i] = algrand()*2.-1.;
				z_[i] = algrand()*2.-1.;
				
				// the distance of x,y,z from 0,0,0 is sqrt(x*x+y*y+z*z), but   
				// as usual we can omit the sqrt(), since sqrt(1) == 1 and we   
				// use the 1 as boundary:   
				if ( x_[i]*x_[i] + y_[i]*y_[i] + z_[i]*z_[i] <= 1){break;}// found one valid point inside   
			} 
			x_[i]=offset[0]+semiaxis[0]*x_[i];
			y_[i]=offset[1]+semiaxis[1]*y_[i];
			z_[i]=offset[2]+semiaxis[2]*z_[i];

			mass = Ions[i].Getmass();
			//mass = proton.Getmass();
			// velocities	  
			
			theta_gas = algrand()*2.0*pi;
			phi_gas = acos(algrand()*2.0-1.0);
			//theta_gas = atan(x/y);
			vmb = sqrt(-2.0*(eV_max_boltz*eV_to_Joule)/(mass)*log(1-sqrt(algrand())));
			vx_[i] = vmb*cos(theta_gas)*sin(phi_gas);
			vy_[i] = vmb*sin(theta_gas)*sin(phi_gas);
			vz_[i] = vmb*cos(phi_gas);
			
		} // end loop
      
	} // end if master process 

	//create the cloud
	CreateCloud(nparticles,Ions,_cloud,trap_param);

	tmplogger<<nparticles;tmplogger<<" particles created. Max-Boltz distribution around ";tmplogger<<eV_max_boltz;tmplogger<<" eV. \n";
	tmplogger<<"spheroid cloud with semi-axis's = {";tmplogger<<semiaxis[0];tmplogger<<";";tmplogger<<semiaxis[1];tmplogger<<";";tmplogger<<semiaxis[2];tmplogger<<"}";  
	tmplogger<<" and offset = {";tmplogger<<offset[0];tmplogger<<";";tmplogger<<offset[1];tmplogger<<";";tmplogger<<offset[2];tmplogger<<"} ";
	tmplogger<<" \n";
	//if(myid==0)cout<<"cloud created\n";
} 


void DoSimulation_with_input_file(string sim_input_filename, IonCloud & _cloud,_ode_vars &odev)
{
    // MPI
    int myid =0;
#ifdef __MPI_ON__
  myid = MPI::COMM_WORLD.Get_rank();
  int numprocs = MPI::COMM_WORLD.Get_size();
#endif    
    //file
    ifstream gui_file;
    string comment;

    // PARAMETERS
    // PARTICLES
    int NRPARTICLES=0;
    int NSPARTICLES=0;
    int npart_test =0;
    vector<Ion > Ions_particles;
    vector<Ion > Ions_cloud;
    vector<Ion > Ions_import;

    Ion Ion_tmp;
    vector<string > ion_name;
    vector<string > ion_orig_name;
    vector<double > fraction_of_ion;
    vector<int> ion_charge;
    double frac_ =0 ;
    double fraction=0;
    // Config intitiatilsaion
    string init_config;
    stringstream filename_importdata;
    string filename_prefix_importdata;
    ifstream file_importdata;
    // Cloud
    double semiaxis_cloud[3];
    double offset_cloud[3];
    double energy=0;
    // GAS
    double p_buff_mbar=0;
    // ODE
    int ODEORDER=0;
    double timestep=0;
    int adaptative_timestep=0;
    bool adaptative_stepsize;
    int Coulomb=0;
    bool Coulomb_enable;
    double Coulomb_scale_factor=0;
    // output file
    string filename_prefix_;
    double printinterval_timestep=0;
    unsigned int seperateparticlefile=0;
    bool seperateparticlefile_bool;
    // trap
    int n_file_map=0;
    int tmp_int=0;
    string trap_config;
    string file_mapEr,file_mapEz,file_mapB;
    ifstream file_map;
    double parameter1,parameter2=0;
    // OPERATION
    string operation;
    string tmp_name;
    string subline;
    stringstream sstring_tmp;
    int _countword=0;
    unsigned int int_buffer_gas=0;
    string str_buffer_gas;
    double tmp_double=0;
    double tmp_amplitude=0;
    double tmp_frequency=0;
    double tmp_amplitude2=0;
    double tmp_frequency2=0;
    double tmp_amplitude3=0;
    double tmp_frequency3=0;
    double tmp_amplitude4=0;
    double tmp_frequency4=0;
    double tmp_f_bias=0;
    double tmp_f_bias2=0;
    double tmp_f_bias3=0;
    double tmp_f_bias4=0;
    double tmp_order=0;
    double relec_tmp=0;
    double Ud2_tmp=0;
    double B_tmp=0;
    //double tmp_frequency;
    string tmp_string;
    string orig_name;
    double time_operation=0;
    Operation Ope_tmp;
    Operations Ope;
   

    //IonTable Table("../libraries/ame/mass.mas12");
    if(myid == 0){
	cout << endl;
	cout << "\t*****************" << endl;
	cout << "\t* S I M B U C A *" << endl;
	cout << "\t*****************" << endl;
    }
    if(myid==0){
    	cout << "\nInitialization with Simbuca file" << endl;
	 cout << "*********************************************" << endl;
    }
    SimParser sparser(sim_input_filename.c_str());
    //   sparser.PrintParams();
 
    //----------------------------------------------
    // CLOUD
    if(sparser.mycreatecloud.flag){
      NRPARTICLES = sparser.mycreatecloud.cloudsize;
      if(sparser.mycloudparts.flag == false){
	cout << "Error: No cloud fractions provided!" << endl;
	exit(-1);
      }else{
	for(int i=0; i < sparser.mycloudparts.cloudfracs.size(); i++){
	  tmp_string = sparser.mycloudparts.cloudfracs[i].first;
	  orig_name = tmp_string;
	  frac_ = sparser.mycloudparts.cloudfracs[i].second;

	  //Read Ion properties
	  bool speciesfound = false;
	  double charge = 0;
	  
	  //check if the charge a charge is given (f.e. 35Ar3+) -> lastcharacter is a + or -
	  char lastchar = tmp_string[tmp_string.size()-1];
	  if(lastchar == '+' || lastchar == '-'){ 
	    unsigned int loopi;
	    for(loopi=tmp_string.size()-2; loopi > 0; loopi--){
	      if(!isdigit(tmp_string[loopi])) break;
	    }
	    string charge_string = tmp_string.substr(loopi+1,tmp_string.size()-loopi-2);
	    charge = atoi(charge_string.c_str());
	    if(lastchar == '-') charge *= -1;
	    tmp_string = tmp_string.substr(0,loopi+1);	
	    ion_charge.push_back(charge);			
	  }
	  else{ ion_charge.push_back(1);charge=1;} //no charge is given so charge +1
          
	  ion_name.push_back(tmp_string);
	  ion_orig_name.push_back(orig_name);
	  if((Table.TestIonName(tmp_string))){
	    tmp_double = Table.ResearchMass(tmp_string);
	    Ion_tmp.SetParameters(tmp_double,100.,charge,odev.forcev.trap_param);
	    Ion_tmp.SetName(tmp_string);
	    speciesfound = true;
	  }
	  
	  if(pdgTable.TestPDGName(orig_name)){
	    int pdgid = pdgTable.GetPDGId(orig_name);
	    tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	    Ion_tmp.SetParameters(tmp_double,100.,pdgTable.GetPDGCharge(pdgid),odev.forcev.trap_param);
	    Ion_tmp.SetName(orig_name);
  	    speciesfound = true;
	  }

	  if(!speciesfound){
	    if(myid==0)
	      cout << "ERROR in the names of mass species" << endl; 
	    exit(EXIT_FAILURE);	
	  }   
	  else{//particle found
	      //if(myid==0){
	      //  cout << Ion_tmp.GetName() << " " << Ion_tmp.Getmass()/amu << " ";
	      //  if (Ion_tmp.Getcharge() > 0)cout<<Ion_tmp.Getcharge()<<"+";
	      //  else cout<<Ion_tmp.Getcharge()*-1<<"-";
	     // cout<< " " << endl;
	     //}
	}

	  // add particles
	  int n = frac_*NRPARTICLES/100.;
	  for( int j= 0; j < n; j++){
		Ions_cloud.push_back(Ion_tmp);
	  }
	  ion_name.push_back(tmp_string);
	  ion_orig_name.push_back(orig_name);
	  fraction_of_ion.push_back(frac_);
	  fraction+=frac_;

	odev.SetwithCharge(false);
	for(int i=0;i<ion_charge.size();i++)
	  {
	    if(ion_charge[i]!=1)
	      {
		odev.SetwithCharge(true);
	      }
	    
	  }
	if(fraction>100)
	  {
	    if(myid==0)
	      cout << "ERROR in the composition of the cloud, fraction <" << fraction << ">" << endl; 
	    exit(EXIT_FAILURE);	
	    //return;
	  }
 	}//end for loop  
      }//end if loop      
      if(sparser.mycloudcoord.flag == false){
	if(myid == 0)
	  cout << "Error: No cloud coordinates provided!" << endl;
	exit(-1);
      }else{

	// Cloud parameters
	semiaxis_cloud[0] = sparser.mycloudcoord.semiaxis_cloud[0];
	semiaxis_cloud[1] = sparser.mycloudcoord.semiaxis_cloud[1];
	semiaxis_cloud[2] = sparser.mycloudcoord.semiaxis_cloud[2];
	offset_cloud[0]   = sparser.mycloudcoord.offset_cloud[0];
	offset_cloud[1]   = sparser.mycloudcoord.offset_cloud[1];
	offset_cloud[2]   = sparser.mycloudcoord.offset_cloud[2];
      }
      if(sparser.mytemp.flag == false){
	if(myid == 0)cout << "Error: No temperature provided!" << endl;
	exit(-1);
      }else{
	energy = sparser.mytemp.temp;
      }
    }

    //----------------------------------------------
    // IMPORTDATA
    if(sparser.myimportdata.flag){
      odev.SetwithCharge(false);
      filename_prefix_importdata = sparser.myimportdata.prev_simu_file;
#ifdef __MPI_ON__
#else // __MPI_ON__		
      int numprocs = 1;
#endif
      for(int k=0;k<numprocs;k++)
	{
	  filename_importdata << filename_prefix_importdata  << "_pAll.txt";
	  file_importdata.open(filename_importdata.str().c_str(),ios::in);
	  if(!file_importdata)
	    {
	      if(myid==0)
		cout << "ImportData file : " << filename_importdata.str() << " doesn't exist" << endl; 
	      exit(EXIT_FAILURE);
	      //return;
	    }
	  file_importdata.close();
	  filename_importdata.str("");
	}
    }

    //----------------------------------------------
    // PARTICLES
    if(sparser.mycreateparticles.flag){
      NSPARTICLES = sparser.mycreateparticles.nparts;
      if( sparser.myparticles.particles.size() != NSPARTICLES){
		cout<<"number of initialized particles in CREATEPARTICLES doesn`t correspond to number of input particles"<<endl;
		exit(-1);
      }
      if(sparser.myparticles.flag == false){
       if(myid == 0)
	  cout << "Error: no particles given!" << endl;
	exit(-1);
      }else{
	x_.resize(NSPARTICLES);
	y_.resize(NSPARTICLES);
	z_.resize(NSPARTICLES);
	vx_.resize(NSPARTICLES);
	vy_.resize(NSPARTICLES);
	vz_.resize(NSPARTICLES);
	int charge;
	bool speciesfound = false;
	
	for(int i = 0; i < sparser.myparticles.particles.size(); i ++){
	  string name = sparser.myparticles.particles[i].first;

	  vector<double> myvec = sparser.myparticles.particles[i].second;
	  speciesfound = false;
	  tmp_string = name;
	  orig_name = tmp_string;
	  
	  //check if the charge a charge is given (f.e. 35Ar3+) -> lastcharacter is a + or -
	  char lastchar = tmp_string[tmp_string.size()-1];
	  if(lastchar == '+' || lastchar == '-'){ 
	    unsigned int loopi;
	    for(loopi=tmp_string.size()-2; loopi > 0; loopi--){
	      if(!isdigit(tmp_string[loopi])) break;
	    }
	    string charge_string = tmp_string.substr(loopi+1,tmp_string.size()-loopi-2);
	    charge = atoi(charge_string.c_str());
	    if(lastchar == '-') charge *= -1;
	    tmp_string = tmp_string.substr(0,loopi+1);	
	    ion_charge.push_back(charge);			
	  }
	  else{ ion_charge.push_back(1);charge=1;} //no charge is given so charge +1
          
	  ion_name.push_back(tmp_string);
	  ion_orig_name.push_back(orig_name);
	  if((Table.TestIonName(tmp_string))){
	    tmp_double = Table.ResearchMass(tmp_string);
	    Ion_tmp.SetParameters(tmp_double,100.,charge,odev.forcev.trap_param);
	    Ion_tmp.SetName(tmp_string);
	    Ions_particles.push_back(Ion_tmp);
	    speciesfound = true;
	  }
	  
	  if(pdgTable.TestPDGName(orig_name)){
	    int pdgid = pdgTable.GetPDGId(orig_name);
	    tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	    Ion_tmp.SetParameters(tmp_double,100.,pdgTable.GetPDGCharge(pdgid),odev.forcev.trap_param);
	    Ion_tmp.SetName(orig_name);
	    Ions_particles.push_back(Ion_tmp);
	    speciesfound = true;
	  }
	  if(!speciesfound){
	    if(myid==0)
	      cout << "ERROR in the names of mass species" << endl; 
	    exit(EXIT_FAILURE);	
	  }   
	  else{ //particle found
	      //if(myid==0){
	      //  cout << Ion_tmp.GetName() << " " << Ion_tmp.Getmass()/amu << " ";
	      //  if (Ion_tmp.Getcharge() > 0)cout<<Ion_tmp.Getcharge()<<"+";
	      //  else cout<<Ion_tmp.Getcharge()*-1<<"-";
	      //cout<< " " << endl;
	     //}
	  }
	  x_[i] = myvec[0];
	  y_[i] = myvec[1];
	  z_[i] = myvec[2];
	  vx_[i] = myvec[3];
	  vy_[i] = myvec[4];
	  vz_[i] = myvec[5];
	  	  // cout<<x_[i] <<" "<< y_[i] <<" "<< z_[i] <<" "<< vx_[i] <<" "<< vy_[i] <<" "<< vz_[i] <<" "<< Ions_particles[i]<<endl;cin.get();
	}
		
      }
    }

    //----------------------------------------------
    // BUFFER GAS
    if(sparser.mybuffer.flag){
      int_buffer_gas = sparser.mybuffer.index;
      p_buff_mbar = sparser.mybuffer.pressure;
    }

    //----------------------------------------------
    //ODE
    if(sparser.myode.flag){
      ODEORDER = sparser.myode.order;
      adaptative_timestep = sparser.myode.adaptive;
      if(adaptative_timestep==0) 
	adaptative_stepsize = false;
      else
	adaptative_stepsize = true;

      timestep = sparser.myode.stepsize;
    }else{
      if(myid == 0)
	cout << "Error: no ODE configuration provided" << endl;
      exit(-1);
    }

    //----------------------------------------------
    // Coulomb
    if(sparser.mycoulomb.flag){
      Coulomb = sparser.mycoulomb.index;
      Coulomb_scale_factor = sparser.mycoulomb.weight;
      if(Coulomb==1) 
 	{
	  Coulomb_enable = true;
 	}
      else
	Coulomb_enable = false;
    }else{
      if(myid == 0)
	cout << "Error: no Coulomb configuration provided" << endl;
      exit(-1);
    }

    //----------------------------------------------
    //OUTPUT FILE
    if(sparser.myoutput.flag){
      filename_prefix_ = sparser.myoutput.outfile;
      printinterval_timestep = sparser.myoutput.sample_time;
      seperateparticlefile = sparser.myoutput.separate;
      if(seperateparticlefile == 1){
	seperateparticlefile_bool = true;}
      else{
	seperateparticlefile_bool = false;}
      if(seperateparticlefile == 2)
	{
	  seperateparticlefile_bool = false;
	  odev.SetPrintAfterOperation(true);
	}
        if(seperateparticlefile == 3)
        {
            seperateparticlefile_bool = false;
            odev.SetPrintatZpos_bool(true);
            odev.SetPrintZpos(sparser.myoutput.sample_time);
        }
        
    }else{
      if(myid ==0)
	cout << "Error: no Output file configuration provided!" << endl;
      exit(EXIT_FAILURE);
    } 



    if((sparser.mycreatecloud.flag == false)&&(sparser.myimportdata.flag == false)&&(sparser.mycreateparticles.flag == false)){
      if(myid==0)
	cout << "Error : Wrong cloud initialisation config " << endl;
      exit(EXIT_FAILURE);
      //return;
    }

    //----------------------------------------------
    // IDEAL TRAP
     if(sparser.myidealtrap.flag){
       odev.forcev.trap_param.r_0 = sparser.myidealtrap.r_electrode;
       odev.forcev.trap_param.Ud2 = sparser.myidealtrap.Ud2;
       odev.forcev.trap_param.B0 = sparser.myidealtrap.B;
     }       


    //----------------------------------------------
    // NONIDEAL TRAP
     if(sparser.mynonitrap.flag){
       sstring_tmp.str(sparser.mynonitrap.trapconfigfile);
       tmp_int = odev.forcev.trap_param.ReadFile(sstring_tmp.str());
       if(tmp_int==1){
	 exit(EXIT_FAILURE);
       }
       
       if(myid==0)
	 odev.forcev.trap_param.Print();
     }

    //----------------------------------------------
    // REAL TRAP
     if(sparser.myrealtrap.flag){
       odev.forcev.trap_param.trap_config = 3;
       n_file_map = sparser.myrealtrap.nfiles;
       if(sparser.myrealtrap.filenames.size() != n_file_map){
	 if(myid == 0)
	   cout << "Error: n of files and actual listed files are not consistent!" << endl;
	 exit(1);
       }
       if(n_file_map==2){	
	 file_mapEz = sparser.myrealtrap.filenames[0];
	 file_mapB = sparser.myrealtrap.filenames[1];	
	 file_map.open(file_mapEz.c_str(),ios::in);
	 if(!file_map){
	   if(myid == 0)
	     cout << "Erz map file : " << file_mapEz << " doesn't not exist" << endl; 
	   exit(1);
	 }
	 file_map.close();	
	 file_map.open(file_mapB.c_str(),ios::in);
	 if(!file_map){
	   if(myid == 0)
	     cout << "B map file : " << file_mapB << " doesn't not exist" << endl; 
	   exit(1);
	 }
	 file_map.close();						
       }

       if(n_file_map==3){
	 file_mapEr = sparser.myrealtrap.filenames[0];
	 file_mapEz = sparser.myrealtrap.filenames[1];
	 file_mapB = sparser.myrealtrap.filenames[2];
	 file_map.open(file_mapEr.c_str(),ios::in);
	 if(!file_map){
	   if(myid == 0)
	     cout << "Er map file : " << file_mapEr << "doesn't not exist" << endl; 
	   exit(1);
	 }
	 file_map.close();	
	 file_map.open(file_mapEz.c_str(),ios::in);
	 if(!file_map){
	   if(myid == 0)
	     cout << "Ez map file : " << file_mapEz << "doesn't not exist" << endl; 
	   exit(1);
	 }
	 file_map.close();	
	 file_map.open(file_mapB.c_str(),ios::in);
	 if(!file_map){
	   if(myid == 0)
	     cout << "B map file : " << file_mapB << "doesn't not exist" << endl; 
	   exit(1);
	 }
	 file_map.close();
       }
       if((n_file_map!=2)&(n_file_map!=3)){
	 if(myid == 0)
	   cout << "Error : Wrong number of EM map files" << endl; 
	 exit(1);
       }
     }

    //----------------------------------------------
    // CALC TRAP
     if(sparser.mycalctrap.flag){
       odev.forcev.trap_param.trap_config = 2;
       parameter1 = sparser.mycalctrap.param1;
       parameter2 = sparser.mycalctrap.param2;
       InitFieldCalc(parameter1,parameter2);
       CalcElectricField();
       CalcMagneticField();
     }
    //----------------------------------------------
    // POT TRAP
    if(sparser.mypottrap.flag){
        odev.forcev.trap_param.trap_config = 4;
        odev.forcev.Potential_map.ReadEPotential(sparser.mypottrap.potmap_filename);
        odev.forcev.trap_param.B0 = sparser.mypottrap.Bz;
        NoIdealTrap();
    }
    
     // OUTPUT
     if(myid==0)
       {
	 cout << "   Simulation Parameters:" << endl;
	 if(sparser.mycreateparticles.flag){
		cout << NSPARTICLES << " particles created"<< endl;
	   for(unsigned int i=0; i<Ions_particles.size();i++){
	       Ion_tmp=Ions_particles[i];
	       cout << "\t"<<Ion_tmp.GetName() << " " << Ion_tmp.Getmass()/amu << " ";
	       if (Ion_tmp.Getcharge() > 0) cout<<Ion_tmp.Getcharge()<<"+\n";
		 else cout<<Ion_tmp.Getcharge()*-1<<"-\n";
	   }
	}

     if(sparser.mycreatecloud.flag){
	cout << NRPARTICLES << " particles for cloud"<< endl;
       if(myid==0)
	 vector<int> nfractions;
        //Print out ion name, charge and fraction
	int index=0;
       for(unsigned int i=0; i<sparser.mycloudparts.cloudfracs.size();i++){
	   index +=sparser.mycloudparts.cloudfracs[i].second/100.*Ions_cloud.size();
	   Ion_tmp=Ions_cloud[index-1];
	   cout << "\t"<<Ion_tmp.GetName() << " " << Ion_tmp.Getmass()/amu << " ";
	   if (Ion_tmp.Getcharge() > 0) cout<<Ion_tmp.Getcharge()<<"+";
		 else cout<<Ion_tmp.Getcharge()*-1<<"-";
	  cout<< " " << sparser.mycloudparts.cloudfracs[i].second << " %" << endl;
	}
       
	 cout << "Cloud semiaxis: " << semiaxis_cloud[0] << " " << semiaxis_cloud[1] << " " << semiaxis_cloud[2] << " m" << endl;
	 cout << "Cloud center position: " << offset_cloud[0] << " " << offset_cloud[1] << " " << offset_cloud[2] << " m" << endl;
	 cout << "MB energy distribution with <E>= " << energy << " eV" << endl;
       }//if(myid==0)
     }

     if(sparser.myimportdata.flag)
       if(myid==0)
	 cout << "DATA imported from: " << filename_prefix_importdata << endl;
     if(myid==0)
       {       
	 if( int_buffer_gas == 1) {
	   cout << "with buffer gas, pressure: " << p_buff_mbar << " mbar" << endl;}
	 else{
	   cout << "without buffer gas"<< endl;}
	 
	 if(sparser.myidealtrap.flag)
	   {
	     cout << "IDEAL TRAP configuration, parameters in trapparameters.h" << endl;
	   }
	 if(sparser.myrealtrap.flag)
	   {
	     cout << "REAL TRAP with EM field maps :"  << endl;
	     if(n_file_map==2)
	       cout <<"\t " << file_mapEz << " , " <<file_mapB << endl;
	     if(n_file_map==3)
	       cout <<"\t " << file_mapEr << " , " << file_mapEz << " , " <<file_mapB << endl;
	   }
	 if(sparser.mycalctrap.flag)
	   {
	     cout << "Trap electric and magnetic field calculated with equations defined in fieldcalc.cpp" << endl;
	   }
	 cout << "ODE order : " << ODEORDER << ", timestep : " << timestep <<  " s " << endl;
	 if(adaptative_timestep) 
	   cout << "with adaptative time step" << endl;
	 if(Coulomb==1)
	   cout << "with Coulomb force calculation, Coulomb scale factor : " << Coulomb_scale_factor << endl;
	 else
	   cout << "without Coulomb force calculation" << endl;
	 cout << "filename prefix: " << filename_prefix_ << endl;
	 cout << "data saved every " << printinterval_timestep << " s" << endl;
	 
       }    
	
//initialize the prefix of the filenames
	size_t found=filename_prefix_.find_last_of("/\\");
	string folder = "./";	
	if(filename_prefix_.substr(0,found).length() != filename_prefix_.length()){  	
		folder = filename_prefix_.substr(0,found);
		cout << " folder: " << folder << endl;
#ifdef __linux__
	    	mkpath(folder,0755);
#endif // __linux__
 	}



#ifdef __MPI_ON__
       char filename_prefix[250];
       stringstream ssltemp2;ssltemp2<< filename_prefix_ << "-" << myid;
       //stringstream ssltemp2;ssltemp2<<"test" << myid;
       ssltemp2>>filename_prefix;
#else // __MPI_ON__  
       const  char * filename_prefix = filename_prefix_.c_str();
       //stringstream ssltemp2;
       //ssltemp2<< filename_prefix_;
       //ssltemp2>>filename_prefix;
#endif // __MPI_ON__   
	
//initialize the random seed
       unsigned long seed;	
       seed = time(0);

       // open logger 
       stringstream ssltemp;ssltemp<<filename_prefix<<"_logfile.txt";
       char filename_logger[250];ssltemp>>filename_logger;
       tmplogger.open(filename_logger);
	
#ifdef __NBODY_ON__	
 	tmplogger << "Using NBODY code \n";
#endif
#ifdef __CUNBODY_ON__
 	tmplogger <<"Using CUNBODY code\n";
#endif // __NBODY_ON__	

 	//Print out the particle information in
 	//true: each particle has his own file (default)
 	//false: print out all information in one .._pAll.txt
 	use_particle_file(seperateparticlefile_bool,_cloud);
	
	
 	//Initialize the ODE. First the filename_prefix, then the Order of the ODE can be 1,4 or 5, next is the initial timestep (should be around 1e-9), next is a boulean
 	//true= with adaptive stepsize, false is without adaptive stepsize.
 	InitIonFly(filename_prefix,ODEORDER,timestep,adaptative_stepsize,_cloud,odev);
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
 	if(sparser.myrealtrap.flag)
	  {
	    if(n_file_map==2)
	      NoIdealTrap(&file_mapEz[0],&file_mapB[0]);
	    if(n_file_map==3)
	      NoIdealTrap(&file_mapEr[0],&file_mapEz[0],&file_mapB[0]);
	  }
 	if(sparser.mycalctrap.flag)
	  {
	    NoIdealTrap();
	  }

	//Import the data from previous simulations
 	//Or create the ions according to a maxwell boltzmann distribution with the maxima around 1.5eV 
 	//MPI_ImportData("Rb512-dip");	


 	//log the input charges and fractions. CHANGE XXX
 	//for(unsigned int i=0; i<ion_name.size();i++){
	//  tmplogger<<Ions_cloud[i].GetName();tmplogger<<" ";tmplogger<<Ions_cloud[i].Getmass()/amu;tmplogger<<" ";			
	//  if (Ions_cloud[i].Getcharge() > 0) {tmplogger<<Ions_cloud[i].Getcharge();tmplogger<<"+";}
	//  else {tmplogger<<Ions_cloud[i].Getcharge()*-1;tmplogger<<"-";}
	//  tmplogger<< " ";tmplogger<<fraction_of_ion[i] ;tmplogger<<" %\n";
 	//}

	if(sparser.mycreateparticles.flag){
	  CreateCloud(NSPARTICLES,Ions_particles,_cloud,odev.forcev.trap_param);
 	}

	if(sparser.mycreatecloud.flag){
	  InitCloud(NRPARTICLES,energy,seed,semiaxis_cloud,offset_cloud, Ions_cloud,_cloud,odev.forcev.trap_param);      
	}

	if(sparser.myimportdata.flag){
	  if(sparser.mycreateparticles.flag && !sparser.mycreatecloud.flag){
	    npart_test = ImportData(filename_prefix_importdata.c_str(),Table,pdgTable, _cloud,odev.forcev.trap_param);
	    NRPARTICLES = npart_test;
	    CreateCloud(NSPARTICLES,Ions_particles,_cloud,odev.forcev.trap_param);
	  }else if(!sparser.mycreateparticles.flag && sparser.mycreatecloud.flag){
	    npart_test = ImportData(filename_prefix_importdata.c_str(),Table,pdgTable, _cloud,odev.forcev.trap_param);
	    InitCloud(NRPARTICLES,energy,seed,semiaxis_cloud,offset_cloud, Ions_cloud,_cloud,odev.forcev.trap_param);            

	  }else{

	  }

	}


#ifdef __MPI_ON__ 
 	MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__

    //----------------------------------------------
    // OPERATION
	
	if(sparser.myne.flag){
	  time_operation = sparser.myne.time;

	  Ope_tmp.SetName("NE");
	  Ope_tmp.SetTime(time_operation);
	  Ope_tmp.SetBuff(p_buff_mbar);
	  if(int_buffer_gas == 1){Ope_tmp.SetBuffBool(true);}
	  else{ Ope_tmp.SetBuffBool(false);}
	  Ope.AddOperation(Ope_tmp);	
            
	}

	if(sparser.mydew.flag || sparser.mydeb.flag || sparser.myqew.flag || sparser.myqeb.flag || sparser.myoew.flag || sparser.myoeb.flag){
	  // DEW
	  if(sparser.mydew.flag){
	    time_operation = sparser.mydew.time;
	    tmp_string = sparser.mydew.EigenLett;
	    if(tmp_string == "f"){
	      tmp_frequency = sparser.mydew.freq_bias;
	      tmp_amplitude = sparser.mydew.amp;
	    }else{
	      tmp_name = sparser.mydew.Element;
	      tmp_f_bias = sparser.mydew.freq_bias;
	      tmp_amplitude = sparser.mydew.amp;
	      if(Table.TestIonName(tmp_name)){
		tmp_double = Table.ResearchMass(tmp_name);
	      }else if(pdgTable.TestPDGName(tmp_name)){
		int pdgid = pdgTable.GetPDGId(orig_name);
		tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	      }else{
		if(myid == 0)cout << "Error: can't calculate particle mass for DEW!" << endl;
		exit(-1);
	      }
	      Ope_tmp.SetIonName(tmp_name);
	      Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	      if(tmp_string == "c")tmp_frequency = Ion_tmp.Getwc()+tmp_f_bias;
	      if(tmp_string == "m")tmp_frequency = Ion_tmp.Getwmin()+tmp_f_bias;
	      if(tmp_string == "p")tmp_frequency = Ion_tmp.Getwplus()+tmp_f_bias;
	      if(myid==0)
		{
		  //cout << tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude/(2.*(Ion_tmp.Getwplus()-Ion_tmp.Getwmin())) << endl;
		}
	      }
	      if(tmp_string!="c"&&tmp_string!="f"&&tmp_string!="m"&&tmp_string!="p"){
		if(myid==0)
		  cout << "Error in operation label" << endl;
		return;
	      }

	      Ope_tmp.SetName("DE");
	      Ope_tmp.SetTime(time_operation);		    
	      Ope_tmp.SetAmplitude(tmp_amplitude);			
	      Ope_tmp.SetFrequency(tmp_frequency);
	      Ope_tmp.SetFrequencyBias(tmp_f_bias);
	      Ope_tmp.SetBuffBool(false);
	      Ope_tmp.SetBuff(0.);
	      Ope.AddOperation(Ope_tmp);
	  }

	  // DEB
	  if(sparser.mydeb.flag){
	    time_operation = sparser.mydeb.time;
	    tmp_string = sparser.mydeb.EigenLett;
	    if(tmp_string == "f"){
	      tmp_frequency = sparser.mydeb.freq_bias;
	      tmp_amplitude = sparser.mydeb.amp;
	    }else{
	      tmp_name = sparser.mydeb.Element;
	      tmp_f_bias = sparser.mydeb.freq_bias;
	      tmp_amplitude = sparser.mydeb.amp;
	      if(Table.TestIonName(tmp_name)){
		tmp_double = Table.ResearchMass(tmp_name);
	      }else if(pdgTable.TestPDGName(tmp_name)){
		int pdgid = pdgTable.GetPDGId(orig_name);
		tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	      }else{
		if(myid == 0)cout << "Error: can't calculate particle mass for DEB!" << endl;
		exit(-1);
	      }
	      Ope_tmp.SetIonName(tmp_name);
	      Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	      if(tmp_string == "c")tmp_frequency = Ion_tmp.Getwc()+tmp_f_bias;
	      if(tmp_string == "m")tmp_frequency = Ion_tmp.Getwmin()+tmp_f_bias;
	      if(tmp_string == "p")tmp_frequency = Ion_tmp.Getwplus()+tmp_f_bias;
	      if(myid==0)
		{
		  //cout << tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude/(2.*(Ion_tmp.Getwplus()-Ion_tmp.Getwmin())) << endl;
		}
	      }
	      if(tmp_string!="c"&&tmp_string!="f"&&tmp_string!="m"&&tmp_string!="p"){
		if(myid==0)
		  if(myid == 0)cout << "Error in operation label" << endl;
		return;
	      }

	      Ope_tmp.SetName("DE");
	      Ope_tmp.SetTime(time_operation);		    
	      Ope_tmp.SetAmplitude(tmp_amplitude);			
	      Ope_tmp.SetFrequency(tmp_frequency);
	      Ope_tmp.SetFrequencyBias(tmp_f_bias);
	      Ope_tmp.SetBuffBool(true);
	      Ope_tmp.SetBuff(p_buff_mbar);
	      Ope.AddOperation(Ope_tmp);
	  }

	  // QEW
	  if(sparser.myqew.flag){
	    time_operation = sparser.myqew.time;
	    tmp_string = sparser.myqew.EigenLett;
	    if(tmp_string == "f"){
	      tmp_frequency = sparser.myqew.freq_bias;
	      tmp_amplitude = sparser.myqew.amp;
	    }else{
	      tmp_name = sparser.myqew.Element;
	      tmp_f_bias = sparser.myqew.freq_bias;
	      tmp_amplitude = sparser.myqew.amp;
	      if(Table.TestIonName(tmp_name)){
		tmp_double = Table.ResearchMass(tmp_name);
	      }else if(pdgTable.TestPDGName(tmp_name)){
		int pdgid = pdgTable.GetPDGId(orig_name);
		tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	      }else{
		if(myid == 0)cout << "Error: can't calculate particle mass for QEW!" << endl;
		exit(-1);
	      }
	      Ope_tmp.SetIonName(tmp_name);
	      Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	      if(tmp_string == "c")tmp_frequency = Ion_tmp.Getwc()+tmp_f_bias;
	      if(tmp_string == "m")tmp_frequency = Ion_tmp.Getwmin()+tmp_f_bias;
	      if(tmp_string == "p")tmp_frequency = Ion_tmp.Getwplus()+tmp_f_bias;
	      if(myid==0)
		{
		  //cout << tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude/(2.*(Ion_tmp.Getwplus()-Ion_tmp.Getwmin())) << endl;
		}
	      }
	      if(tmp_string!="c"&&tmp_string!="f"&&tmp_string!="m"&&tmp_string!="p"){
		if(myid==0)
		  cout << "Error in operation label" << endl;
		return;
	      }

	      Ope_tmp.SetName("QE");
	      Ope_tmp.SetTime(time_operation);		    
	      Ope_tmp.SetAmplitude(tmp_amplitude);			
	      Ope_tmp.SetFrequency(tmp_frequency);
	      Ope_tmp.SetFrequencyBias(tmp_f_bias);
	      Ope_tmp.SetBuffBool(false);
	      Ope_tmp.SetBuff(0.);
	      Ope.AddOperation(Ope_tmp);
	  }

	  // QEB
	  if(sparser.myqeb.flag){
	    time_operation = sparser.myqeb.time;
	    tmp_string = sparser.myqeb.EigenLett;
	    if(tmp_string == "f"){
	      tmp_frequency = sparser.myqeb.freq_bias;
	      tmp_amplitude = sparser.myqeb.amp;
	    }else{
	      tmp_name = sparser.myqeb.Element;
	      tmp_f_bias = sparser.myqeb.freq_bias;
	      tmp_amplitude = sparser.myqeb.amp;
	      if(Table.TestIonName(tmp_name)){
		tmp_double = Table.ResearchMass(tmp_name);
	      }else if(pdgTable.TestPDGName(tmp_name)){
		int pdgid = pdgTable.GetPDGId(orig_name);
		tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	      }else{
		if(myid == 0)cout << "Error: can't calculate particle mass for QEB!" << endl;
		exit(-1);
	      }
	      Ope_tmp.SetIonName(tmp_name);
	      Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	      if(tmp_string == "c")tmp_frequency = Ion_tmp.Getwc()+tmp_f_bias;
	      if(tmp_string == "m")tmp_frequency = Ion_tmp.Getwmin()+tmp_f_bias;
	      if(tmp_string == "p")tmp_frequency = Ion_tmp.Getwplus()+tmp_f_bias;
	      if(myid==0)
		{
		  //cout << tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude/(2.*(Ion_tmp.Getwplus()-Ion_tmp.Getwmin())) << endl;
		}
	      }
	      if(tmp_string!="c"&&tmp_string!="f"&&tmp_string!="m"&&tmp_string!="p"){
		if(myid==0)
		  cout << "Error in operation label" << endl;
		return;
	      }

	      Ope_tmp.SetName("QE");
	      Ope_tmp.SetTime(time_operation);		    
	      Ope_tmp.SetAmplitude(tmp_amplitude);			
	      Ope_tmp.SetFrequency(tmp_frequency);
	      Ope_tmp.SetFrequencyBias(tmp_f_bias);
	      Ope_tmp.SetBuffBool(true);
	      Ope_tmp.SetBuff(p_buff_mbar);
	      Ope.AddOperation(Ope_tmp);
	  }

	  // OEW
	  if(sparser.myoew.flag){
	    time_operation = sparser.myoew.time;
	    tmp_string = sparser.myoew.EigenLett;
	    if(tmp_string == "f"){
	      tmp_frequency = sparser.myoew.freq_bias;
	      tmp_amplitude = sparser.myoew.amp;
	    }else{
	      tmp_name = sparser.myoew.Element;
	      tmp_f_bias = sparser.myoew.freq_bias;
	      tmp_amplitude = sparser.myoew.amp;
	      if(Table.TestIonName(tmp_name)){
		tmp_double = Table.ResearchMass(tmp_name);
	      }else if(pdgTable.TestPDGName(tmp_name)){
		int pdgid = pdgTable.GetPDGId(orig_name);
		tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	      }else{
		if(myid == 0)cout << "Error: can't calculate particle mass for OEW!" << endl;
		exit(-1);
	      }
	      Ope_tmp.SetIonName(tmp_name);
	      Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	      if(tmp_string == "c")tmp_frequency = Ion_tmp.Getwc()+tmp_f_bias;
	      if(tmp_string == "m")tmp_frequency = Ion_tmp.Getwmin()+tmp_f_bias;
	      if(tmp_string == "p")tmp_frequency = Ion_tmp.Getwplus()+tmp_f_bias;
	      if(myid==0)
		{
		  //cout << tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude/(2.*(Ion_tmp.Getwplus()-Ion_tmp.Getwmin())) << endl;
		}
	      }
	      if(tmp_string!="c"&&tmp_string!="f"&&tmp_string!="m"&&tmp_string!="p"){
		if(myid==0)
		  cout << "Error in operation label" << endl;
		return;
	      }

	      Ope_tmp.SetName("OE");
	      Ope_tmp.SetTime(time_operation);		    
	      Ope_tmp.SetAmplitude(tmp_amplitude);			
	      Ope_tmp.SetFrequency(tmp_frequency);
	      Ope_tmp.SetFrequencyBias(tmp_f_bias);
	      Ope_tmp.SetBuffBool(false);
	      Ope_tmp.SetBuff(0.);
	      Ope.AddOperation(Ope_tmp);
	  }

	  // OEB
	  if(sparser.myoeb.flag){
	    time_operation = sparser.myoeb.time;
	    tmp_string = sparser.myoeb.EigenLett;
	    if(tmp_string == "f"){
	      tmp_frequency = sparser.myoeb.freq_bias;
	      tmp_amplitude = sparser.myoeb.amp;
	    }else{
	      tmp_name = sparser.myoeb.Element;
	      tmp_f_bias = sparser.myoeb.freq_bias;
	      tmp_amplitude = sparser.myoeb.amp;
	      if(Table.TestIonName(tmp_name)){
		tmp_double = Table.ResearchMass(tmp_name);
	      }else if(pdgTable.TestPDGName(tmp_name)){
		int pdgid = pdgTable.GetPDGId(orig_name);
		tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	      }else{
		if(myid == 0)cout << "Error: can't calculate particle mass for OEB!" << endl;
		exit(-1);
	      }
	      Ope_tmp.SetIonName(tmp_name);
	      Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	      if(tmp_string == "c")tmp_frequency = Ion_tmp.Getwc()+tmp_f_bias;
	      if(tmp_string == "m")tmp_frequency = Ion_tmp.Getwmin()+tmp_f_bias;
	      if(tmp_string == "p")tmp_frequency = Ion_tmp.Getwplus()+tmp_f_bias;
	      if(myid==0)
		{
		  //cout << tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude << endl;
		  //cout << Ion_tmp.Getkq_div_u()*tmp_amplitude/(2.*(Ion_tmp.Getwplus()-Ion_tmp.Getwmin())) << endl;
		}
	      }
	      if(tmp_string!="c"&&tmp_string!="f"&&tmp_string!="m"&&tmp_string!="p"){
		if(myid==0)
		  cout << "Error in operation label" << endl;
		return;
	      }

	      Ope_tmp.SetName("OE");
	      Ope_tmp.SetTime(time_operation);		    
	      Ope_tmp.SetAmplitude(tmp_amplitude);			
	      Ope_tmp.SetFrequency(tmp_frequency);
	      Ope_tmp.SetFrequencyBias(tmp_f_bias);
	      Ope_tmp.SetBuffBool(true);
	      Ope_tmp.SetBuff(p_buff_mbar);
	      Ope.AddOperation(Ope_tmp);
	  }	
	} // end if DEW, DEB, QEW, QEB, OEW, OEB

	if(sparser.myrww.flag || sparser.myrwb.flag || sparser.myaww.flag || sparser.myawb.flag){
	  // RWW
	  if(sparser.myrww.flag){
	    time_operation = sparser.myrww.time;
	    tmp_order = sparser.myrww.order;
	    tmp_frequency = sparser.myrww.frequency;
	    tmp_amplitude = sparser.myrww.amp;

	    Ope_tmp.SetName("RW");
	    Ope_tmp.SetTime(time_operation);
	    Ope_tmp.SetOrder(tmp_order);		    
	    Ope_tmp.SetAmplitude(tmp_amplitude);			
	    Ope_tmp.SetFrequency(tmp_frequency);
	    Ope_tmp.SetBuffBool(false);
	    Ope_tmp.SetBuff(0.);
	    Ope.AddOperation(Ope_tmp);	
	  }	
	  // RWB
	  if(sparser.myrwb.flag){
	    time_operation = sparser.myrwb.time;
	    tmp_order = sparser.myrwb.order;
	    tmp_frequency = sparser.myrwb.frequency;
	    tmp_amplitude = sparser.myrwb.amp;

	    Ope_tmp.SetName("RW");
	    Ope_tmp.SetTime(time_operation);
	    Ope_tmp.SetOrder(tmp_order);		    
	    Ope_tmp.SetAmplitude(tmp_amplitude);			
	    Ope_tmp.SetFrequency(tmp_frequency);
	    Ope_tmp.SetBuffBool(true);
	    Ope_tmp.SetBuff(p_buff_mbar);
	    Ope.AddOperation(Ope_tmp);	
	  }
	  // AWW
	  if(sparser.myaww.flag){
	    time_operation = sparser.myaww.time;
	    tmp_order = sparser.myaww.order;
	    tmp_frequency = sparser.myaww.frequency;
	    tmp_amplitude = sparser.myaww.amp;

	    Ope_tmp.SetName("AW");
	    Ope_tmp.SetTime(time_operation);
	    Ope_tmp.SetOrder(tmp_order);		    
	    Ope_tmp.SetAmplitude(tmp_amplitude);			
	    Ope_tmp.SetFrequency(tmp_frequency);
	    Ope_tmp.SetBuffBool(false);
	    Ope_tmp.SetBuff(0.);
	    Ope.AddOperation(Ope_tmp);	
	  }	
	  // AWB
	  if(sparser.myawb.flag){
	    time_operation = sparser.myawb.time;
	    tmp_order = sparser.myawb.order;
	    tmp_frequency = sparser.myawb.frequency;
	    tmp_amplitude = sparser.myawb.amp;

	    Ope_tmp.SetName("AW");
	    Ope_tmp.SetTime(time_operation);
	    Ope_tmp.SetOrder(tmp_order);		    
	    Ope_tmp.SetAmplitude(tmp_amplitude);			
	    Ope_tmp.SetFrequency(tmp_frequency);
	    Ope_tmp.SetBuffBool(true);
	    Ope_tmp.SetBuff(p_buff_mbar);
	    Ope.AddOperation(Ope_tmp);	
	  }
	}

	//AC
	if(sparser.myac.flag){
	  time_operation = sparser.myac.time;
	  tmp_order = sparser.myac.order;
	  tmp_string = sparser.myac.Element;
	  tmp_amplitude = sparser.myac.amp;
	  if(Table.TestIonName(tmp_string)){
	    tmp_double = Table.ResearchMass(tmp_string);
	  }else if(pdgTable.TestPDGName(tmp_string)){
	    int pdgid = pdgTable.GetPDGId(tmp_string);
	    tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	  }else{
	    if(myid == 0)cout << "Error: can't calculate particle mass for AC!" << endl;
	    exit(-1);
	  }
	  Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	  Ope_tmp.SetName("AC");
	  Ope_tmp.SetTime(time_operation);
	  //Ope_tmp.SetOrder(tmp_order);
	  Ope_tmp.SetAmplitude(tmp_amplitude);
	  //Ope_tmp.SetFrequency(tmp_f_bias);
	  if(tmp_order==0)
	    Ope_tmp.SetFrequency(Ion_tmp.Getwplus()+sqrt(Ion_tmp.Getwz2()));
	  if(tmp_order==1)
	    Ope_tmp.SetFrequency(Ion_tmp.Getwplus()-sqrt(Ion_tmp.Getwz2()));
	  if(tmp_order==2)
	    Ope_tmp.SetFrequency(Ion_tmp.Getwmin()-sqrt(Ion_tmp.Getwz2()));
	  if(tmp_order==3)
	    Ope_tmp.SetFrequency(Ion_tmp.Getwmin()+sqrt(Ion_tmp.Getwz2()));
	  Ope_tmp.SetBuff(0.);Ope_tmp.SetBuffBool(false);
	  Ope.AddOperation(Ope_tmp);
	}  
	//SC
	if(sparser.mysc.flag){
	  time_operation = sparser.mysc.time;
          tmp_order = sparser.mysc.order;
	  if(tmp_order==0) // frequencies and amplitudes are given
	    {
	      tmp_frequency = sparser.mysc.freq1;
	      tmp_amplitude = sparser.mysc.amp1;
	      tmp_frequency2 = sparser.mysc.freq2;
	      tmp_amplitude2 = sparser.mysc.amp2;
	    }
	  else // name of the isotope is given
	    {
	      tmp_name = sparser.mysc.Element;
	      tmp_f_bias = sparser.mysc.freq1;
	      tmp_amplitude = sparser.mysc.amp1;
	      tmp_f_bias2 = sparser.mysc.freq2;
	      tmp_amplitude2 = sparser.mysc.amp2;
	      if(Table.TestIonName(tmp_name)){
		tmp_double = Table.ResearchMass(tmp_name);
	      }else if(pdgTable.TestPDGName(tmp_name)){
		int pdgid = pdgTable.GetPDGId(tmp_name);
		tmp_double = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
	      }else{
		if(myid == 0)cout << "Error: can't calculate particle mass for SC!" << endl;
		exit(-1);
	      }
	      Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
	      tmp_frequency = Ion_tmp.Getwmin() + tmp_f_bias;
	      tmp_frequency2 = Ion_tmp.Getwc() + tmp_f_bias2;
	    }
            
             Ope_tmp.SetName("SC");
             Ope_tmp.SetTime(time_operation);
             Ope_tmp.SetBuff(0.);Ope_tmp.SetBuffBool(false);
             Ope_tmp.SetAmplitude(tmp_amplitude);
             Ope_tmp.SetFrequency(tmp_frequency);
             Ope_tmp.SetAmplitude2(tmp_amplitude2);
             Ope_tmp.SetFrequency2(tmp_frequency2);
             Ope.AddOperation(Ope_tmp);
	}
      //AR
      if(sparser.myar.flag){
	time_operation = sparser.myar.time;
       	tmp_frequency = sparser.myar.freq1;
      	tmp_amplitude = sparser.myar.amp1;
      	tmp_frequency2 = sparser.myar.freq2;
        Ope_tmp.SetName("AR");
        Ope_tmp.SetTime(time_operation);
        Ope_tmp.SetBuff(0.);Ope_tmp.SetBuffBool(false);
        Ope_tmp.SetAmplitude(tmp_amplitude);
        Ope_tmp.SetFrequency(tmp_frequency);
        Ope_tmp.SetFrequency2(tmp_frequency2);
        Ope.AddOperation(Ope_tmp);
     }



//         if(subline=="AQ")
//         {
// 	        gui_file >> time_operation;
//             gui_file >> tmp_order;
//             if(tmp_order==0) // frequencies and amplitudes are given
//             {
//                 gui_file >> tmp_frequency;
//                 gui_file >> tmp_amplitude;
//                 gui_file >> tmp_frequency2;
//                 gui_file >> tmp_amplitude2;
//             }
//             else // name of the isotope is given
//             {
                
//                 gui_file >> tmp_name;
//                 gui_file >> tmp_f_bias;
//                 gui_file >> tmp_amplitude;
//                 gui_file >> tmp_f_bias2;
//                 gui_file >> tmp_amplitude2;
//                 gui_file >> tmp_f_bias3;
//                 gui_file >> tmp_amplitude3;
//                 gui_file >> tmp_f_bias4;
//                 gui_file >> tmp_amplitude4;
//                 tmp_double = Table.ResearchMass(tmp_name); // mass of the element
//                 Ion_tmp.SetParameters(tmp_double,100.,odev.forcev.trap_param);
//                 tmp_frequency2 =  Ion_tmp.Getwc()  + tmp_f_bias;
//                 tmp_frequency3 = Ion_tmp.Getwplus()-sqrt(Ion_tmp.Getwz2()) + tmp_f_bias2;
//                 tmp_frequency = Ion_tmp.Getwmin();
//                 tmp_frequency4 = sqrt(Ion_tmp.Getwz2());

//             }
            
//             Ope_tmp.SetName(subline);
//             Ope_tmp.SetTime(time_operation);
//             Ope_tmp.SetBuff(0.);Ope_tmp.SetBuffBool(false);
//             Ope_tmp.SetAmplitude(tmp_amplitude);
//             Ope_tmp.SetFrequency(tmp_frequency);
//             Ope_tmp.SetAmplitude2(tmp_amplitude2);
//             Ope_tmp.SetFrequency2(tmp_frequency2);
//             Ope_tmp.SetAmplitude3(tmp_amplitude3);
//             Ope_tmp.SetFrequency3(tmp_frequency3);
//             Ope_tmp.SetAmplitude4(tmp_amplitude4);
//             Ope_tmp.SetFrequency4(tmp_frequency4);
//             Ope.AddOperation(Ope_tmp);
//         }
	
 	if(myid==0)
 	{
 		cout << "OPERATIONS" << endl;
 		Ope.Write();
 		cout << "*********************************************" << endl;
 	}
			
#ifdef __MPI_ON__ 
  	MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__		
 	Ope.Launch(_cloud,odev);
	
 	ExitIonFly(_cloud);

 	/* PostProces the data */
 	if(seperateparticlefile_bool)
 	  BundleData(filename_prefix,NRPARTICLES+NSPARTICLES);
 	//email_it("G100_bundled.txt" , "blabla@blablab.com");
 	//see http://cc.byexamples.com/20070120/print-color-string-without-ncurses/
	if(myid==0)
	  printf("%c[%d;%dmProgram Ended%c[%dm\n",27,4,31,27,0);   
	
	/* clear GPU memory */
    
#ifdef __MPI_ON__ 
 	MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__		 
   

	
	
	
 	/***************  finalize  MPI ****************/
#ifdef __MPI_ON__
  	MPI::COMM_WORLD.Barrier();		
 	// Concatenate data	
 	if(!seperateparticlefile_bool)
 	  MPI_Concatenate_Data(NRPARTICLES+NSPARTICLES,filename_prefix);
				
  	MPI::COMM_WORLD.Barrier();	
#endif // __MPI_ON__



    return;
}

int closest_sup_power_of_two(int n)
{
    int p=1;
    while(p<n)
    {
        p*=2;
    }
    return p;
}



