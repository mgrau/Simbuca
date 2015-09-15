#include "MPI_simbuca.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator> 
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include "globals.h"
#include "trapparameters.h"
#include "ion_fly.h"
#include "ion.h"
#include "mtrand.h"
#include "logfile.h"
#include "inireader.h"
#include "parser.h"
#include "IonTable.h"
#include "PDGTable.h"
#include "operation.h"
#include "initialization.h"

#ifdef __linux
#include <sys/stat.h>
#include "errno.h"
#endif

#define MAX_OPERATIONS 32
#define is_power_of_two(x) ( ((x&(x-1)) == 0) && (x!=0) )
#define is_in(sections,section) (find(sections.begin(),sections.end(),section) != sections.end())

LogFile tmplogger;

MTRand algrand;
IonTable Table("../libraries/ame/mass.mas12");
PDGTable pdgTable("../libraries/pdg/pdg_table.txt");

vector< double > x_;
vector< double > y_;
vector< double > z_;
vector< double > vx_;
vector< double > vy_;
vector< double > vz_;

#ifdef __linux__
void mkpath(std::string s,mode_t mode) {
    size_t pre=0,pos;
    std::string dir;
    int mdret;

    if(s[s.size()-1]!='/')
        s+='/';// force trailing / so we can handle everything in loop

    while((pos=s.find_first_of('/',pre))!=std::string::npos){
        dir=s.substr(0,pos++);
        pre=pos;
        if(dir.size()==0) continue; // if leading / first time is 0 length
        if((mdret=mkdir(dir.c_str(),mode)) && errno!=EEXIST)
            cout<<"error?"; 
    }
}
#endif

void CreateCloud(int nparticles, vector<Ion> Ions,IonCloud &_cloud, _trap_param & trap_param) {
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
    for(int i=0;i<numprocs;i++) {
        counts[i] = nparticles/numprocs;
        if(i<last_part)
            counts[i] += 1;
    }

    for(int i=1;i<numprocs;i++) {
        displ[i] = displ[i-1] +counts[i-1];
    }

    for(int i=0;i < counts[myid] ;i++) {
        ii = i + displ[myid];
        AddParticle(x_[ii],y_[ii],z_[ii],vx_[ii],vy_[ii],vz_[ii],Ions[ii],_cloud);
    }

    delete [] counts;
    delete [] displ;

    _cloud.UpdateIDs(nparticles);

#else // __MPI_ON__
    frac_ = 0;
    j =0;
    for(int i=0;i<nparticles;i++)
    {
        // cout<<"particle added\t"<<i<<"\t"<<x_[i] <<" "<< y_[i] <<" "<< z_[i] <<" "<< vx_[i] <<" "<< vy_[i] <<" "<< vz_[i] <<" "<< endl;//Ions[i]<<endl;
        AddParticle(x_[i],y_[i],z_[i],vx_[i],vy_[i],vz_[i],Ions[i],_cloud);
    }
#endif // __MPI_ON__
}

void InitCloud(int nparticles, vector<double > eV_max_boltz,int seed,double semiaxis[3], double offset[3], vector<Ion> Ions, IonCloud &_cloud, _trap_param & trap_param) {
    // Creation of the ions with a Maxwell Boltzmann distribution with ion positions in an ellipsoid 
    // eV_max_boltz is the Energy in eV where there is the maximum of the Maxwell Boltzmann distibution.
    // note: The total speed 'vmb' is calculated with a Maxwell Boltzmann distribution (more or less. Check
    // Petterson thesis about ISCool RFQ. 

    // the number of particles must be divible by numprocs!
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
        cout << "InitCloud: nparticles="<<nparticles<<endl;
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


            theta_gas = algrand()*2.0*pi;
            phi_gas = acos(algrand()*2.0-1.0);
            //theta_gas = atan(x/y);
            if(eV_max_boltz.size()==2)
            {
                vmb = sqrt(-2.0*(eV_max_boltz[0]*eV_to_Joule)/(mass)*log(1-sqrt(algrand())));
                vx_[i] = vmb*cos(theta_gas)*sin(phi_gas);
                vy_[i] = vmb*sin(theta_gas)*sin(phi_gas);
                vmb = sqrt(-2.0*(eV_max_boltz[1]*eV_to_Joule)/(mass)*log(1-sqrt(algrand())));
                vz_[i] = vmb*cos(phi_gas);
            }
            else
            {
                vmb = sqrt(-2.0*(eV_max_boltz[0]*eV_to_Joule)/(mass)*log(1-sqrt(algrand())));
                vx_[i] = vmb*cos(theta_gas)*sin(phi_gas);
                vy_[i] = vmb*sin(theta_gas)*sin(phi_gas);
                vz_[i] = vmb*cos(phi_gas);
            }
        } // end loop
    } // end if master process 

    //create the cloud
    CreateCloud(nparticles,Ions,_cloud,trap_param);

    tmplogger<<nparticles;tmplogger<<" particles created. Max-Boltz distribution around ";tmplogger<<eV_max_boltz[0];tmplogger<<" eV. \n";
    tmplogger<<"spheroid cloud with semi-axis's = {";tmplogger<<semiaxis[0];tmplogger<<";";tmplogger<<semiaxis[1];tmplogger<<";";tmplogger<<semiaxis[2];tmplogger<<"}";  
    tmplogger<<" and offset = {";tmplogger<<offset[0];tmplogger<<";";tmplogger<<offset[1];tmplogger<<";";tmplogger<<offset[2];tmplogger<<"} ";
    tmplogger<<" \n";
}

// void DoSimulation(SimParser & parser, IonCloud & _cloud,_ode_vars &odev) throw(const char *) {
void DoSimulation(INIReader & parser, IonCloud & _cloud,_ode_vars &odev) throw(const char *) {
    // MPI
    int myid =0;
#ifdef __MPI_ON__
    myid = MPI::COMM_WORLD.Get_rank();
    int numprocs = MPI::COMM_WORLD.Get_size();
#endif    
    ifstream gui_file;
    string comment;

    set<string> sections = parser.GetSections();
    SLogger slogger("Initialization");
    SMsgType type = INFO;
    SLogWriter::Instance()->SetMinType( type );

    try{
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
        vector<double> energy;
        // GAS
        double p_buff_mbar=0;
        // ODE
        int ODEORDER=0;
        double timestep=0;
        int adaptive_timestep=0;
        bool adaptive_stepsize;
        int Coulomb=0;
        bool Coulomb_enable;
        double Coulomb_scale_factor=0;
        // output file
        string filename_prefix_;
        double printinterval_timestep=0;
        unsigned int separateparticlefile=0;
        bool separateparticlefile_bool;
        // trap
        int n_file_map=0;
        int tmp_int=0;
        string trap_config;
        string file_mapEr,file_mapEz,file_mapB;
        ifstream file_map;
        double parameter1,parameter2=0;
        // OPERATION
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
        string tmp_string;
        string orig_name;
        double time_operation=0;
        operation op;
        operations ops;

        if(myid == 0){
            cout << endl;
            cout << "\t*****************" << endl;
            cout << "\t* S I M B U C A *" << endl;
            cout << "\t*****************" << endl;
            cout << "\nInitialization with Simbuca file" << endl;
            cout << "*********************************************" << endl;
        }
        // CLOUD
        if (is_in(sections,"cloud")) {
            NRPARTICLES = parser.GetReal("cloud","NParticles",1024);
            energy.assign(3,parser.GetReal("cloud","temp",1e-3));
            semiaxis_cloud[0] = parser.GetReal("cloud","sx",0.001);
            semiaxis_cloud[1] = parser.GetReal("cloud","sy",0.001);
            semiaxis_cloud[2] = parser.GetReal("cloud","sz",0.001); 
            offset_cloud[0]   = parser.GetReal("cloud","x0",0.0); 
            offset_cloud[1]   = parser.GetReal("cloud","y0",0.0); 
            offset_cloud[2]   = parser.GetReal("cloud","z0",0.0);

            cout << "[cloud]" << endl;
            cout << "    Nparticles = " << NRPARTICLES << endl;
            cout << "    Temp = " << energy[0] << "eV" << endl;
            cout << "    offset = (" << offset_cloud[0] << "," << offset_cloud[1] << "," << offset_cloud[2] << ")" << endl;
            cout << "    size = (" << semiaxis_cloud[0] << "," << semiaxis_cloud[1] << "," << semiaxis_cloud[2] << ")" << endl;
        }

        // CLOUD COMPOSITON
        odev.with_charge = false;
        if (is_in(sections,"composition")) {
            cout << "[composition]" << endl;
            set<string> components = parser.GetFields("composition");
            for (set<string>::iterator component = components.begin();
                    component != components.end(); component++) {
                string name = *component;
                string original_name = name;
                fraction = parser.GetReal("composition",*component,1.0);

                bool species_found = false;
                double charge = 0;
                double mass = 0;
                // char lastchar = name.pop_back();
                char lastchar = *name.rbegin();
                name = name.substr(0,name.length()-1);
                if (lastchar == '+' || lastchar == '-') {
                    size_t last_index = name.find_last_not_of("0123456789");
                    string charge_string = name.substr(last_index + 1);
                    if (charge_string.empty())
                        charge = 1;
                    else
                        charge = atoi(charge_string.c_str());
                    name = name.substr(0,last_index+1);
                    if (lastchar == '-')
                        charge *= -1;
                }
                else {
                    name.push_back(lastchar);
                    charge = 1;
                }

                if (Table.TestIonName(name)) {
                    mass = Table.ResearchMass(name);
                    species_found = true;
                    // cout << "Found " << name << " in Table with mass = " << mass << endl;
                }
                else if (pdgTable.TestPDGName(name)) {
                    int pdg_id = pdgTable.GetPDGId(name);
                    mass = pdgTable.GetPDGMass(pdg_id)*1000.0/mass_Mev; //PDG table mass in GeV, convert to amu;
                    species_found = true;
                    // cout << "Found " << name << " in PDG Table with mass = " << mass << endl;
                }
                else {
                    size_t first_index = name.find_first_not_of("0123456789");
                    mass = atof(name.substr(0,first_index).c_str());
                    species_found = true;
                }

                if (charge!=1)
                    odev.with_charge = true;

                if(fraction>1)
                    if(myid==0)
                        throw("wrong composition of the cloud in fraction");

                if(!species_found)
                    if(myid==0)
                        throw("Don't know what species you mean");

                cout << "    " << name << ": " << fraction;
                cout << ", charge=" << charge << ", mass=" << mass;
                Ion_tmp = Ion(mass,charge);
                Ion_tmp.SetName(original_name);

                // add particles
                int n = fraction*NRPARTICLES;
                if (std::distance(component,components.end())==1)
                    n = NRPARTICLES-Ions_cloud.size();

                if (n > NRPARTICLES-Ions_cloud.size())
                    n = NRPARTICLES-Ions_cloud.size();

                cout << ", adding " << n << " " << Ion_tmp << endl;

                for(int i=0; i < n; i++)
                    Ions_cloud.push_back(Ion_tmp);

            }
        }

        //ODE
        if (is_in(sections,"ode")) {
            cout << "[ode]" << endl;
            string method = parser.Get("ode","method","default");
            if (!method.compare("DP") || !method.compare("Dormand-Prince") || !method.compare("5")) {
                ODEORDER = 5;
                cout << "    Method = Dormand-Prince" << endl;
            }
            else if (!method.compare("RK") || !method.compare("Runge-Kutta") || !method.compare("4")) {
                ODEORDER = 4;
                cout << "    Method = Runge-Kutta" << endl;
            }
            else if (!method.compare("GE") || !method.compare("Gear") || !method.compare("1")) {
                ODEORDER = 1;
                cout << "    Method = Gear" << endl;
            }
            else {
                ODEORDER = 0;
                cout << "    Method = None!" << endl;
            }
            timestep = parser.GetReal("ode","timestep",1e-9);
            cout << "    timestep = " << timestep << "s" << endl;
            adaptive_stepsize = parser.GetBoolean("ode","adaptive",false);
            cout << "    adaptive = " << adaptive_stepsize << endl;
        }
        else
            if (myid == 0)
                throw("no ODE configuration provided.");

        // Coulomb
        if (is_in(sections,"coulomb")) {
            cout << "[coulomb]" << endl;
            Coulomb_enable = parser.GetBoolean("coulomb","enable",false);
            cout << "    enabled = " << Coulomb_enable << endl;
            Coulomb_scale_factor = parser.GetReal("coulomb","scalefactor",1);
            cout << "    scale_factor = " << Coulomb_scale_factor << endl;
        }

        // IDEAL TRAP
        if (is_in(sections,"trap")) {
            cout << "[trap]" << endl;
            string type = parser.Get("trap","type","ideal");
            if (!type.compare("ideal")) {
                odev.forcev.trap_param.trap_type = 0;
                odev.forcev.trap_param.Vrf = parser.GetReal("trap","Vrf",50); 
                odev.forcev.trap_param.Vdc = parser.GetReal("trap","Vdc",5); 
                odev.forcev.trap_param.kappa = parser.GetReal("trap","kappa",0.025); 
                odev.forcev.trap_param.r0 = parser.GetReal("trap","r0",0.04); 
                odev.forcev.trap_param.z0 = parser.GetReal("trap","z0",0.04); 
                odev.forcev.trap_param.freq_rf = parser.GetReal("trap","freq_rf",0.04); 
                cout << "    type = ideal" << endl;
                cout << "    Vrf = " << odev.forcev.trap_param.Vrf << endl;
                cout << "    Vdc = " << odev.forcev.trap_param.Vdc << endl;
                cout << "    kappa = " << odev.forcev.trap_param.kappa << endl;
                cout << "    r0 = " << odev.forcev.trap_param.r0 << endl;
                cout << "    z0 = " << odev.forcev.trap_param.z0 << endl;
                cout << "    freq_rf = " << odev.forcev.trap_param.freq_rf << endl;
            }
        }

        // IMPORTDATA
        if (is_in(sections,"import")) {
            odev.with_charge = false;
            filename_prefix_importdata = parser.Get("import","prefix","out"); 
#ifndef __MPI_ON__
            int numprocs = 1;
#endif
            for(int k=0;k<numprocs;k++)
            {
                filename_importdata << filename_prefix_importdata  << "_pAll.txt";
                file_importdata.open(filename_importdata.str().c_str(),ios::in);
                if(!file_importdata)
                {
                    if(myid==0)
                        throw("Import file doesn`t exist");
                }
                file_importdata.close();
                filename_importdata.str("");
            }
            npart_test = ImportData(filename_prefix_importdata.c_str(),Table,pdgTable, _cloud,odev.forcev.trap_param);
            NRPARTICLES = npart_test;
        }

        //OUTPUT FILE
        if (is_in(sections,"output")) {
            cout << "[output]" << endl;
            filename_prefix_ = parser.Get("output","prefix","out");
            cout << "    prefix = " << filename_prefix_ << endl;
            printinterval_timestep = parser.GetReal("output","timestep",1e-6);
            cout << "    timestep = " << printinterval_timestep << endl;
            separateparticlefile_bool = parser.GetBoolean("output","separate",false);
            cout << "    separate_files = " << separateparticlefile_bool << endl;
            bool print_after_operation = parser.GetBoolean("output","print_after_operation",false);
            cout << "    print_after_operation = " << print_after_operation << endl;
            bool print_at_x = parser.GetBoolean("output","print_at_x",false);
            cout << "    print_at_x = " << print_at_x << endl;
            double print_x = parser.GetReal("output","print_x",false);
            cout << "    print_x = " << print_x << endl;

            odev.print_after_operation = print_after_operation;
            odev.print_at_x = print_at_x;
            odev.print_x = print_x;
        }
        else
            if(myid ==0)
                throw("no Output file configuration provided");

        //OPERATIONS
        for (int i=0; i<MAX_OPERATIONS; i++) {
            ostringstream op_num;
            op_num << "operation" << i;
            if (is_in(sections,op_num.str())) {
                cout << "[" << op_num.str() << "]" << endl;
                string op_type = parser.Get(op_num.str(),"type","none");
                op.time = parser.GetReal(op_num.str(),"duration",1e-3);
                op.print_time = parser.GetReal(op_num.str(),"print_timestep",0.0);
                op.timestep = parser.GetReal(op_num.str(),"ode_timestep",0.0);
                cout << "    type = " << op_type << endl;
                if (op_type == "normal" || op_type=="none" || op_type=="n" || op_type=="trap") {
                    op.name = "normal";
                    ops.add(op);
                }
                else if (op_type == "tof" || op_type=="TOF") {
                    op.name = "tof";
                    odev.forcev.trap_param.E_kick = parser.GetReal(op_num.str(),"Ekick",0.0);
                    ops.add(op);
                }
                else if (op_type == "trap ramp" || op_type == "trap_ramp" ) {
                    op.name = "trap_ramp";
                    odev.forcev.trap_param.newVrf = parser.GetReal(op_num.str(),"Vrf",0.0); 
                    odev.forcev.trap_param.newVdc = parser.GetReal(op_num.str(),"Vdc",0.0); 
                    ops.add(op);
                }
                else if (op_type == "dissociation") {
                    op.name = "dissociation";
                    op.dissociation_fraction = parser.GetReal(op_num.str(),"fraction",1.0);
                    op.dissociation_reactant = parser.Get(op_num.str(),"reactant","none");
                    op.dissociation_product = parser.Get(op_num.str(),"product","none");
                    op.product_mass = parser.GetReal(op_num.str(),"mass",0.0);
                    op.product_energy = parser.GetReal(op_num.str(),"energy",0.0);
                    ops.add(op);
                }
                else if (op_type == "nop" || op_type=="none") {
                    op.name = "nop";
                    ops.add(op);
                }
            }
        }


        if(myid==0) {
            cout << endl << "OPERATIONS" << endl;
            ops.write();
            cout << "*********************************************" << endl;
        }
        // OUTPUT
        if(myid==0) {
            if(true){
                cout << NRPARTICLES << " particles for cloud"<< endl;
                if(myid==0)
                    vector<int> nfractions;
            }
        }

        //initialize the prefix of the filenames
        size_t found=filename_prefix_.find_last_of("/\\");
        string folder = "./";
        if(filename_prefix_.substr(0,found).length() != filename_prefix_.length()){
            folder = filename_prefix_.substr(0,found);
            cout << "outputfolder: " << folder << endl;
#ifdef __linux__
            mkpath(folder,0755);
#endif // __linux__
        }

#ifdef __MPI_ON__
        char filename_prefix[250];
        stringstream ssltemp2;ssltemp2<< filename_prefix_ << "-" << myid;
        ssltemp2>>filename_prefix;
#else // __MPI_ON__
        const  char * filename_prefix = filename_prefix_.c_str();
#endif // __MPI_ON__

        //initialize the random seed
        unsigned long seed;
        seed = 0;//time(0);
        // open logger 
        stringstream ssltemp;ssltemp<<filename_prefix<<"_logfile.txt";
        char filename_logger[250];ssltemp>>filename_logger;
        tmplogger.open(filename_logger);

#ifdef __NBODY_ON__	
        tmplogger << "Using NBODY single precision code \n";
#endif
#ifdef __CUNBODY_ON__
        tmplogger <<"Using CUNBODY code\n";
        tmplogger <<"\t *** WARNING CUNBODY doesn`t include the charge state of the particle,\n";
        tmplogger<<"\t *** it assumes the same charge when calculating the Coulomb interaction\n";
        cout <<"\t *** WARNING CUNBODY doesn`t include the charge state of the particle,\n";
        cout<<"\t *** it assumes the same charge when calculating the Coulomb interaction\n";
#endif

        //Print out the particle information in
        //true: each particle has his own file (default)
        //false: print out all information in one .._pAll.txt
        _cloud.use_particles_files(separateparticlefile_bool);

        //Initialize the ODE. First the filename_prefix, then the Order of the ODE can be 1,4 or 5, next is the initial timestep (should be around 1e-9), next is a boulean
        //true= with adaptive stepsize, false is without adaptive stepsize.
        InitIonFly(filename_prefix,ODEORDER,timestep,adaptive_stepsize,_cloud,odev);
#ifdef __MPI_ON__
        MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__  

        //print out the ions informations to the file every certain timestep
        SetPrintInterval(printinterval_timestep);

        //with or without Coulomb Interaction
        SetCoulomb(Coulomb_enable);

        //scaled Coulomb Factor, if used
        UseScaledCoulomb(Coulomb_scale_factor);

        InitCloud(NRPARTICLES,energy,seed,semiaxis_cloud,offset_cloud, Ions_cloud,_cloud,odev.forcev.trap_param);      

#ifdef __MPI_ON__ 
        MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__


#ifdef __MPI_ON__ 
        MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__	

        ops.launch(_cloud,odev);	

        ExitIonFly(_cloud);

        if(separateparticlefile_bool)
            BundleData(filename_prefix,NRPARTICLES+NSPARTICLES);

        if(myid==0) {
#ifdef __linux
            printf("%c[%d;%dmProgram Ended%c[%dm\n",27,4,31,27,0);
#else
            cout<<"Program ended"<<endl;
#endif
        }
#ifdef __MPI_ON__ 
        MPI::COMM_WORLD.Barrier();

        // Concatenate data	
        if(!separateparticlefile_bool)
            MPI_Concatenate_Data(NRPARTICLES+NSPARTICLES,filename_prefix);

        MPI::COMM_WORLD.Barrier();	
#endif // __MPI_ON__
    } //
    catch(const char * error) {
        cout<<"Error: "<<error<<endl;
        slogger << ERROR << error << SLogger::endmsg;
    }
}

void Run(INIReader & parser) {
#ifdef __MPI_ON__
    double startwtime = 0.0, endwtime;
    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI::Status status;
    MPI::Init();
    int myid = MPI::COMM_WORLD.Get_rank();
    int numprocs = MPI::COMM_WORLD.Get_size();
    MPI::Get_processor_name(processor_name,namelen);
    MPI::COMM_WORLD.Barrier();
    if (myid == 0)
        startwtime = MPI_Wtime();
#endif // __MPI_ON__

    IonCloud Cloud;
    _ode_vars Ode_vars;

#ifdef __MPI_ON__
    Ode_vars.create_mpiv();
#endif // _MPI_ON__
    DoSimulation(parser, Cloud, Ode_vars);
#ifdef __MPI_ON__
    MPI::COMM_WORLD.Barrier();
    if (myid == 0){
        endwtime = MPI_Wtime();
        std::cout << "wall clock time = " <<  endwtime-startwtime << std::endl;
    }
    MPI::Finalize();
#endif // __MPI_ON__
}
