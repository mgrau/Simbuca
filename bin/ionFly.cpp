//ionfly.cpp
#include "ionFly.h"
#include "MPI_simbuca.h"
#include "force.h"
bool buffergas = false;
bool imagecharges = false;
bool include_boundaries = true;
bool firstoperation = true;
LogFile ilogger;
double print_interval=0.01e-3; //in sec
int nr_interval = 0;
bool coulombinteraction;
double scaledCoulombFactor;

vector<double> im_q;
vector<double> im_i;
double im_i_sum = 0.0;
double t_prev = 0.0;
ofstream outfile;

const char * filename_begin;

#ifdef __GUI_ON__
Counter *percentage = new Counter;
#endif

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////////////////////Main Stuff///////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

void IncludeElectrodeBoundaries(bool _bool){
    include_boundaries = _bool;
    if(_bool){
        ilogger<<"With Electrode Boundaries\n";}
    else{ilogger<<"Without Electrode Boundaries\n";}
}



void InitIonFly(const char * _filenamebegin, int _ode_order, double _timestep, bool _adaptive_stepsize,IonCloud &_cloud, _ode_vars & odev){ //IonFly stands above Ode and Coll...
	ilogger<<"In case of ideal trap: B = ";ilogger<<odev.forcev.trap_param.B0;ilogger<<" T.U0/d2 = ";ilogger<< odev.forcev.trap_param.Ud2;ilogger<<"\n";
  	//InitOde(_ode_order, _timestep, _adaptive_stepsize);
	InitColl();
    	odev.Init_ode( _ode_order,  _timestep,  _adaptive_stepsize);
    	_cloud.Create(_filenamebegin);
	filename_begin=_filenamebegin;
}


void MoveParticles(double _time_movement,IonCloud &_cloud,_ode_vars &odev){
	//main lus

  // logger object
  //SLogger m_logger("MoveParticles");
  //SMsgType type = INFO;
  //SLogWriter::Instance()->SetMinType( type );
  // use e.g. : (see SMsgType.h for output levels)
   //m_logger << ERROR << "my error message" << SLogger::endmsg;
   //m_logger << INFO  << "my info message" << SLogger::endmsg;
   //m_logger << DEBUG  << "my debug message" << SLogger::endmsg;


    bool firststep = true;
	double endtime = _time_movement + _cloud.lifetime;
    odev.final_time = endtime;
	double starttime = 0;
    int myid =0;
#ifdef __MPI_ON__
    myid = MPI::COMM_WORLD.Get_rank();
    
#endif
	if(firstoperation)
	{
	    starttime = _cloud.lifetime;
        odev.initial_time = starttime;
	    firstoperation = false;
        //odev
        odev.Initpoolvectors(_cloud.nrparticles);
        odev.forcev.coulombinteraction = coulombinteraction;
        odev.forcev.scaledCoulombFactor = scaledCoulombFactor;
#ifdef __MPI_ON__
        odev.mpiv->ComputeCounts_Disps(_cloud.nrparticles);
#ifdef __CPUNBODY_ON__
        odev.mpiv->Init_GlobalPos(_cloud.nrparticles);
#else
        odev.mpiv->LoadCharge(_cloud.charge);
#endif // CPU
#endif // MPI
        
#ifdef __CUNBODY_ON__
#ifdef __MPI_ON__
        odev.cunbody.InitMPI(odev.mpiv);
#else
        odev.cunbody.AllocArrays(_cloud.nrparticles);
#endif //MPI
#endif // __CUNBODY_ON__
        
#ifdef __NBODY_ON__
#ifdef __MPI_ON__
        odev.nbody.InitMPI(odev.mpiv);
#endif
        odev.nbody.Initialization(_cloud, false); // same charge
#endif // __NBODY_ON__
        
        //cloud
        _cloud.InitializePoolVectors();
        _cloud.CopyParticlesToVectors();
        
#ifdef __CUNBODY_ON__
        if(odev.forcev.coulombinteraction)
        {
            odev.cunbody.begin(_cloud);
            odev.cunbody.end(_cloud);
        }
#endif // __CUNBODY_ON__
        
        
        //PRINT PARTICLES
        if(odev.PrintAfterOperation)
        {
            _cloud.CopyVectorsToParticles();
           // _cloud.PrintParticles();
        }
        
        //INIT IMAGE CHARGES
        if(imagecharges){
            im_q.resize(_cloud.nrparticles);
        	im_i.resize(_cloud.nrparticles);
            t_prev = 0.0 ;
       	    stringstream ssltemp;ssltemp<<filename_begin<<"_imagecharges.txt";
       	    char imagecharge_stream[250];ssltemp>>imagecharge_stream;
            outfile.open(imagecharge_stream);
            outfile<<"# (1)t\t\t\t(2)z(mm)\t\t\t(3)vz\t\t\t(4)im_q\t\t\t(5)im_i\t\t\t(6)EED\t\t\t(7)BesselCharge[q]\t\t\t(8)BesselCurrent[q/s]\t(9)BesselCharge[C]\t\t(10)BesselCurrent[A]\n";
            
            outfile.width(12);
            outfile.precision(10);
            outfile.setf(std::ios::scientific, std::ios::floatfield);
        }
	}
    
	bool collision;
	//_cloud.PrintParticles();
    odev.time_ini_ope =_cloud.lifetime;
    // PRINT
    
    double Energy_step;
	while (_cloud.lifetime < endtime) {
        //balint
        unsigned int nsteps = (unsigned int)((odev.total_time-starttime)/GetTimeStep(odev));
        unsigned int cstep = (unsigned int)((_cloud.lifetime-starttime)/GetTimeStep(odev));
        if(myid==0)
        {

	  loadbar(cstep, nsteps, 20);

        }
		if(_cloud.nrparticles == 0) break;
		//move particles
        
        if(odev.PrintatZpos_bool)
        {
            for(int i=0;i<_cloud.nrparticles;i++)
            {
                _cloud.old_z[i] = _cloud.pos[i][2];
            }
        }
        
        
		step(_cloud,odev);
        
		if(buffergas){
			for(int j=0; j< _cloud.nrparticles ; j++){
				collision = gas_interact(_cloud.vel[j][0],_cloud.vel[j][1],_cloud.vel[j][2],(*_cloud.ions)[j].Getmass(),GetTimeStep(odev));
				if(collision == true)
                {
                    _cloud.nrcoll[j]++;
                }
			}
		}//--buffergas
        
        
		//if CALCPOTENTIAL
		if(odev.forcev.trap_param.trap_config == 2){
			ChangePotential(_cloud.lifetime);
		}
        
		if(imagecharges){
			int n = _cloud.nrparticles;
            im_i_sum = 0.0 ;
			double z_sum = 0.0;
			double vz_sum = 0.0;
			double im_q_sum = 0.0;
			double EED = 0.0;
			double dt = _cloud.lifetime - t_prev;
            
            static double imageCharge = 0.0; //static so its scope is program runtime, therefore it holds its value over all integrations!! (need that to compute the current (see below))
            double imageCharge_new = 0.0;
            double imageCurrent = 0.0;
            
            
            for(int j_=0;j_<n;j_++){ // calculation for each particle of the no	de
                double q_prev = im_q[j_];
                double z_prev = _cloud.pos[j_][2];
                double vz_prev = _cloud.vel[j_][2];
                im_q[j_] = GetImageChargeApprox(_cloud.charge[j_], z_prev);
                
                im_i[j_] = (q_prev-im_q[j_])/dt; // /(z_prev-z);
                EED = im_q[j_]/im_i[j_]*vz_prev;
                im_i_sum += im_i[j_];
                im_q_sum += im_q[j_];
                z_sum += z_prev;
                vz_sum += vz_prev;
                //cout<<"image charge & image current \t"<< im_q[j_]<<"\t"<<im_i[j_]<<endl;
                
                //for each particle add the deposited image charge to get the overall image charge
                imageCharge_new += GetImageChargeBessel(_cloud.charge[j_], sqrt(pow(_cloud.pos[j_][0],2.0)+pow(_cloud.pos[j_][1],2.0)), _cloud.pos[j_][2]);
                
                
            }
			t_prev = _cloud.lifetime;
            
            
			//image Current is the time derivative of the image charge, so its (q_old-q_new)/delta_t
			imageCurrent = (imageCharge_new - imageCharge)/GetTimeStep(odev);
			//update imageCharge to computed new value
			imageCharge = imageCharge_new;
			//add image current to the cloud
			_cloud.AddImageCurrent(1.60217657E-19*imageCurrent);
            
			if(!firststep && ((_cloud.lifetime-starttime )> (print_interval*nr_interval)) ){
                //cout<<" total charge \t"<<im_i_sum<<endl;cin.get();
                outfile<<t_prev<<"\t"<<z_sum*1000.<<"\t"<<vz_sum<<"\t"<<im_q_sum<<"\t"<<im_i_sum<<"\t"<<EED<< "\t" << imageCharge << "\t" << imageCurrent << "\t" << 1.60217657E-19*imageCharge << "\t" << 1.60217657E-19*imageCurrent <<endl;
			}
			firststep = false;
		}//--imagecharges
        
        // PRINT PARTICLES
        if(!odev.PrintAfterOperation&&!odev.PrintatZpos_bool)
        {
            if ( (_cloud.lifetime-starttime )> (print_interval*nr_interval)){
                //cout << _cloud.lifetime << endl;
                _cloud.CopyVectorsToParticles();
                _cloud.PrintParticles();
                nr_interval++;
                //printf("%f %f\n", starttime , _cloud.lifetime);
            }
        }
        if(odev.PrintatZpos_bool)
        {
            for(int i=0;i<_cloud.nrparticles;i++)
            {
                if((_cloud.old_z[i]<odev.PrintZpos)&&(_cloud.pos[i][2]>odev.PrintZpos))
                {
                    _cloud.PrintParticle(i);
                }
            }
        }
        
		if(include_boundaries){
            //check if particles are out of boundary
            // WARNING THE LOOP HAS TO BE DECREMENTAL
            
            int n_temp = _cloud.nrparticles-1;
            for(int j= n_temp; j >=0 ;j--){   //check if radial size is ok,
                if((_cloud.pos[j][2] > 0.094) && (_cloud.pos[j][2] < 0.147) && (sqrt(_cloud.pos[j][0]*_cloud.pos[j][0]+_cloud.pos[j][1]*_cloud.pos[j][1]) > 0.001)){
                    cout<<"Particle hit pumping diaphragm with:\nindex = "<<j<<"\nx = "<<_cloud.pos[j][0]<<"\ny = "<<_cloud.pos[j][1]<<"\nz = "<<_cloud.pos[j][2]<<endl;
                    _cloud.DelParticle(j,"particle hit Pumping diaphragm");
                    if(j!=_cloud.nrparticles) j++;
                }
                else if (fabs(_cloud.pos[j][0]) > 0.039 || fabs(_cloud.pos[j][1]) > 0.039 || sqrt(_cloud.pos[j][0]*_cloud.pos[j][0]+_cloud.pos[j][1]*_cloud.pos[j][1]) > 0.039) {
                    cout<<"Particle hit electrode walls with:\nindex = "<<j<<"\nx = "<<_cloud.pos[j][0]<<"\ny = "<<_cloud.pos[j][1]<<"\nz = "<<_cloud.pos[j][2]<<endl;
                    _cloud.DelParticle(j,"particle hit electrode walls");
                    if(j!=_cloud.nrparticles) j++;
                }
                else if (_cloud.pos[j][2] < -0.29) {
                    cout<<"Particle hit bottom of the trap with:\nindex = "<<j<<"\nx = "<<_cloud.pos[j][0]<<"\ny = "<<_cloud.pos[j][1]<<"\nz = "<<_cloud.pos[j][2]<<endl;
                    _cloud.DelParticle(j,"particle at bottom of the trap.");
                    if(j!=_cloud.nrparticles) j++;
                }
                else if (_cloud.pos[j][2] > -0.001) {
                    cout<<"Particle hit top of the trap with:\nindex = "<<j<<"\nx = "<<_cloud.pos[j][0]<<"\ny = "<<_cloud.pos[j][1]<<"\nz = "<<_cloud.pos[j][2]<<endl;
                    _cloud.DelParticle(j,"particle at top of trap");
                    if(j!=_cloud.nrparticles) j++;
                }
            }
#ifdef __MPI_ON__
            // MPI variables must be re calculated if particles can be deleted
            odev.mpiv->ComputeCounts_Disps(_cloud.nrparticles);
#ifndef __CPUNBODY_ON__
            odev.mpiv->LoadCharge(_cloud.charge);
#endif
#endif // __MPI_ON__
        }
        
	}
	//printf("number of step : %d\n", GetCountStep() );
    if(odev.PrintAfterOperation)
    {
        _cloud.CopyVectorsToParticles();
        _cloud.PrintParticles();
    }
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////





void AddParticle(double _x, double _y, double _z, double _vx, double _vy, double _vz, Ion& _Ion, IonCloud &_cloud){
    Particle tmpparticle(_x,_y,_z,_vx,_vy,_vz);
    _cloud.AddParticle(tmpparticle, _Ion);
}

void DelParticle(int _index,IonCloud &_cloud ){
    _cloud.DelParticle(_index, "nothing special, XXX");
}


// ** FUNCTION SetCoulomb **//
// ** FUNCTION SetCoulomb **//
// ** FUNCTION SetCoulomb **//
void SetCoulomb(bool _coulombinteraction){
    coulombinteraction = _coulombinteraction;
    if(_coulombinteraction){
        ilogger<<"with Coulomb Interaction";
        ilogger<<"\n";
    }
    else{
        ilogger<<"without Coulomb Interaction";
        ilogger<<"\n";
    }
}
// ** FUNCTION UseScaledCoulomb **//
// ** FUNCTION UseScaledCoulomb **//
// ** FUNCTION UseScaledCoulomb **//
void UseScaledCoulomb(double _ScaledCoulombFactor){
    scaledCoulombFactor = _ScaledCoulombFactor;
    ilogger<<"with scaled coulomb Interaction: factor ";
    ilogger<<scaledCoulombFactor;
    ilogger<<".\n";
}





void DoNoExcitation(double _time_movement, bool _buffergas, double _p_buffergas, IonCloud &_cloud, _ode_vars & odev){
    buffergas = _buffergas; SetPressure(_p_buffergas);
    odev.InitExcitationVars(0,false, false,0.0, 0.0);
    ilogger<<"No excitation for ";ilogger<<_time_movement;ilogger<<" sec,\n";
    if(_buffergas){ilogger<<"with buffergas: ";ilogger<<_p_buffergas*1000.0;ilogger<<" mbar\n";}
    else{ilogger<<"with no buffergas\n";}
    
    MoveParticles(_time_movement,_cloud,odev);
}
void ChangeEfieldmap(char * _trapErz){
	ilogger<<"Read new electric fieldmap: ";ilogger<<_trapErz;ilogger<<"\n";
	ChangeEfield(_trapErz);
}

void DoTransfer(double _time_movement, bool _buffergas, double _p_buffergas, char * _transferErz, char * _trapErz, IonCloud &_cloud, _ode_vars & odev){
    buffergas = _buffergas; SetPressure(_p_buffergas);
    odev.InitExcitationVars(0,false, false,0.0, 0.0);
    ilogger<<"Doing Trap Transfer";ilogger<<_time_movement;ilogger<<" sec\n";
    ChangeEfield(_transferErz);
    MoveParticles(_time_movement,_cloud,odev);
    ChangeEfield(_trapErz);
}


void DoQuadrupoleExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, IonCloud &_cloud, _ode_vars & odev){
    buffergas = false;
    odev.InitExcitationVars(2,false,false,_U_exc, _exc_w);
    ilogger<<"Quadrupole no buffergas for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
}

void DoQuadrupoleExcitationWithBuffergas(double _time_movement, long double _exc_w, double _U_exc, double _p_buffergas, IonCloud &_cloud, _ode_vars & odev){
    buffergas = true; SetPressure(_p_buffergas);
    odev.InitExcitationVars(2,false,false,_U_exc, _exc_w);
    ilogger<<"Quadrupole and Buffergas ";ilogger<<_p_buffergas*1000.0;ilogger<<" mbar for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
    
}

void DoDipoleExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, IonCloud &_cloud, _ode_vars & odev){
    buffergas = false;
    odev.InitExcitationVars(1,false,false,_U_exc, _exc_w);
    ilogger<<"Dipole no Buffergas for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
}
void DoDipoleExcitationWithBuffergas(double _time_movement, long double _exc_w, double _U_exc, double _p_buffergas, IonCloud &_cloud, _ode_vars & odev){
    buffergas = true;  SetPressure(_p_buffergas);
    odev.InitExcitationVars(1, false,false,_U_exc, _exc_w);
    ilogger<<"Dipole and Buffergas ";ilogger<<_p_buffergas*1000.0;ilogger<<" mbar for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
}



void DoOctupoleExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, IonCloud &_cloud, _ode_vars & odev){
    buffergas = false;
    odev.InitExcitationVars(3,false,false,_U_exc, _exc_w);
    ilogger<<"Octupole no Buffergas for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
}
void DoOctupoleExcitationWithBuffergas(double _time_movement, long double _exc_w, double _U_exc, double _p_buffergas, IonCloud &_cloud, _ode_vars & odev){
    buffergas = true;  SetPressure(_p_buffergas);
    odev.InitExcitationVars(3,false,false,_U_exc, _exc_w);
    ilogger<<"Octupole and Buffergas ";ilogger<<_p_buffergas*1000.0;ilogger<<" mbar for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
}




void DoAxialCouplingExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, IonCloud &_cloud, _ode_vars & odev){
	buffergas = false;
	odev.InitExcitationVars(4,false,false,_U_exc,_exc_w);
	ilogger<<"Axial Coupling excitation no Buffergas for ";ilogger<<_time_movement;ilogger<<" sec\n";
	ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
    
}

void DoAxialQuadCoulpingExcitationWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, long double _exc_w2, double _U_exc2,
                                                   long double _exc_w3, double _U_exc3,long double _exc_w4, double _U_exc4, IonCloud &_cloud, _ode_vars & odev){
	buffergas = false;
	odev.InitExcitationVars(8,false,false,_U_exc,_exc_w);
    odev.U_exc2 = _U_exc2;
    odev.U_exc3 = _U_exc3;
    odev.U_exc4 = _U_exc4;
    odev.w_exc2 = _exc_w2;
    odev.w_exc3 = _exc_w3;
    odev.w_exc4 = _exc_w4;
    
	ilogger<<"SIMCO Axial Coupling excitation no Buffergas for ";ilogger<<_time_movement;ilogger<<" sec\n";
	ilogger<<"\t with axial quad frequency ";ilogger<<_exc_w2;ilogger<<" Hz and amplitude ";ilogger<<_U_exc2;ilogger<<"V\n";
    ilogger<<"\t with quadrupolar frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
    
}

void DoSIMCOWithoutBuffergas(double _time_movement, long double _exc_w, double _U_exc, long double _exc_w2, double _U_exc2, IonCloud &_cloud, _ode_vars & odev){
	buffergas = false;
	odev.InitExcitationVars(7,false,false,_U_exc,_exc_w);
    odev.U_exc2 = _U_exc2;
    odev.w_exc2 = _exc_w2;
	ilogger<<"SIMCO no Buffergas for ";ilogger<<_time_movement;ilogger<<" sec\n";
	ilogger<<"\t with dipolar frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    ilogger<<"\t with quadrupolar frequency ";ilogger<<_exc_w2;ilogger<<" Hz and amplitude ";ilogger<<_U_exc2;ilogger<<"V\n";
    MoveParticles(_time_movement,_cloud,odev);
}

void DoRotatingWall(int _order,double _time_movement, long double _exc_w, double _U_exc, bool _buffergas, double _p_buffergas, IonCloud &_cloud, _ode_vars & odev){
    buffergas = _buffergas;  SetPressure(_p_buffergas);
    switch ( _order )
    {  case 0:
            odev.InitExcitationVars(1,true,false,_U_exc, _exc_w);
            ilogger<<" No Excitation RW ";
            break;
        case 1:
            odev.InitExcitationVars(1,true,false,_U_exc, _exc_w);
            ilogger<<"Dipole RW ";
            break;
        case 2:
            odev.InitExcitationVars(2,true,false,_U_exc, _exc_w);
            ilogger<<"Quadrupole RW ";
            break;
        case 3:
            odev.InitExcitationVars(3,true,false,_U_exc, _exc_w);
            ilogger<<"Octupole RW ";
            break;
        default:
            odev.InitExcitationVars(2,true,false,_U_exc, _exc_w);
            SLogger slogger("ionfly");
            slogger << ERROR << "choose a correct order for Rotating wall." << SLogger::endmsg;exit(1);
    }
    
    if(buffergas){  ilogger<<"and Buffergas ";ilogger<<_p_buffergas*1000.0;ilogger<<" mbar for ";}
    else {ilogger<<"without Buffergas for ";}
    
    ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    
    MoveParticles(_time_movement,_cloud,odev);
}


void DoAntiRW(int _order,double _time_movement, long double _exc_w, double _U_exc, bool _buffergas, double _p_buffergas, IonCloud &_cloud, _ode_vars & odev){
    buffergas = _buffergas;  SetPressure(_p_buffergas);
    switch ( _order )
    {
        case 1:
            odev.InitExcitationVars(1,false,true,_U_exc, _exc_w);
            ilogger<<"anti-RW (Dipole)";
            break;
        case 2:
            odev.InitExcitationVars(2,false,true,_U_exc, _exc_w);
            ilogger<<"anti-RW (Quadrupole)";
            break;
        default:
            odev.InitExcitationVars(2,false,true,_U_exc, _exc_w);
            SLogger slogger("ionfly");
            slogger << ERROR << "choose a correct order for anti-rotating wall." << SLogger::endmsg;exit(1);
    }
    
    if(buffergas){  ilogger<<"and Buffergas ";ilogger<<_p_buffergas*1000.0;ilogger<<" mbar for ";}
    else {ilogger<<"without Buffergas for ";}
    
    ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with frequency ";ilogger<<_exc_w;ilogger<<" Hz and amplitude ";ilogger<<_U_exc;ilogger<<"V\n";
    
    MoveParticles(_time_movement,_cloud,odev);
}


void DoARexcitation(double _time_movement, double _U_exc,double _exc_w, double _exc_w2, IonCloud & _cloud,_ode_vars & odev){
    // in case of AR: set sweeprate = _exc_wc2
    buffergas = false;
    odev.InitExcitationVars(9,false,false,_U_exc,_exc_w);
    odev.w_exc2 = (_exc_w2-_exc_w)/_time_movement; //sweeprate
    InitARexcitation("../../fieldmaps/E_fieldmap_U8_1V.txt",filename_begin);
    ilogger<<"AR excitation for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with amplitude ";ilogger<<_U_exc;ilogger<<" V\n";
    ilogger<<"\t frequency sweep ";ilogger<<_exc_w;ilogger<<" -> ";ilogger<<_exc_w2;ilogger<<" Hz\n ";
    MoveParticles(_time_movement,_cloud,odev);
}

void DoFBexcitation(double _time_movement, double _U_exc,double _exc_w, double _exc_w2, IonCloud & _cloud,_ode_vars & odev){
    buffergas = false;
    imagecharges = true;
    odev.InitExcitationVars(10,false,false,_U_exc,_exc_w);
    odev.w_exc2 = _exc_w2; 
    InitFBexcitation("../../fieldmaps/E_fieldmap_U8_1V.txt",filename_begin);
    ilogger<<"FB excitation for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with amplitude ";ilogger<<_U_exc;ilogger<<" V\n";
    ilogger<<"\t frequency 1 ";ilogger<<_exc_w;ilogger<<" Hz\n ";
    ilogger<<"\t frequency 2 ";ilogger<<_exc_w2;ilogger<<" Hz\n ";
    MoveParticles(_time_movement,_cloud,odev);
}

/*
void DoFBexcitation(double _time_movement, IonCloud & _cloud,_ode_vars & odev){
    buffergas = false;
    imagecharges = true;
    odev.InitExcitationVars(10,false,false,_U_exc,_exc_w);
    odev.w_exc2 = _exc_w2;
    InitFBexcitation("../../fieldmaps/E_fieldmap_U8_1V.txt",filename_begin);
    ilogger<<"FB excitation for ";ilogger<<_time_movement;ilogger<<" sec\n";
    ilogger<<"\t with amplitude ";ilogger<<_U_exc;ilogger<<" V\n";
    ilogger<<"\t frequency 1 ";ilogger<<_exc_w;ilogger<<" Hz\n ";
    ilogger<<"\t frequency 2 ";ilogger<<_exc_w2;ilogger<<" Hz\n ";
    MoveParticles(_time_movement,_cloud,odev);
}
*/

void DoFrequencyScan(double _time_movement, double _U_exc, int _order, bool _rotatingwall, double _p_buffergas, double _central_freq, double _freq_dev, int _nr_steps, IonCloud &_cloud, _ode_vars & odev){
    buffergas = true; SetPressure(_p_buffergas);
    if(_order > 3){ilogger<<"ERROR order should be < 3 , You shall choose the correct order!\n";}
    odev.InitExcitationVars(_order, _rotatingwall, false, _U_exc, _central_freq);
    
    
    //file initialiseren
	stringstream ss;
	ss<<"Frequencyscan.txt";
	ofstream myfile;
	myfile.open (ss.str().c_str(), ios::trunc);
	print_interval=1000.0; //so never print to particle file.
    double sim_start=time(0);
    
    time_t rawtime2;
    time ( &rawtime2 );
    myfile<<"#Simulation started: "<<ctime (&rawtime2);
    myfile<<"#------------------------------------------------------------------------------#\n";
    myfile<<"#                Penning Trap Simulation Program by Simon Van Gorp             #\n";
    myfile<<"#      Dormand-Prince Runga-Kutta with Proportional Integrating controller     #\n" ;
    myfile<<"#------------------------------------------------------------------------------#\n";
    myfile<<"# Doing ";
    switch ( _order )
    {  case 1:
            ilogger<<"Dipole-";
            break;
        case 2:
            ilogger<<"Quadrupole-";
            break;
        case 3:
            ilogger<<"Octupole-";
            break;
        default:
            ilogger<<"ERROR in ionFly.cpp in DoFrequencyScan: choose a correct order for Rotating wall Frequency scans!!! ";
    }
    myfile<<"excitatie. Met amplitude "<< _U_exc<<endl;
    myfile<<"#Around "<< _central_freq <<" in "<<_nr_steps<<" steps off "<<(_freq_dev*2) /_nr_steps<<" Hz."<<endl;
    myfile<<"#buffergaspressure = "<<_p_buffergas*1000.0<< " mbar. Simulation for "<<_time_movement<<" sec."<<endl;
    myfile << "# \t ";
    for(int i=0;i<_cloud.nrparticles;i++){
        myfile<<" \t "<<(*_cloud.ions)[i];
    }
    myfile<<endl<< "#wq(2pi*Hz) \t R_particle(mm)"<<endl;
    ilogger<<"Frequency Scan around central frequency ";
    ilogger<< _central_freq/(10e6);ilogger<<" Mhz.";ilogger<<" In ";ilogger<<_nr_steps;
    ilogger<<" steps off ";ilogger<<(_freq_dev*2) /_nr_steps;ilogger<<" Hz.\n";
    
    
    double w_min = _central_freq - _freq_dev;
    double w_max = _central_freq + _freq_dev;
    double w_stepsize = (_freq_dev*2) /_nr_steps;
    
    /*for loop over all the frequencys*/
    for(double w_step=w_min; w_step < w_max; w_step +=w_stepsize){
        ilogger<<"Taking step ";ilogger<<w_step;ilogger<<"\n";
        
        odev.w_exc = w_step;
        /*take step*/
        MoveParticles(_time_movement,_cloud,odev);
        
        /*write output*/
        myfile<<w_step-_central_freq<<"\t"; //freqshift
        for(int i=0; i<_cloud.nrparticles;i++){
            myfile<<sqrt(_cloud.pos[i][0]*_cloud.pos[i][0]+_cloud.pos[i][1]*_cloud.pos[i][1])<<"\t";
            myfile<<" ";
        }
        myfile<<endl;
        /*set Ioncloud back to inital positsion */
        _cloud.Reset();
    }
    
    
    time_t rawtime;struct tm * timeinfo;time ( &rawtime );timeinfo = localtime ( &rawtime );
    myfile<< "******************************************************"<<endl;
    myfile<<"#simulation ended after "<<time(0)-sim_start<<" s\n";
    myfile<<"#End time of simulation: "<<asctime (timeinfo)<<endl;
    
    myfile.close();
}

void SetPrintInterval(double _print_interval){
    print_interval = _print_interval;
}

void ExitIonFly(IonCloud &_cloud){
    ilogger<<"Close everything...\n";
    //ExitForce();
    _cloud.Delete();
}

void Print_particles(IonCloud &_cloud){
    _cloud.PrintParticles();
}

void use_particle_file(bool _bool,IonCloud &_cloud){
    _cloud.use_particles_files(_bool);
}

// Set the lifetime of the ioncloud
void SetLifetime(double _lifetime,IonCloud &_cloud)
{
    _cloud.lifetime = _lifetime;
}

void SetFileNamePrefix(char * filenameprefix){

}
void SetTotalTime_of_Simu(double t_,_ode_vars & odev){
    odev.total_time = t_;
}
void SetBufferGas(bool b_){
    buffergas = b_;
}
static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    if ( (x != n) && (x % (n/100) != 0) ) return;
    
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
#ifdef __linux__
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
#endif

#ifdef __GUI_ON__
  percentage->setValue((int)(ratio*100));
#endif
}

#ifdef __GUI_ON__
void SetPercentagePointer(Counter *c){
  percentage = c;
}
#endif
