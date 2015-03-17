//ionfly.cpp
#include "ionFly.h"
#include "MPI_simbuca.h"
#include "force.h"
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

void InitIonFly(const char * _filenamebegin, int _ode_order, double _timestep, bool _adaptive_stepsize,IonCloud &_cloud, _ode_vars & odev){ //IonFly stands above Ode and Coll...
    // ilogger<<"In case of ideal trap: B = ";ilogger<<odev.forcev.trap_param.B0;ilogger<<" T.U0/d2 = ";ilogger<< odev.forcev.trap_param.Ud2;ilogger<<"\n";
    ilogger << "Ideal Trap: Vrf = ";
    ilogger << odev.forcev.trap_param.Vrf;
    ilogger << ", Vdc = ";
    ilogger << odev.forcev.trap_param.Vdc;
    ilogger << "\n";
    odev.Init_ode( _ode_order,  _timestep,  _adaptive_stepsize);
    _cloud.Create(_filenamebegin);
    filename_begin=_filenamebegin;
}


void MoveParticles(double _time_movement,IonCloud &_cloud,_ode_vars &odev){
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
        if(myid==0) {
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

        // PRINT PARTICLES
        if(!odev.PrintAfterOperation&&!odev.PrintatZpos_bool)
        {
            if ( (_cloud.lifetime-starttime )> (print_interval*nr_interval)){
                _cloud.CopyVectorsToParticles();
                _cloud.PrintParticles();
                nr_interval++;
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
    }
    //printf("number of step : %d\n", GetCountStep() );
    if(odev.PrintAfterOperation)
    {
        _cloud.CopyVectorsToParticles();
        _cloud.PrintParticles();
    }
}


void AddParticle(double _x, double _y, double _z, double _vx, double _vy, double _vz, Ion& _Ion, IonCloud &_cloud){
    Particle tmpparticle(_x,_y,_z,_vx,_vy,_vz);
    _cloud.AddParticle(tmpparticle, _Ion);
}

void DelParticle(int _index,IonCloud &_cloud ){
    _cloud.DelParticle(_index, "nothing special, XXX");
}


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
void UseScaledCoulomb(double _ScaledCoulombFactor){
    scaledCoulombFactor = _ScaledCoulombFactor;
    ilogger<<"with scaled coulomb Interaction: factor ";
    ilogger<<scaledCoulombFactor;
    ilogger<<".\n";
}





void DoNoExcitation(double _time_movement, bool _buffergas, double _p_buffergas, IonCloud &_cloud, _ode_vars & odev){
    ilogger<<"No excitation for ";ilogger<<_time_movement;ilogger<<" sec,\n";

    MoveParticles(_time_movement,_cloud,odev);
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
static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    if ( (x != n) && (n<100 || (x % (n/100) != 0)) ) return;

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
