#include "ion_fly.h"
#include "MPI_simbuca.h"
#include "force.h"
bool firstoperation = true;
LogFile ilogger;
double print_interval=0.01e-3; //in sec
int nr_interval = 0;
bool coulomb_interaction;
double coulomb_scale;

vector<double> im_q;
vector<double> im_i;
double im_i_sum = 0.0;
double t_prev = 0.0;
ofstream outfile;

void move_particles(double _time_movement,IonCloud &_cloud,_ode_vars &odev) {
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
        odev.forcev.coulomb_interaction = coulomb_interaction;
        odev.forcev.coulomb_scale = coulomb_scale;
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
        cout << "Init NBODY" << endl;
#ifdef __MPI_ON__
        odev.nbody.InitMPI(odev.mpiv);
#endif
        odev.nbody.Initialization(_cloud, false); // same charge
        cout << "Done Init NBODY" << endl;
#endif // __NBODY_ON__
        //cloud
        _cloud.InitializePoolVectors();
        _cloud.CopyParticlesToVectors();
#ifdef __CUNBODY_ON__
        if(odev.forcev.coulomb_interaction) {
            odev.cunbody.begin(_cloud);
            odev.cunbody.end(_cloud);
        }
#endif // __CUNBODY_ON__


        //PRINT PARTICLES
        if(odev.print_after_operation)
            _cloud.CopyVectorsToParticles();
    }

    odev.time_ini_ope =_cloud.lifetime;
    // PRINT

    double Energy_step;
    while (_cloud.lifetime < endtime) {
        unsigned int nsteps = (unsigned int)((odev.total_time-starttime)/GetTimeStep(odev));
        unsigned int cstep = (unsigned int)((_cloud.lifetime-starttime)/GetTimeStep(odev));
        if(myid==0)
            loadbar(cstep, nsteps, 20);
        if(_cloud.nrparticles == 0)
            break;

        if(odev.print_at_x)
            for(int i=0; i<_cloud.nrparticles; i++)
                _cloud.old_x[i] = _cloud.pos[i][0];

        step(_cloud,odev);

        // PRINT PARTICLES
            if ( (_cloud.lifetime-starttime )> (print_interval*nr_interval)){
                _cloud.CopyVectorsToParticles();
                _cloud.PrintParticles();
                _cloud.PrintCloud();
                nr_interval++;
            }
        if(odev.print_at_x)
            for(int i=0; i<_cloud.nrparticles; i++)
                if(((_cloud.old_x[i]<odev.print_x)&&(_cloud.pos[i][0]>odev.print_x)) || ((_cloud.old_x[i]>odev.print_x)&&(_cloud.pos[i][0]<odev.print_x))) {
                    _cloud.PrintX(i);
                }
    }
    if(odev.print_after_operation)
    {
        _cloud.CopyVectorsToParticles();
        _cloud.PrintParticles();
        _cloud.PrintCloud();
    }
}

void InitIonFly(const char * _filenamebegin, int _ode_order, double _timestep, bool _adaptive_stepsize,IonCloud &_cloud, _ode_vars & odev) {
    ilogger << "Ideal Trap: Vrf = ";
    ilogger << odev.forcev.trap_param.Vrf;
    ilogger << ", Vdc = ";
    ilogger << odev.forcev.trap_param.Vdc;
    ilogger << "\n";
    odev.Init_ode( _ode_order,  _timestep,  _adaptive_stepsize);
    _cloud.Create(_filenamebegin);
}

void ExitIonFly(IonCloud &_cloud) {
    ilogger<<"Close everything...\n";
    _cloud.Delete();
}

void AddParticle(double _x, double _y, double _z, double _vx, double _vy, double _vz, Ion& _Ion, IonCloud &_cloud) {
    Particle tmpparticle(_x,_y,_z,_vx,_vy,_vz);
    _cloud.AddParticle(tmpparticle, _Ion);
}

void normal_operation(double time, IonCloud &cloud, _ode_vars & odev){
    ilogger<<"Normal trap operation for ";ilogger<<time;
    ilogger<<" sec,\n";
    move_particles(time,cloud,odev);
}

void SetCoulomb(bool _coulomb_interaction){
    coulomb_interaction = _coulomb_interaction;
    if(_coulomb_interaction){
        ilogger<<"with Coulomb Interaction";
        ilogger<<"\n";
    }
    else{
        ilogger<<"without Coulomb Interaction";
        ilogger<<"\n";
    }
}
void UseScaledCoulomb(double _scale_factor){
    coulomb_scale = _scale_factor;
    ilogger<<"with scaled coulomb Interaction: factor ";
    ilogger<<coulomb_scale;
    ilogger<<".\n";
}
void SetPrintInterval(double _print_interval){
    print_interval = _print_interval;
}

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50) {
    if ( (x != n) && (n<100 || (x % (n/100) != 0)) ) return;
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
#ifdef __linux__
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
#endif
}
