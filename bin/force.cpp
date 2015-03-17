#include "force.h"
#include "ode.h"
#ifndef __CPUNBODY_ON__
#include "builtin_types.h"
#endif // __CPUNBODY_ON__
#include "MPI_simbuca.h"

#ifdef __NBODY_ON__
#include "nbody.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
#endif // _NBODY_ON__

const bool useTree = false;
//if useTree = false -> Cunbody is used.

//Pool of objects so constructor doesn`t need to be calles so often!
//located here, because returned by reference
vector< vector<double> > pool_force;
vector<double> e_field(3);
vector<double> B_scalar(3);
bool poolvectorsinitialized = false;

ofstream tmpstream;
double tmpinterval;

LogFile flogger;

//    Adjust the settings of these matrices: rows,columns
//    For example for the WITCH experiment: E:(1000,1000), B:(701,26)
//    For the ASACUSA experiment: E:(1000,1000), B:(500,500)
// fieldmap Efieldmap(1000,1000);
// fieldmap Bfieldmap(1000,1000);
// fieldmap Electrodefieldmap(1000,1000); //electrode to do the excitation on

#ifdef __OCTGRAV_ON__
#define NMAX (2048) // NMAX = total number of particles in the cloud
static octgrav tree;
std::vector< float4 > pos(NMAX);
std::vector< float4 > acc(NMAX);
std::vector< float3 > vel(NMAX);
std::vector< float3 > vel1(NMAX);
#endif // __OCTGRAV_ON__

#ifdef __MPI_ON__
int numprocs; // number of nodes
int myid; // id of the node
#endif //__MPI_ON__

int npart; //  number of particles in the sub cloud at the beginning
int npart_tot; // total number of particle in the cloud at the beginning
int jj_;

bool asynchro = true;

void Initpoolvector(int nrParticles,const IonCloud &_cloud,_force_vars &forcev){
    pool_force.resize(nrParticles);

    poolvectorsinitialized = true;
#ifdef __OCTGRAV_ON__
    pos.resize(nrParticles);
    acc.resize(nrParticles);
    vel.resize(nrParticles);
    vel1.resize(nrParticles);
    //initialize the tree
    const double eps = 1e-23;
    const float theta = 0.2;
    tree.set_softening(eps);
    tree.set_opening_angle(theta);

    flogger<<"Using the Octtree code. \n";
#endif // __OCTGRAV_ON__
}

void force(const IonCloud &_cloud,_ode_vars &odev){
    double x,y,z,vx,vy,vz,r2,r3,r5,z2,r4,z4, qDmassa,qqDmassaXke, el2ke;
    double t = _cloud.lifetime;
    _force_vars *forcev = &odev.forcev;
    int trap_type = (*forcev).trap_param.trap_type;

    double kVdc = forcev->trap_param.kappa * forcev->trap_param.Vdc;
    double Vrf = forcev->trap_param.Vrf;
    double r02 = forcev->trap_param.r0 * forcev->trap_param.r0;
    double z02 = forcev->trap_param.z0 * forcev->trap_param.z0;
    double omega_rf = forcev->trap_param.freq_rf*2*pi;

    if(!poolvectorsinitialized){Initpoolvector(_cloud.nrparticles,_cloud,odev.forcev);}

    el2ke = el_charge*el_charge*ke; //charges should be taken into account!


    if ((*forcev).coulombinteraction) {        
#ifdef __CUNBODY_ON__
        odev.cunbody.begin(_cloud);
#endif // __CUNBODY_ON__
#ifdef __NBODY_ON__
        odev.nbody.force_begin(_cloud);
#endif // __NBODY_ON__
    }

    for(unsigned j_=0; j_< _cloud.nrparticles; j_++) {
        x=_cloud.pos2[j_][0];
        y=_cloud.pos2[j_][1];
        z=_cloud.pos2[j_][2];
        vx=_cloud.vel2[j_][0];
        vy=_cloud.vel2[j_][1];
        vz=_cloud.vel2[j_][2];

#ifdef __OCTGRAV_ON__
        //init tree stuff
        pos[j_].x=x;
        pos[j_].y=y;
        pos[j_].z=z;
        pos[j_].w=1.0;
#endif // __OCTGRAV_ON__

        // cout << "omega_rf = " << omega_rf << endl;
        // cout << "Vrf*cos(omega_rf*t) = " << Vrf*cos(omega_rf*t) << endl;

        if (trap_type == 0) {
            qDmassa = _cloud.charge[j_]*el_charge/_cloud.mass[j_];
            (*forcev).derivs[j_][0] = qDmassa*(kVdc/z02 - Vrf*cos(omega_rf*t)/r02)*x;
            (*forcev).derivs[j_][1] = qDmassa*(kVdc/z02 + Vrf*cos(omega_rf*t)/r02)*y;
            (*forcev).derivs[j_][2] = qDmassa*(-2.0*kVdc/z02)*z;
        }            
    }//end off loop over particles



    if ((*forcev).coulombinteraction){      
#ifdef __CUNBODY_ON__
        odev.cunbody.end(_cloud);
        for(int j_=0;j_ <_cloud.nrparticles;j_++)
        {
            qqDmassaXke = (*forcev).scaledCoulombFactor*el2ke/ _cloud.mass[j_];
            (*forcev).derivs[j_][0] -= qqDmassaXke*odev.cunbody.ai[j_][0];
            (*forcev).derivs[j_][1] -= qqDmassaXke*odev.cunbody.ai[j_][1];
            (*forcev).derivs[j_][2] -= qqDmassaXke*odev.cunbody.ai[j_][2];
        }
#endif // __CUNBODY_ON__


#ifdef __NBODY_ON__
        odev.nbody.force_end( _cloud);        
#ifdef __MPI_ON__
        //MPI_Barrier(odev.mpiv->comm2d);
        for(int j_=0;j_ <_cloud.nrparticles;j_++)
        {
            qqDmassaXke = (*forcev).scaledCoulombFactor*el2ke/ _cloud.mass[j_];
            (*forcev).derivs[j_][0] -= qqDmassaXke*odev.nbody.host_acc_sub[j_].x;
            (*forcev).derivs[j_][1] -= qqDmassaXke*odev.nbody.host_acc_sub[j_].y;
            (*forcev).derivs[j_][2] -= qqDmassaXke*odev.nbody.host_acc_sub[j_].z;
        }

#else
        for(int j_=0;j_ <_cloud.nrparticles;j_++)
        {
            qqDmassaXke = (*forcev).scaledCoulombFactor*el2ke/ _cloud.mass[j_];
            (*forcev).derivs[j_][0] -= qqDmassaXke*odev.nbody.host_acc[j_].x;
            (*forcev).derivs[j_][1] -= qqDmassaXke*odev.nbody.host_acc[j_].y;
            (*forcev).derivs[j_][2] -= qqDmassaXke*odev.nbody.host_acc[j_].z;
        }
#endif // __MPI_ON__
#endif // __NBODY_ON__

#ifndef __MPI_ON__
#ifdef __OCTGRAV_ON__
        tree.evaluate_gravity(pos, acc);
        for(j_=0;j_ <_cloud.nrparticles;j_++){
            qqDmassaXke = el2ke/ _cloud.mass[j_];
            (*forcev).derivs[j_][0] += -(*forcev).scaledCoulombFactor*qqDmassaXke*mj[j_]*acc[j_].x;
            (*forcev).derivs[j_][1] += -(*forcev).scaledCoulombFactor*qqDmassaXke*mj[j_]*acc[j_].y;
            (*forcev).derivs[j_][2] += -(*forcev).scaledCoulombFactor*qqDmassaXke*mj[j_]*acc[j_].z;
        }
#endif
#endif

        /*** CPU ***/
#ifdef __CPUNBODY_ON__
        double fx,fy,fz,r_ij, kTqiTqjDmj;
        double bjx, bjy,bjz;
        double rx,ry,rz;
        double distSqr;
        double invDist;
        double invDistCube;

#ifdef __MPI_ON__
        MPI_Allgatherv(_cloud.pos2[0], _cloud.nrparticles*3, MPI_DOUBLE, &odev.mpiv->pos_tot[0][0], odev.mpiv->counts_3, odev.mpiv->disps_3, MPI_DOUBLE,odev.mpiv->comm2d);

        for(int j_=0;j_<_cloud.nrparticles;j_++) // calculation for each particle of the node
        {
            fx=fy=fz=0.0;
            r_ij=0.0,
                kTqiTqjDmj=(*forcev).scaledCoulombFactor*ke*_cloud.charge[j_]*el_charge*el_charge/_cloud.mass[j_];
            //new iterator jj_ taking in account particles of previous nodes
            //jj_ = j_ + MPI_displ[myid]; 
            bjx = _cloud.pos2[j_][0];
            bjy = _cloud.pos2[j_][1];
            bjz = _cloud.pos2[j_][2];

            /*calculate Coulomb normal*/
            for(unsigned int i_=0;i_ <npart_tot;i_++) // calculation of the distance between jj_ and the whole particles
            {
                rx = odev.mpiv->pos_tot[i_][0] -bjx;
                ry = odev.mpiv->pos_tot[i_][1] -bjy;
                rz = odev.mpiv->pos_tot[i_][2] -bjz;
                distSqr = rx*rx+ ry*ry + rz*rz + 1.e-40;
                invDist = 1./sqrt(distSqr); //eps included to make sure that if r=0 no error is thrown.
                invDistCube = (invDist*invDist*invDist);
                fx += rx*invDistCube*kTqiTqjDmj;//*(*_cloud.ions)[i_].Getcharge();
                fy += ry*invDistCube*kTqiTqjDmj;//*(*_cloud.ions)[i_].Getcharge();
                fz += rz*invDistCube*kTqiTqjDmj;//*(*_cloud.ions)[i_].Getcharge();
            }
            (*forcev).derivs[j_][0] -=fx;
            (*forcev).derivs[j_][1] -=fy;
            (*forcev).derivs[j_][2] -=fz;
        }
#else // __MPI_ON__
        //cpu calculation of the coulomb force
        kTqiTqjDmj=0.0;
        for(unsigned j_=0;j_ <_cloud.nrparticles;j_++){
            fx=fy=fz=0.0;
            r_ij=0.0;
            kTqiTqjDmj=(*forcev).scaledCoulombFactor*ke*el_charge*_cloud.charge[j_]*el_charge/_cloud.mass[j_];
            bjx = _cloud.pos2[j_][0];
            bjy = _cloud.pos2[j_][1];
            bjz = _cloud.pos2[j_][2];
            /*calculate Coulomb normal*/
            for(unsigned int i_=0;i_ <_cloud.nrparticles;i_++){
                if(i_ != j_){
                    rx = _cloud.pos2[i_][0] -bjx;
                    ry = _cloud.pos2[i_][1] -bjy;
                    rz = _cloud.pos2[i_][2] -bjz;
                    distSqr = rx*rx+ ry*ry + rz*rz;
                    invDist = 1./sqrt(distSqr); //eps included to make sure that if r=0 no error is thrown.
                    invDistCube = (invDist*invDist*invDist);
                    fx += rx*invDistCube*kTqiTqjDmj*(*_cloud.ions)[i_].Getcharge();
                    fy += ry*invDistCube*kTqiTqjDmj*(*_cloud.ions)[i_].Getcharge();
                    fz += rz*invDistCube*kTqiTqjDmj*(*_cloud.ions)[i_].Getcharge();
                }
            }
            (*forcev).derivs[j_][0] -=fx;
            (*forcev).derivs[j_][1] -=fy;
            (*forcev).derivs[j_][2] -=fz;
        }
#endif // _MPI_ON__
#endif // __CPUNBODY_ON__
    } //end if forcev.coulombinteraction
}
