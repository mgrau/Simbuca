#ifdef __CUNBODY_ON__
#include <iostream>
#include <math.h>
#include "cunbody.h"
#include "ioncloud.h"
#include <string.h>
using namespace std;


_cunbody::_cunbody() {
    eps= 1e-23;
}

_cunbody::~_cunbody() {
    if(ai!=0)
        delete [] ai;
    if(m!=0)
        delete [] m;
    if(pos_tmp!=0)
        delete [] pos_tmp;
}

#ifdef __MPI_ON__
void _cunbody::InitMPI(_mpi_vars * _mpiv) {
    mpiv = _mpiv;
    AllocArrays(mpiv->n_part_tot);
    coords0 = mpiv->coords[0];
}
#endif //MPI

void _cunbody::AllocArrays(int nrparticles) {
    pos_tmp = new double[nrparticles][3];
    m = new double[nrparticles];
    for(int i=0;i<nrparticles;i++)
        m[i] =1;
#ifdef __MPI_ON__
    pos_tot = new double[nrparticles][3];
    ai_temp = new double[nrparticles][3];
#endif
    ai = new double[nrparticles][3];
}

void _cunbody::begin(const IonCloud &_cloud) {
#ifdef __MPI_ON__
    for(int dev=0;dev<mpiv->nb_gpu;dev++)
        MPI_Gatherv(_cloud.pos2[0], _cloud.nrparticles*3, MPI_DOUBLE, &pos_tot[0][0], mpiv->counts_3, mpiv->disps_3, MPI_DOUBLE,mpiv->ID_node_GPU[dev],mpiv->comm2d);

    if(mpiv->nb_gpu==1) {
        if(mpiv->coords[1]==0) {
            std::copy(pos_tot,pos_tot+mpiv->n_part_tot,pos_tmp);
            cunbody1_force_begin(pos_tot, m, pos_tmp, eps, ai_temp, mpiv->n_part_tot, mpiv->n_part_tot);
        }
    }
    else // some GPUs
    {
        // copy only the handled part 
        if(mpiv->coords[1]==0)
        {
            std::copy(pos_tot+mpiv->disps_Arow[coords0],pos_tot+mpiv->disps_Arow[coords0]+mpiv->counts_Arow[coords0],pos_tmp);
            cudaSetDevice(coords0);
            cunbody1_force_begin(pos_tmp, m, pos_tot, eps, ai, mpiv->n_part_tot, mpiv->counts_Arow[coords0]);
        }
    }
#else
    std::copy(_cloud.pos2,_cloud.pos2+_cloud.nrparticles,pos_tmp);
    cunbody1_force_begin(_cloud.pos2, m, pos_tmp, eps, ai, _cloud.nrparticles, _cloud.nrparticles);
#endif //MPI
}

void _cunbody::end(const IonCloud &_cloud) {
#ifdef __MPI_ON__
    if(mpiv->nb_gpu==1)
    {
        if(mpiv->coords[1]==0)
            cunbody1_force_end(pos_tot, m, pos_tmp, eps, ai_temp, mpiv->n_part_tot, mpiv->n_part_tot);
        MPI_Scatterv(&ai_temp[0][0],mpiv->counts_3, mpiv->disps_3, MPI_DOUBLE,ai[0],_cloud.nrparticles*3,MPI_DOUBLE,0,mpiv->comm2d);
    }
    else // some GPUs
    {
        if(mpiv->coords[1]==0)
        {
            cunbody1_force_end(pos_tmp, m, pos_tot, eps, ai, mpiv->n_part_tot, mpiv->counts_Arow[coords0]);
            MPI_Allreduce(&ai[0][0],ai_temp[0],mpiv->n_part_tot*3,MPI_DOUBLE,MPI_SUM,mpiv->comm_gpu);
        }
        MPI_Scatterv(&ai_temp[0][0],mpiv->counts_row, mpiv->disps_row_return, MPI_DOUBLE,ai[0],_cloud.nrparticles*3,MPI_DOUBLE,0,mpiv->comm1d_h);
    }
#else
    cunbody1_force_end(_cloud.pos2, m, pos_tmp, eps, ai, _cloud.nrparticles, _cloud.nrparticles);
#endif //MPI
}

#endif // __CUNBODY_ON__
