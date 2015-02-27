/*
 *  test.cpp
 *  gather xi for cunbody
 *  each process has an array xj[NTOT/NUMPROCS][3], the master process gather them in a xtotj[NTOT][3]
 *  Created by pierre dupre on 03/07/12.
 *  Copyright 2012 CEA/DSM/IRFU/SPP. All rights reserved.
 *
 */





#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "main.h"
/* Derived datatype example program. Sending subset of entries in an array */






using namespace std;

extern _mpi_vars mpi_vars;
extern _cloud cloud;
int main( int argc, char *argv[])
{
    init_mpi(argc,argv);
    
    int NRPARTICLES = atoi(argv[1]);
    
    Initialization(NRPARTICLES);

    // AT THIS POINT ALL NODES HAS AN ARRAY OF POSITIONS SUCH AS SIMBUCA
    MultiGPUs_CUNBODY_begin();
    
    MultiGPUs_CUNBODY_end();
    

    CompaCPU();
    
    clean_mpi();
}

