#include <iostream>
#include <fstream>
#include <sstream>

#include "initialization.h"
#include "SLogger.h"

#if defined(__MPI_ON__)
#include <mpi.h>
#endif

#ifdef __NBODY_ON__
#include "nbody.h"
#endif //__NBODY_ON__

int main(int argc,  char* argv[] )
{
    if(argc==2)
    {
        string simFile = argv[1];
        if(simFile.rfind(".sim")!=string::npos)
        {
            SimParser sparser(simFile.c_str());
            Run(argc,argv,sparser);
        }
    }
    else
    {
        SLogger slogger("trapparameters");
        slogger << ERROR << "Simbuca requires an inputfile\n" << SLogger::endmsg;exit(1);
    }
    return 0;
}
