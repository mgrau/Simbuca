#include <iostream>
#include <fstream>
#include <sstream>

#include "initialization.h"
#include "inireader.h"
#include "SLogger.h"

#if defined(__MPI_ON__)
#include <mpi.h>
#endif

#ifdef __NBODY_ON__
#include "nbody.h"
#endif //__NBODY_ON__

int main(int argc, char* argv[])
{
    if(argc>1)
    {
        try {
            INIReader parser(argv[1]);
            // SimParser sparser();
            Run(parser);
        }
        catch (const char * error) {
            cout << "could not open configuration file " << argv[1] << ". Error: " << error << endl;
        }
    }
    else
    {
        SLogger slogger("trapparameters");
        slogger << ERROR << "Simbuca requires an inputfile\n" << SLogger::endmsg;
        exit(1);
    }
    return 0;
}
