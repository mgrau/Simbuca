#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include "ion.h"
#include "trapparameters.h"
#include "logfile.h"
#include "ionFly.h"
#include "IonTable.h"
#include "PDGTable.h"
#include <stdio.h>
#include <string.h>
#include "SLogger.h"

#ifdef __MPI_ON__
void MPI_Concatenate_Data(int NRPARTICLES, char * filename_prefix);
void Rename_part_file(const char* filenamebegin,int nrparticle);
#endif // __MPI_ON__
int ImportData(const char * filename_prefix,IonTable Table,PDGTable pdgTable, IonCloud &cloud, _trap_param & trap_param); // return the total number of particle
void BundleData(const char* filenamebegin, int nrparticles);
void PrintIonCloudinfo(const char *beginstream, int nrparticles, double diaphragm_radius_mm);
void PrintIonCloudGaussEvo(const char *beginstream, int nrparticles, double diaphragm_radius_mm, bool skipfirstparticle);
