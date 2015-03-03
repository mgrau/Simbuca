//CreateRecoilDistr.h
#ifndef CREATERECOILDISTR_H
#define CREATERECOILDISTR_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "mtrand.h"
#include <math.h>
#include "globals.h"
#include <time.h>

void InitRecoilDistr();

void CreateRecoilSpectraDifferential(char * filename);
//writes a recoil spectrum to a file 
//and generates (If I feel like it) smalla gnuplot file 

void GenerateRandomEnergy(char * filename,int N);
//writes out a file with N random energies written to it.

double MonteCarloEnergy();
//throws back a certain energy value in eV. Accordinly to the recoil energy spectra.

#endif
