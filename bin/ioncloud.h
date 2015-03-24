#ifndef IONCLOUD_H // header guards
#define IONCLOUD_H

#include "globals.h"
#include "ion.h"
#include "particle.h"
#include "trapparameters.h"
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <sstream>
#include "logfile.h"
#include <time.h>

struct IonCloud
{      
    IonCloud();
    ~IonCloud();

    vector<Particle*> particles;
    vector<Particle> initialvalues; //to be able to reset the cloud p.e. for quadrupolescans!
    vector<Ion> * ions; 
    set<string> ion_types; 


    double (*pos)[3];
    double (*pos2)[3];
    double (*vel)[3];
    double (*vel2)[3];
    vector<double > old_z; // old z for Zpos print option
    vector < double > mass;
    vector < int > charge;
    int nrparticles;
    int initial_nrparticles;
    void InitializePoolVectors();
    void CopyVectorsToParticles();
    void CopyParticlesToVectors();

    double lifetime;

    std::ofstream* globalstream;
    std::vector<std::ofstream*> cloud_stream;
    std::vector<std::ofstream*> streamvector; //so #streams = #particles
    double sim_time_start;
    const char * filenamebegin;

    bool particles_files; //standard true
    vector<int> IDs;

    //structors!
    void Create(const char* _filename); 
    void Delete();   //also close all the files he!
    void Reset();
    void AddParticle(Particle _p, Ion _i);
    void PrintParticle(int &k);
    void PrintParticles();
    void PrintCloud();
    void PrintMembers();
    void DelParticle(int _index, char* _reason);
    double GetNrParticles();

    void use_particles_files(bool _bool);
    pair<double,double> Temperature();
    void CloseFile(int _pindex, char* _reason);
    void CreateFile();
#ifdef __MPI_ON__
    void UpdateIDs(int nrparticles);
#endif // __MPI_ON__
};

ostream& operator<<(ostream& os,IonCloud _cloud);

#endif
