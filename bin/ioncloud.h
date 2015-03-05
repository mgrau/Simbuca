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
#include <sstream>
#include "logfile.h"
#include <time.h>


class IonCloud
{      
    public:
        IonCloud();
        ~IonCloud();

<<<<<<< HEAD
        //global vars 
=======
>>>>>>> temp
        vector<Particle*> particles;
        vector<Particle> initialvalues; //to be able to reset the cloud p.e. for quadrupolescans!
        vector<Ion> * ions; 
        vector< pair<double,double> > images; //Images with x: time(s) and y: induced current(A)

        std::vector<std::ofstream*> streamvector; //so #streams = #particles

        double (*pos)[3];
        double (*pos2)[3];
        double (*vel)[3];
        double (*vel2)[3];
        vector<double > old_z; // old z for Zpos print option
        // mass charge
        vector < double > mass;
        vector < int > charge;
        vector < double > wc;
        vector < double > wz2;
        int nrparticles;
        int initial_nrparticles;
        void InitializePoolVectors();
        void CopyVectorsToParticles();
        void CopyParticlesToVectors();

        // Update eigen frequency
<<<<<<< HEAD
        void UpdateIonParameters(_trap_param & trap_param);
        vector<int> nrcoll;
        double lifetime;
        //double time; //is the same for everyone he!
=======
        vector<int> nrcoll;
        double lifetime;
>>>>>>> temp

        //structors!
        void Create(const char* _filename); 
        void Delete();   //also close all the files he!
        void Reset();
<<<<<<< HEAD
        //accessable functions
=======
>>>>>>> temp
        void AddParticle(Particle _p, Ion _i);
        void PrintParticle(int &k);
        void PrintParticles();
        void PrintMembers();
        void DelParticle(int _index, char* _reason);
        double GetNrParticles();

        void use_particles_files(bool _bool);
        pair<double,double> Temperature();
        std::ofstream* globalstream;
        double sim_time_start;
        const char * filenamebegin;
        void CloseFile(int _pindex, char* _reason);

<<<<<<< HEAD
        void AddImageCurrent(double _image);
        double GetImageCurrent(double delaytime) const;
        //to be implemented      int nrStepChanges;
        //to be implemented      bool alive;
        bool particles_files; //standard true
        //if false: print everything in one file.
=======
        bool particles_files; //standard true
>>>>>>> temp
        vector<int> IDs;
#ifdef __MPI_ON__
        void UpdateIDs(int nrparticles);
#endif // __MPI_ON__

    private:
        void CreateFile();
};

ostream& operator<<(ostream& os,IonCloud _cloud);

#endif
