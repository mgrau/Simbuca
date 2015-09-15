#include "ioncloud.h"
#include "MPI_simbuca.h"
LogFile clogger;

//to do with printing
//char *stream =new char[256];
char line[256];

IonCloud::~IonCloud() {}

IonCloud::IonCloud() {
    //definition off static things...
    sim_time_start = 0;
    globalstream = NULL;
    print_x_stream = NULL;

    filenamebegin = ("xxx");

    particles_files = true;
    IDs.resize(0);
    nrparticles= 0;

    cloud_stream.resize(0); //so #streams = #particles
    streamvector.resize(0); //so #streams = #particles
    lifetime = 0.0;
}

void IonCloud::Create(const char* _filename) {
    lifetime=0.0;
    filenamebegin = _filename;
    sim_time_start = time(0);
    ions = new vector<Ion> ;
}

void IonCloud::Delete() {
    time_t tm;
    tm = time(NULL);

    PrintParticles();
    if(particles_files){
        for(int k=0;k < particles.size(); k++) {
            (*streamvector[k])<< "#*******************************************************************************"<<endl;
            (*streamvector[k])<<"#simulation ended after "<<time(0)-sim_time_start<<" s\n";
            (*streamvector[k])<<"# #Particles left= "<<particles.size()<<endl;
            (*streamvector[k])<<"#End time of simulation: "<<ctime(&tm);
            (*streamvector[k]).flush();
            (*streamvector[k]).close();
        }
    }
    else {
        (*globalstream)<< "#*******************************************************************************"<<endl;
        (*globalstream)<<"#simulation ended after "<<time(0)-sim_time_start<<" s\n";
        (*globalstream)<<"# #Particles left= "<<particles.size()<<endl;
        (*globalstream)<<"#End time of simulation: "<<ctime(&tm);
        (*globalstream).flush();
        (*globalstream).close();
    }

    for (int k=0; k < ion_types.size(); k++) {
        (*cloud_stream[k]).flush();
        (*cloud_stream[k]).close();
    }

    if (print_x_stream != NULL) {
        (*print_x_stream).flush();
        (*print_x_stream).close();
    }

    clogger<<"#simulation ended after ";clogger<<time(0)-sim_time_start;clogger<<" s\n";
    cout<<"#simulation ended after "<<time(0)-sim_time_start<<" s ";
#ifdef __MPI_ON__
    cout<<"on CPU "<<MPI::COMM_WORLD.Get_rank()<<endl;
#else // __MPI_ON__
    cout<<endl;
#endif // __MPI_ON__
}

void IonCloud::Reset() {
    lifetime=0.0;
    for(int k=0;k < particles.size(); k++) {
        *particles[k] = initialvalues[k];
    }
}

void IonCloud::PrintParticle(int &k) {
    double Energy = (*ions)[k].Getmass()*(vel[k][0]*vel[k][0]+vel[k][1]*vel[k][1]+vel[k][2]*vel[k][2])*0.5*Joule_to_eV;
    if(particles_files) {
        (*streamvector[k])<<IDs[k]<<" ";
        (*streamvector[k])<<(*ions)[k];
        (*streamvector[k])<<" "<< lifetime*1000.0 ;
        (*streamvector[k])<<" "<<pos[k][0]*1000.0;
        (*streamvector[k])<<" "<<pos[k][1]*1000.0;
        (*streamvector[k])<<" "<<pos[k][2]*1000.0;
        (*streamvector[k])<<" "<< vel[k][0];
        (*streamvector[k])<<" "<< vel[k][1];
        (*streamvector[k])<<" "<< vel[k][2];
        (*streamvector[k])<<" "<<Energy;
        (*streamvector[k])<<endl;
    }
    else {
        (*globalstream)<<IDs[k]<<" ";
        (*globalstream)<<(*ions)[k];
        (*globalstream)<<" "<< lifetime*1000.0 ;
        (*globalstream)<<" "<<pos[k][0]*1000.0;
        (*globalstream)<<" "<<pos[k][1]*1000.0;
        (*globalstream)<<" "<<pos[k][2]*1000.0;
        (*globalstream)<<" "<< vel[k][0];
        (*globalstream)<<" "<< vel[k][1];
        (*globalstream)<<" "<< vel[k][2];
        (*globalstream)<<" "<<Energy;
        (*globalstream)<<endl;
    }
}

void IonCloud::PrintParticles() {
    for(int k=0;k < particles.size(); k++)
        PrintParticle(k);
}

void IonCloud::PrintX(int k) {
    if (print_x_stream == NULL) { 
        string print_x_filename(filenamebegin);
        print_x_filename.append("_x.txt");
        print_x_stream = new ofstream(print_x_filename.c_str());
        (*print_x_stream)<<"#Print at X position";
        (*print_x_stream)<<"#------------------------------------------------------------------------------#\n";
        (*print_x_stream)<<"#                          JILA EDM Paul Trap Simulation                       #\n";
        (*print_x_stream)<<"#      Dormand-Prince Runga-Kutta with Proportional Integrating controller     #\n" ;
        (*print_x_stream)<<"#------------------------------------------------------------------------------#\n";
        (*print_x_stream)<<"#index \t type \t #particules \t t (ms) \t x0 y0 z0 (mm) sx sy sz (m/s) \t Energy(eV)\n";
        (*print_x_stream)<<"#*******************************************************************************\n";
        (*print_x_stream).setf(std::ios::scientific, std::ios::floatfield);
    }
    else {
        (*print_x_stream)<<IDs[k]<<" ";
        (*print_x_stream)<<(*ions)[k];
        (*print_x_stream)<<" "<< lifetime*1000.0 ;
        (*print_x_stream)<<" "<<pos[k][0]*1000.0;
        (*print_x_stream)<<" "<<pos[k][1]*1000.0;
        (*print_x_stream)<<" "<<pos[k][2]*1000.0;
        (*print_x_stream)<<" "<< vel[k][0];
        (*print_x_stream)<<" "<< vel[k][1];
        (*print_x_stream)<<" "<< vel[k][2];
        (*print_x_stream)<<endl;
    }
}

void IonCloud::PrintCloud() {
    for (set<string>::iterator type = ion_types.begin(); type!=ion_types.end(); ++type) {
        int id = distance(ion_types.begin(),type);
        double x0 = 0, y0 = 0, z0 = 0;
        double x1 = 0, y1 = 0, z1 = 0;
        double energy = 0.0;
        for (int i=0; i<particles.size(); i++)
            if ((*ions)[i].name == *type) {
                x0 += pos[i][0];
                y0 += pos[i][1];
                z0 += pos[i][2];
                energy += 0.5*Joule_to_eV*(*ions)[i].mass*(vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][1]);
            }
        x0 /= particles.size();
        y0 /= particles.size();
        z0 /= particles.size();
        for (int i=0; i<particles.size(); i++)
            if ((*ions)[i].name == *type) {
                x1 += (pos[i][0] - x0)*(pos[i][0] - x0);
                y1 += (pos[i][1] - y0)*(pos[i][1] - y0);
                z1 += (pos[i][2] - z0)*(pos[i][2] - z0);
            }
        x1 = sqrt(x1/particles.size());
        y1 = sqrt(y1/particles.size());
        z1 = sqrt(z1/particles.size());

        (*cloud_stream[id]) << id <<" ";
        (*cloud_stream[id]) << *type << " ";
        (*cloud_stream[id]) << particles.size() <<" ";
        (*cloud_stream[id]) << lifetime*1e3 << " ";
        (*cloud_stream[id]) << 1000.0*x0 << " " << 1000.0*y0 << " " << 1000.0*z0 << " ";
        (*cloud_stream[id]) << 1000.0*x1 << " " << 1000.0*y1 << " " << 1000.0*z1 << " ";
        (*cloud_stream[id]) << energy;
        (*cloud_stream[id]) << endl;
    }
}

void IonCloud::use_particles_files(bool _bool) {
    particles_files = _bool;
}

void IonCloud::PrintMembers() {
    for(int j=0; j< ions->size(); j++){
        cout<<"ion ";cout<<j;clogger<<" ";cout<<(*ions)[j];cout<<"\n";
    }
    clogger<<"cloud lifetime = ";clogger<<lifetime;clogger<<"\n";
}

void IonCloud::CreateFile() {
    if (particles.size()>0) {
        stringstream ss (stringstream::in | stringstream::out);
        if (particles_files) {
            int particleIndex = (particles.size()-1);
            ss<<filenamebegin<<"_p"<<(particleIndex+1)<<".txt";
            char filename[100]; //different from _filename (parameter) he!!
            ss >>filename;
            streamvector.push_back(new ofstream(filename));
            if(!(*streamvector[particleIndex]).is_open())
            {
                cout << "problem to open the file of particle " <<  particleIndex+1 << endl;
            }
            time_t rawtime2;
            time ( &rawtime2 );
            (*streamvector[particleIndex])<<"#Simulation started: "<<ctime (&rawtime2);
            (*streamvector[particleIndex])<<"#------------------------------------------------------------------------------#\n";
            (*streamvector[particleIndex])<<"#                          JILA EDM Paul Trap Simulation                       #\n";
            (*streamvector[particleIndex])<<"#      Dormand-Prince Runga-Kutta with Proportional Integrating controller     #\n" ;
            (*streamvector[particleIndex])<<"#------------------------------------------------------------------------------#\n";
            (*streamvector[particleIndex])<<"#index \t mass \t t (ms) \t x y z (mm) vx vy vz (m/s) \t Energy(eV)\n";
            (*streamvector[particleIndex])<<"#*******************************************************************************\n";
            (*streamvector[particleIndex]).setf(std::ios::scientific, std::ios::floatfield);
        }
        else {
            ss<<filenamebegin<<"_pAll.txt";
            char filename[100]; //different from _filename (parameter) he!!
            ss >>filename;
            if(globalstream == NULL) {
                globalstream = new ofstream(filename);
                time_t rawtime2;
                time ( &rawtime2 );
                (*globalstream)<<"#Simulation started: "<<ctime (&rawtime2);
                (*globalstream)<<"#------------------------------------------------------------------------------#\n";
                (*globalstream)<<"#                          JILA EDM Paul Trap Simulation                       #\n";
                (*globalstream)<<"#      Dormand-Prince Runga-Kutta with Proportional Integrating controller     #\n" ;
                (*globalstream)<<"#------------------------------------------------------------------------------#\n";
                (*globalstream)<<"#index \t mass \t t (ms) \t x y z (mm) vx vy vz (m/s) \t Energy(eV)\n";
                (*globalstream)<<"#*******************************************************************************\n";
                (*globalstream).setf(std::ios::scientific, std::ios::floatfield);
            }
        }
        if (cloud_stream.size() < ion_types.size()) {
            stringstream ss (stringstream::in | stringstream::out);
            string cloud_filename(filenamebegin);
            int ion_type = ion_types.size() - 1;
            ss << "_cloud" << ion_type << ".txt";
            cloud_filename.append(ss.str());
            cloud_stream.push_back(new ofstream(cloud_filename.c_str()));
            time_t rawtime2;
            time ( &rawtime2 );
            (*cloud_stream[ion_type])<<"#Simulation started: " << ctime(&rawtime2);
            (*cloud_stream[ion_type])<<"#------------------------------------------------------------------------------#\n";
            (*cloud_stream[ion_type])<<"#                          JILA EDM Paul Trap Simulation                       #\n";
            (*cloud_stream[ion_type])<<"#      Dormand-Prince Runga-Kutta with Proportional Integrating controller     #\n" ;
            (*cloud_stream[ion_type])<<"#------------------------------------------------------------------------------#\n";
            (*cloud_stream[ion_type])<<"#index \t type \t #particules \t t (ms) \t x0 y0 z0 (mm) sx sy sz (m/s) \t Energy(eV)\n";
            (*cloud_stream[ion_type])<<"#*******************************************************************************\n";
            (*cloud_stream[ion_type]).setf(std::ios::scientific, std::ios::floatfield);
        }
    }
}

void IonCloud::CloseFile(int _pindex, char* _reason) {
    //print particle for the last time he.
    (*streamvector[_pindex])<<"#******************************************************************************\n";
    (*streamvector[_pindex])<<"#Particle lost: "<<_reason<<endl;
    time_t rawtime;struct tm * timeinfo;time ( &rawtime );timeinfo = localtime ( &rawtime );
    (*streamvector[_pindex])<<"# simulation ended after "<<time(0)-sim_time_start<<" s\n";
    (*streamvector[_pindex])<<"# #Particles= "<<particles.size()<<endl;
    (*streamvector[_pindex])<<"#End time of simulation: "<<asctime (timeinfo)<<endl;
    (*streamvector[_pindex])<<"#simulation ended after "<<time(0)-sim_time_start<<" s\n";
}

void IonCloud::AddParticle(Particle _p, Ion _i) {
    IDs.push_back(particles.size()); //dus ID = 0,1,2,3,4,5,...
    Particle * tmpptr = new Particle;
    *tmpptr = _p;
    particles.push_back(tmpptr);

    initialvalues.push_back(_p);
    ions->push_back(_i);
    ion_types.insert(_i.name);
    nrparticles++;
    CreateFile();
}

void IonCloud::DelParticle(int _index, char* _reason) {
    if (particles_files == true)
    {
        PrintParticle(_index);
        CloseFile(_index,_reason); //this function first, prints the particle for the last time.
    }

    particles.erase(particles.begin()+_index);//remove the index-1'th particle

    ions->erase(ions->begin()+_index);
    charge.erase(charge.begin()+_index);
    mass.erase(mass.begin()+_index);
    streamvector.erase(streamvector.begin()+_index);

    if(_index!=0)
    {
        std::copy(pos,pos+_index,pos2);
        std::copy(vel,vel+_index,vel2);
    }
    //copy of the second part in buf
    std::copy(pos+_index+1,pos+nrparticles,pos2+_index);
    std::copy(vel+_index+1,vel+nrparticles,vel2+_index);
    //recopy buf in main array
    nrparticles--;
    std::copy(pos2,pos2+nrparticles,pos);
    std::copy(vel2,vel2+nrparticles,vel);

    IDs.erase(IDs.begin()+_index);

    clogger<<"Deleted particle nr ";clogger<<initialvalues.size()-particles.size();clogger<<" with ID: ";
    clogger<<IDs[_index];clogger<<" because ";clogger<<_reason;clogger<<" @";clogger<<lifetime;clogger<<" sec.\n";
}

pair<double,double> IonCloud::Temperature() {
    pair<double,double> mean_energy(0.0,0.0);

    if (particles.size() != 0){
        for(int k=0; k<particles.size();k++){
            mean_energy.first += (*ions)[k].Getmass()*(vel[k][0]*vel[k][0]+vel[k][1]*vel[k][1]);
            mean_energy.second += (*ions)[k].Getmass()*(vel[k][2]*vel[k][2]);
        }
        mean_energy.first *=(0.5/particles.size());
        mean_energy.second *=(0.5/particles.size());
        pair<double,double> Temperature = mean_energy;
        Temperature.first *= 2.0/(3.0*kb);
        Temperature.second *= 2.0/(3.0*kb);
        //    cout<<"mean energy is: "<<mean_energy<<endl;
        //         cout<<"Temperature is: "<<Temperature<<endl;
        return Temperature;
    }
    else
        return mean_energy;
}
#ifdef __MPI_ON__
void IonCloud::UpdateIDs(int nrparticles) {
    int myid = MPI::COMM_WORLD.Get_rank();
    int numprocs = MPI::COMM_WORLD.Get_size();
    int * counts = new int[numprocs];
    int * displ = new int[numprocs];
    int last_part = nrparticles % numprocs;
    displ[0] = 0;
    for(int i=0;i<numprocs;i++)
    {
        counts[i] = nrparticles/numprocs;
        if(i<last_part)
            counts[i] += 1;
    }
    for(int i=1;i<numprocs;i++)
    {
        displ[i] = displ[i-1] +counts[i-1];
    }

    for(int k=0;k < particles.size(); k++)
    {
        IDs[k] = IDs[k] + displ[myid];
    }
    delete [] counts;
    delete [] displ;

    return;
}
#endif // __MPI_ON__

ostream& operator<<(ostream& os,IonCloud _cloud) {
    for(int k=0; k<_cloud.particles.size();k++){
        //index off the particle in the particles vector
        os<<k<<" ";
        //write away ion
        os<<" "<<(*_cloud.ions)[k];
        //write away x,y,z (mm)
        os<<" "<<_cloud.pos[k][0]*1000;
        os<<" "<<_cloud.pos[k][1]*1000;
        os<<" "<<_cloud.pos[k][2]*1000;
        //time in (ms)
        os<<" "<< _cloud.lifetime*1000 ;
        os<<endl;
    }
    return os;
}

void IonCloud::InitializePoolVectors() {
    int n=particles.size();
    mass.resize(n);
    charge.resize(n);
    pos = new double[n][3];
    pos2 = new double[n][3];
    vel = new double[n][3];
    vel2 = new double[n][3];
    old_x.resize(n);
}

void IonCloud::CopyParticlesToVectors() {
    int n=particles.size();
    for(unsigned i=0;i<n;i++)
    {
        pos[i][0]  = particles[i]->px;
        pos[i][1]  = particles[i]->py;
        pos[i][2]  = particles[i]->pz;
        pos2[i][0] = particles[i]->px;
        pos2[i][1] = particles[i]->py;
        pos2[i][2] = particles[i]->pz;
        vel[i][0]  = particles[i]->pvx;
        vel[i][1]  = particles[i]->pvy;
        vel[i][2]  = particles[i]->pvz;
        vel2[i][0] = particles[i]->pvx;
        vel2[i][1] = particles[i]->pvy;
        vel2[i][2] = particles[i]->pvz;

        mass[i] = (*ions)[i].Getmass();
        charge[i] = (*ions)[i].Getcharge();
    }
}

void IonCloud::CopyVectorsToParticles() {
    int n=particles.size();
    for(unsigned i=0;i<n;i++)
    {
        particles[i]->px = pos[i][0];
        particles[i]->py = pos[i][1];
        particles[i]->pz = pos[i][2];

        particles[i]->pvx = vel[i][0] ;
        particles[i]->pvy = vel[i][1] ;
        particles[i]->pvz = vel[i][2] ;
    }
}
