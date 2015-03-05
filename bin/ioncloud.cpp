#include "ioncloud.h"
#include "MPI_simbuca.h"
LogFile clogger;

//to do with printing
//char *stream =new char[256];
char line[256];

IonCloud::~IonCloud() {
    //   for(int k=0;k < particles.size(); k++){
    //         delete streamvector[k];
    //   }
}

IonCloud::IonCloud() {
    //definition off static things...
    //private:
    sim_time_start = 0;
    globalstream = NULL;

    filenamebegin = ("xxx");

    particles_files = true;
    IDs.resize(0);
    nrparticles= 0;
    //particles.resize(0); //vector off size 0 he
    //initialvalues.resize(0);

    streamvector.resize(0); //so #streams = #particles
    nrcoll.resize(0);
    lifetime = 0.0;
}

void IonCloud::Create(const char* _filename) {
    lifetime=0.0;
    filenamebegin = _filename;
    sim_time_start = time(0);
    ions = new vector<Ion> ;
    images.push_back(pair< double,double> (0,0));
}

void IonCloud::Delete() {
    //time_t rawtime;
    //time ( &rawtime );
    time_t tm;
    tm = time(NULL);
    //printf(ctime(&tm));

    PrintParticles();
    if(particles_files){
        for(int k=0;k < particles.size(); k++){
            //IMPORTANT, Don`t delete this line, filesplitter uses the first * on the line to break the whole sequence
            (*streamvector[k])<< "#*******************************************************************************"<<endl;
            (*streamvector[k])<<"#simulation ended after "<<time(0)-sim_time_start<<" s\n";
            (*streamvector[k])<<"# #Particles left= "<<particles.size()<<endl;
            (*streamvector[k])<<"#total number of collisions = "<< nrcoll[k]<<endl;
            (*streamvector[k])<<"#End time of simulation: "<<ctime(&tm);
            (*streamvector[k]).flush();
            (*streamvector[k]).close();

        }
    }else{

        double tmpdouble=0.0;
        for(int k=0;k < particles.size(); k++){
            tmpdouble += nrcoll[k];
            //      (*streamvector[k]).flush();
        }
        //(*globalstream)<< globalstringstream.str();
        (*globalstream)<< "#*******************************************************************************"<<endl;
        (*globalstream)<<"#simulation ended after "<<time(0)-sim_time_start<<" s\n";
        (*globalstream)<<"# #Particles left= "<<particles.size()<<endl;
        (*globalstream)<<"#total number of collisions of all particles = "<< tmpdouble<<endl;
        (*globalstream)<<"#End time of simulation: "<<ctime(&tm);
        (*globalstream).flush();
        (*globalstream).close();
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

    //particles =  initialvalues
    for(int k=0;k < particles.size(); k++){
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

void IonCloud::use_particles_files(bool _bool) {
    particles_files = _bool;
}

void IonCloud::PrintMembers() {
    for(int j=0; j< ions->size(); j++){
        cout<<"ion ";cout<<j;clogger<<" ";cout<<(*ions)[j];cout<<"\n";
        cout<<"nrcoll = ";cout<<nrcoll[j];cout<<"\n";
    }
    clogger<<"cloud lifetime = ";clogger<<lifetime;clogger<<"\n";
}

void IonCloud::CreateFile() {
    //file initialiseren
    //IMPORTANT, Don`t delete this line, filesplitter uses the first * on the line to break the whole sequence



    int particleIndex = (particles.size()-1);
    stringstream ss (stringstream::in | stringstream::out);
    if (particles_files == true){
        ss<<filenamebegin<<"_p"<<(particleIndex+1)<<".txt";
        char filename[100]; //different from _filename (parameter) he!!
        ss >>filename;
        //cout<<filename<<endl; //just a check wich files are created...
        streamvector.push_back(new ofstream(filename));
        if(!(*streamvector[particleIndex]).is_open())
        {
            cout << "problem to open the file of particle " <<  particleIndex+1 << endl;
        }
        time_t rawtime2;
        time ( &rawtime2 );
        (*streamvector[particleIndex])<<"#Simulation started: "<<ctime (&rawtime2);
        (*streamvector[particleIndex])<<"#------------------------------------------------------------------------------#\n";
        (*streamvector[particleIndex])<<"#                Penning Trap Simulation Program by Simon Van Gorp             #\n";
        (*streamvector[particleIndex])<<"#      Dormand-Prince Runga-Kutta with Proportional Integrating controller     #\n" ;
        (*streamvector[particleIndex])<<"#------------------------------------------------------------------------------#\n";
        (*streamvector[particleIndex])<<"#index mass x y z (mm) vx vy vz (m/s) \t r+(mm) \t r-(mm) \t R(mm) \t Energy(eV) \t Temperature(K) \t t(ms) \n";
        (*streamvector[particleIndex])<<"#*******************************************************************************\n";
        //cout<<"stream for particles"<<particleIndex+1<<endl;
    }
    else{
        ss<<filenamebegin<<"_pAll.txt";
        char filename[50]; //different from _filename (parameter) he!!
        ss >>filename;
        if(globalstream == NULL){
            globalstream = new ofstream(filename);
            time_t rawtime2;
            time ( &rawtime2 );
            streamvector.push_back(new ofstream(filename));
            streamvector[particleIndex] = globalstream;
            (*streamvector[particleIndex])<<"#Simulation started: "<<ctime (&rawtime2);
            (*streamvector[particleIndex])<<"#------------------------------------------------------------------------------#\n";
            (*streamvector[particleIndex])<<"#                Penning Trap Simulation Program by Simon Van Gorp             #\n";
            (*streamvector[particleIndex])<<"#      Dormand-Prince Runga-Kutta with Proportional Integrating controller     #\n" ;
            (*streamvector[particleIndex])<<"#------------------------------------------------------------------------------#\n";
            (*streamvector[particleIndex])<< "#index mass x y z (mm) \t r+(mm) \t r-(mm) \t R(mm) \t Energy(eV) \t Temperature(K) \t t(ms) \n";
            (*streamvector[particleIndex])<<"#*******************************************************************************\n";

        }
        else{streamvector.push_back(globalstream);
        }

        //cout<<filename<<endl; //just a check wich files are created...
        //cout<<"no stream for particles"<<particleIndex+1<<endl;
    }
    //(*streamvector[particleIndex]).precision(5);
    (*streamvector[particleIndex]).setf(std::ios::scientific, std::ios::floatfield);
    //(*streamvector[particleIndex]).setf(std::ios::showpos);
}

void IonCloud::CloseFile(int _pindex, char* _reason) {
    //print particle for the last time he.
    //close the file

    (*streamvector[_pindex])<<"#******************************************************************************\n";
    (*streamvector[_pindex])<<"#Particle lost: "<<_reason<<endl;
    //(*streamvector[_pindex])<<"#Last particle position: \n#";
    //(*streamvector[_pindex])<<particles[_pindex]<<endl;
    time_t rawtime;struct tm * timeinfo;time ( &rawtime );timeinfo = localtime ( &rawtime );
    (*streamvector[_pindex])<<"# simulation ended after "<<time(0)-sim_time_start<<" s\n";
    (*streamvector[_pindex])<<"# #Particles= "<<particles.size()<<endl;
    (*streamvector[_pindex])<<"#total number of collisions = "<< nrcoll[_pindex]<<endl;
    (*streamvector[_pindex])<<"#End time of simulation: "<<asctime (timeinfo)<<endl;

    (*streamvector[_pindex])<<"#simulation ended after "<<time(0)-sim_time_start<<" s\n";
    //(*streamvector[_pindex]).close();
    //delete the particle out the streamvector he.
}

void IonCloud::AddParticle(Particle _p, Ion _i) {
    //    Particle tmpparticle(_x,_y,_z,_vx,_vy,_vz, _Ion);
    IDs.push_back(particles.size()); //dus ID = 0,1,2,3,4,5,...
    Particle * tmpptr = new Particle;
    *tmpptr = _p;
    particles.push_back(tmpptr);


    initialvalues.push_back(_p);
    ions->push_back(_i);
    nrcoll.push_back(0);
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
    nrcoll.erase(nrcoll.begin()+_index);
    charge.erase(charge.begin()+_index);
    mass.erase(mass.begin()+_index);
    wc.erase(wc.begin()+_index);
    wz2.erase(wz2.begin()+_index);
    streamvector.erase(streamvector.begin()+_index);

    // 2D arrays
    // here pos2 and vel2 are used as buffer

    // copy of the first part in buf
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
    // E=m*v*v/2
    // gemE = 3/2*kb*T
    //dus vx,vy,vz->E -> alle E's -> Egem -> T
    //if particles.size != 0
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
        //Radial Distance (mm)
        os<<" "<<sqrt((_cloud.pos[k][0]*_cloud.pos[k][0]+_cloud.pos[k][1]*_cloud.pos[k][1]))*1000;
        //time in (ms)
        os<<" "<< _cloud.lifetime*1000 ;

        os<<endl; //flush and endl
    }
    return os;
}

void IonCloud::InitializePoolVectors() {
    int n=particles.size();
    //    pos.resize(n);
    //    pos2.resize(n);
    //    vel.resize(n);
    //    vel2.resize(n);
    //    for(unsigned i=0;i<n;i++)
    //    {
    //        pos[i].resize(3);
    //        pos2[i].resize(3);
    //        vel[i].resize(3);
    //        vel2[i].resize(3);
    //    }
    mass.resize(n);
    charge.resize(n);
    wc.resize(n);
    wz2.resize(n);
    pos = new double[n][3];
    pos2 = new double[n][3];
    vel = new double[n][3];
    vel2 = new double[n][3];
    old_z.resize(n);
}

void IonCloud::CopyParticlesToVectors() {
    int n=particles.size();
    for(unsigned i=0;i<n;i++)
    {
        pos[i][0] = particles[i]->px;
        pos[i][1] = particles[i]->py;
        pos[i][2] = particles[i]->pz;
        pos2[i][0] = particles[i]->px;
        pos2[i][1] = particles[i]->py;
        pos2[i][2] = particles[i]->pz;
        vel[i][0] = particles[i]->pvx;
        vel[i][1] = particles[i]->pvy;
        vel[i][2] = particles[i]->pvz;
        vel2[i][0] = particles[i]->pvx;
        vel2[i][1] = particles[i]->pvy;
        vel2[i][2] = particles[i]->pvz;

        mass[i] = (*ions)[i].Getmass();
        charge[i] = (*ions)[i].Getcharge();
        wc[i] = (*ions)[i].Getwc();  //el_charge*B/temp.Getmass();
        wz2[i] = (*ions)[i].Getwz2();// (el_charge*Ud2)/temp.Getmass();
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

void IonCloud::UpdateIonParameters(_trap_param & trap_param) {
    int n=particles.size();
    for(unsigned i=0;i<n;i++)
    {
        (*ions)[i].SetParameters(trap_param);
    }
}
