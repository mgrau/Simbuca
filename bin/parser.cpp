#include "parser.h"
#include "MPI_simbuca.h"
#include "IonTable.h"
#include "ion.h"
using namespace std;

LogFile plogger;   

int ImportData(const char * filename_prefix,IonTable Table,PDGTable pdgTable, IonCloud &cloud, _trap_param & trap_param) {
#ifdef  __MPI_ON__  
    // MPI
    int myid = MPI::COMM_WORLD.Get_rank();
    int numprocs = MPI::COMM_WORLD.Get_size();
#else 
    int myid=0;  
#endif // __MPI_ON__  
    //file
    char line[256];
    char c = 'a';
    ifstream input_file;
    stringstream ssfilename_input;
    string filename_temp;
    streampos input_pos;

    // data 
    int npart_tot;
    // containers
    vector< double > x_;
    vector< double > y_;
    vector< double > z_;
    vector< double > vx_;
    vector< double > vy_;
    vector< double > vz_;
    vector< double > mass_;
    vector< string > name_;
    vector<Ion > Ion_;
    Ion Ion_tmp;
    double index,mass,x,y,z,vx,vy,vz,rplus,rmin,R,Energy,Temp,t;
    string name;
    vector<double > mass_frac;
    vector<string > name_frac;
    vector<int > frac;
    bool new_mass;
    fstream tmp_file;
    double tmp;
    int ii;
    if(myid==0)
    {   
        filename_temp = filename_prefix;
#ifdef  __MPI_ON__ 
        //filename_temp.erase(filename_temp.size() - 2, filename_temp.size()-1); // remove the "-0"
#endif // __MPI_ON__
        ssfilename_input << filename_temp << "_pAll.txt";
        input_file.open(ssfilename_input.str().c_str(),ios::in);
        if(!input_file)
        {
            cout << "Problem: file " << ssfilename_input.str() << " doesn't exist" << endl;     
            //return;
            exit(EXIT_FAILURE);
        }
        else
        {
            //cout << ssfilename_input.str() << " is open " << endl;   
        }
        input_file.seekg(0,ios::end); // place the get pointer at the end


        // research of the second '=' from the end
        for(int i=0 ;i<2 ;i++)
        {
            while(c!='=')
            {
                input_file.unget();
                c =input_file.peek();
                //	cout << c;
            }
            //cout  << endl;
            c='a';
            //input_pos = input_file.tellg();
        }
        // extract the number of particle

        input_file.ignore(256,'='); 
        input_file >> npart_tot;

        cout << "number of remained particles = " << npart_tot << endl;

        // go backward during npart_tot lines + 3 lines 
        for(int i=0 ;i<npart_tot + 3  ;i++)
        {
            while(c!='\n')
            {
                input_file.unget();
                c =input_file.peek();
                //cout << c;
            }
            //cout  << endl;
            c='a';
            //input_pos = input_file.tellg();
        }
        // read these lines
        input_file.getline(line,256);
        while((c!='#')&&(!input_file.eof()))
        {
            input_file>>index>>name>>x>>y>>z>>vx>>vy>>vz>>rplus>>rmin>>R>>Energy>>Temp>>t;input_file.getline(line,256);
            //cout <<  index <<" "<<name<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<t<<endl;
            c = input_file.peek();
            // fill the containers;
            x_.push_back(x/1000.);
            y_.push_back(y/1000.);
            z_.push_back(z/1000.);
            vx_.push_back(vx);
            vy_.push_back(vy);
            vz_.push_back(vz);
            name_.push_back(name);
            if(Table.TestIonName(name)){
                mass = Table.ResearchMass(name);
            }else if(pdgTable.TestPDGName(name)){
                int pdgid = pdgTable.GetPDGId(name);
                mass = pdgTable.GetPDGMass(pdgid)*1000.0/mass_Mev; // PDG table: need to convert here GeV to amu
            }else{
                cout << "Particle name in ImportData <" << name << "> not found! Exiting... :(" << endl;
                exit(-1);
            }		
            mass_.push_back(mass);
            Ion_tmp.SetParameters(mass);
            Ion_tmp.SetName(name);
            Ion_.push_back(Ion_tmp);
        }
        input_file.close();

        tmp_file.open("tmp_file_name.txt",ios::out);
        for(int i=0;i<name_.size();i++)
        {
            tmp_file << name_[i] << "\n" ;
        }
        tmp_file.close();
        //for(int i=0;i<npart_tot;i++)
        //{
        //cout <<  i <<" "<<name_[i]<<" "<<x_[i]<<" "<<y_[i]<<" "<<z_[i]<<" "<<vx_[i]<<" "<<vy_[i]<<" "<<vz_[i] <<endl;
        //cout <<  i <<" "<<Ion_[i].Getmass()/amu<<" "<<x_[i]<<" "<<y_[i]<<" "<<z_[i]<<" "<<vx_[i]<<" "<<vy_[i]<<" "<<vz_[i] <<endl;
        //}


        // Research of species
        mass_frac.push_back(mass_[0]);
        name_frac.push_back(name_[0]);
        for(int i=1;i<npart_tot;i++)
        {
            new_mass = false;
            for(unsigned int j=0;j<mass_frac.size();j++)
            {
                if(mass_[i]!=mass_frac[j])
                {
                    new_mass = true;
                }
                else
                {
                    new_mass = false;
                }			
            }
            if(new_mass)
            {
                mass_frac.push_back(mass_[i]);
                name_frac.push_back(name_[i]);
            }	
        }
        // Calculation of the fraction
        //frac.resize(mass_frac.size());
        //for(unsigned int j=0;j<mass_frac.size();j++)
        //{
        //	frac[j] = 0;
        //	for(int i=0;i<npart_tot;i++)
        //	{
        //		if(mass_[i]==mass_frac[j])
        //		{
        //		frac[j]++;
        //		}
        //	}

        //}
        // Set the good fraction	
        //for(int i=0;i<npart_tot;i++)
        //{
        //	for(unsigned int j=0;j<frac.size();j++)
        //	{
        //		if(mass_[i]==mass_frac[j])
        //		{
        //			tmp =  (double )frac[j]/npart_tot*100.;

        //cout << tmp << endl;
        //			Ion_[i].SetFraction(round(tmp));
        //		}
        //	}
        //cout << i << " " << Ion_[i].GetFraction() << endl;
        //}

    } // end myid==0


#ifdef __MPI_ON__ 
    // Bcast npart_tot
    MPI::COMM_WORLD.Bcast(&npart_tot,1,MPI::DOUBLE,0);
    if(myid!=0)
    {
        x_.resize(npart_tot);
        y_.resize(npart_tot);
        z_.resize(npart_tot);
        vx_.resize(npart_tot);
        vy_.resize(npart_tot);
        vz_.resize(npart_tot);
        mass_.resize(npart_tot);
        name_.resize(npart_tot);
        tmp_file.open("tmp_file_name.txt",ios::in);
        if(!tmp_file&&myid==1)
        {
            SLogger slogger("parser");
            slogger << ERROR << "Error in parsing MPI Importdatafunction" << SLogger::endmsg;exit(EXIT_FAILURE);
        }
        else
        {
            for(int i=0;i<name_.size();i++)
            {
                tmp_file >> name_[i];
            }
        }
    }
    // Bcast data
    MPI::COMM_WORLD.Bcast(&x_[0],x_.size(),MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&y_[0],y_.size(),MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&z_[0],z_.size(),MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&vx_[0],vx_.size(),MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&vy_[0],vy_.size(),MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&vz_[0],vz_.size(),MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&mass_[0],mass_.size(),MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&t,1,MPI::DOUBLE,0);
    // Add Particle in each node
    int * counts = new int[numprocs];
    int * displ = new int[numprocs];
    int last_part = npart_tot % numprocs;
    displ[0] = 0;
    for(int i=0;i<numprocs;i++)
    {
        counts[i] = npart_tot/numprocs;
        if(i<last_part)
            counts[i] += 1;
    }
    for(int i=1;i<numprocs;i++)
    {
        displ[i] = displ[i-1] +counts[i-1];
    }

    for(int i=0;i < counts[myid] ;i++)
    {
        ii = i + displ[myid];
        Ion_tmp.SetParameters(mass_[ii]);
        Ion_tmp.SetName(name_[ii]);
        AddParticle(x_[ii],y_[ii],z_[ii],vx_[ii],vy_[ii],vz_[ii],Ion_tmp,cloud);
    }
    cloud.UpdateIDs(npart_tot);
    delete [] counts;
    delete [] displ;
#else // __MPI_ON__
    for(int i=0;i < npart_tot ;i++)
    {		
        AddParticle(x_[i],y_[i],z_[i],vx_[i],vy_[i],vz_[i],Ion_[i],cloud);
    }


#endif // __MPI_ON__  

    // Set lifetime
    cloud.lifetime = t*1e-3;
    // fill logger file
    if(myid==0)
    {
        plogger<<"Data Imported from ";
        plogger<<ssfilename_input.str();
        plogger<<"\n";
        plogger<<npart_tot;
        plogger<< " particles : " ;
        for(unsigned int j=0;j<frac.size();j++)
        {
            plogger << name_frac[j];
            plogger << " " ;
            plogger << round((double )frac[j]/npart_tot*100.);
            plogger << "% " ;
        }
        plogger << "\n";
    }

    if(myid==0)
        remove("tmp_file_name.txt");
#ifdef  __MPI_ON__ 	
    MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__  

    return npart_tot;
}

#ifdef __MPI_ON__
void MPI_Concatenate_Data(int NRPARTICLES, char * filename_prefix) {
    // MPI
    int myid = MPI::COMM_WORLD.Get_rank();
    int numprocs = MPI::COMM_WORLD.Get_size();

    MPI_Datatype strtype;
    MPI_Type_contiguous(256,MPI::CHAR,&strtype);
    MPI_Type_commit(&strtype);
    // FILE 
    stringstream ssfilename;
    string filename_temp;
    ifstream input_file;
    ofstream output_file; 
    stringstream strsendbuff;
    stringstream tmpstr;
    streampos input_pos;
    char c ='a';
    char line_end0[256];
    char line_end1[256];
    char line[256];
    vector<string > vec_line;

    int index =-1;
    int lastindex = -2;

    int npart_node;
    int npart_tot;
    int coll_node;
    int coll_tot;

    // node file
    filename_temp = filename_prefix;
    filename_temp.erase(filename_temp.size() - 2, filename_temp.size()-1);
    ssfilename <<  filename_temp << "-" << myid << "_pAll.txt" ;
    input_file.open(ssfilename.str().c_str(),ios::in);
    // merged file
    ssfilename.str("");

    ssfilename <<  filename_temp << "_pAll.txt";
    remove(ssfilename.str().c_str());
    MPI::COMM_WORLD.Barrier();
    output_file.open(ssfilename.str().c_str(),ios::app);	
    int dum;
    dum=0;
    int dum2 [numprocs];
    char * buf = NULL;
    //Read
    if(!input_file)
    {
        cout << "Problem: output node file not found" << endl;
        return;
        //exit(EXIT_FAILURE);
    }
    // fill head of file
    if(myid==0)
    {
        for(int i=0;i<7;i++)
        {
            input_file.getline(line,256);//1
            output_file << line << "\n";
        }	
    }
    else
    {
        for(int i=0;i<7;i++)
        {
            input_file.getline(line,256);//1
        }
    } 
    MPI::COMM_WORLD.Barrier();


    while( c!='#')
    {
        // read block for each timestep
        lastindex = -2;
        index =-1;
        //cout << "turn" << endl;
        while((lastindex<index))
        {
            c = input_file.peek();
            //cout << "peek \t "<< c << endl;
            if(c=='#')
            {
                //cout << "end of data" << endl;
                break;
            }	
            input_pos = input_file.tellg(); // get pointer
            lastindex =index;
            input_file >> index; // read index	
            tmpstr.str("");
            //tmpstr << index + myid * NRPARTICLES/numprocs;
            //cout << index << endl;
            input_file.seekg(input_pos); // set point at the beginning of the line
            input_file.getline(line,256);// read the line
            tmpstr <<  line;
            vec_line.push_back(tmpstr.str()); // store the line	
            //cout << lastindex << " " << index << endl;
        }
        // as the first line of the next block has been read -> delete the line in the vector et replace the pointer
        // if it is not the last line
        if(c!='#')
        {
            vec_line.pop_back();
            input_file.seekg(input_pos);
        }
        MPI::COMM_WORLD.Barrier();				

        strsendbuff.str(""); // clean buffer
        for(unsigned int i=0;i<vec_line.size();i++)
        {
            strsendbuff  <<vec_line[i] << "\n"; // fill buffer
        }
        //#pragma unroll numprocs
        MPI::COMM_WORLD.Barrier();
        dum = strsendbuff.str().size()*sizeof(char);
        MPI::COMM_WORLD.Gather(&dum,1,MPI::INT,dum2,1,MPI::INT,0);
        /*
           for(int k=0;k<numprocs;k++)
           {
        //printf("%d %d\n",k,dum2[k]);
        if(myid==k)
        {
        printf("%d %d\n",k,strsendbuff.str().size());
        //printf("%d %s\n",k,strsendbuff.str().c_str());
        }
        }
        for(int k=0;k<numprocs;k++)
        {
        //printf("%d %d\n",k,dum2[k]);
        if(myid==0)
        printf("myid %d %d\n",k,dum2[k]);

        }
        */
        MPI::COMM_WORLD.Barrier();
        if(myid==0)
        {
            output_file.seekp(ios::end); // place the pointer at the end
            output_file << strsendbuff.str(); // fill the output file
        }
        for(int k=1;k<numprocs;k++)
        {       
            MPI::COMM_WORLD.Barrier();
            //  			if(myid==k)
            //  			{
            //  				output_file.seekp(ios::end); // place the pointer at the end
            //  				output_file << strsendbuff.str(); // fill the output file
            //  			}
            //  			MPI::COMM_WORLD.Barrier();
            // 			

            if(myid==k)
            {
                MPI::COMM_WORLD.Send(strsendbuff.str().c_str(),(int) (strsendbuff.str().size()+1)*sizeof(char),MPI::CHAR,0,0);
                //printf("%d send: %d\n",myid,strsendbuff.str().size());
                //printf("%d send:\n%s\n",myid,strsendbuff.str().c_str());
            }
            if(myid==0)
            { 
                delete [] buf;
                buf = new char[ dum2[k]+1];
                MPI::COMM_WORLD.Recv(buf,dum2[k]+1,MPI::CHAR,k,0);
                //printf("%d recv: %d\n",myid,dum2[k]);
                //printf("%d recv:\n%s\n",myid,buf);
                output_file << buf ;
            }
            MPI::COMM_WORLD.Barrier();
        }
        // clean vector
        vec_line.clear();
        MPI::COMM_WORLD.Barrier();
    }	



    if(myid==0)
    {
        input_file.getline(line_end0,256);// read the line
        //output_file << line << "\n";
        input_file.getline(line_end1,256);// read the line
        //output_file << line << "\n";
    }
    else
    {
        input_file.getline(line,256);// read the line
        input_file.getline(line,256);// read the line
    }


    input_file.ignore(256,'=');
    input_file >> npart_node;
    MPI::COMM_WORLD.Reduce(&npart_node,&npart_tot,1,MPI::INT,MPI::SUM,0);
    //	if(myid==0)
    //	{
    //		output_file << "# #Particles left= " << npart_tot << "\n"; 
    //	}
    input_file.ignore(256,'\n');	
    input_file.ignore(256,'=');
    input_file >> coll_node;
    MPI::COMM_WORLD.Reduce(&coll_node,&coll_tot,1,MPI::INT,MPI::SUM,0);
    if(myid!=0)
    {
        output_file.close(); // close file
        input_file.close();	
    }
    MPI::COMM_WORLD.Barrier();
    if(myid==0)
    {
        output_file << line_end0 << "\n";
        output_file << line_end1 << "\n";
        output_file << "# #Particles left= " << npart_tot << "\n"; 
        output_file << "#total number of collisions of all particles = " << coll_tot << "\n";
        input_file.ignore(256,'\n');
        input_file.getline(line,256);// read the line
        output_file << line;
        output_file.close(); // close file
        input_file.close();
    }	

    MPI::COMM_WORLD.Barrier();
    if(myid==0)
    {
        for(int k=0;k<numprocs;k++)
        {
            ssfilename.str("");
            ssfilename <<  filename_temp << "-" << k << "_pAll.txt" ;
            remove(ssfilename.str().c_str());
        }
    }
}

void Rename_part_file(const char* filename_prefix,int nrparticles) {
    // MPI
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


    int error = 0;
    stringstream ss_filename_old;
    stringstream ss_filename_new;
    string filename_temp = filename_prefix;
    filename_temp.erase(filename_temp.size() - 2, filename_temp.size()-1); // remove the "-myid"
    for(int i=0;i<counts[myid];i++)
    {
        ss_filename_old.str("");
        ss_filename_new.str("");
        ss_filename_old << filename_prefix << "_p" << i+1 << ".txt";

        ss_filename_new <<filename_temp << "_p" << 1+i + displ[myid] << ".txt";
        cout <<  ss_filename_old.str() << " " << ss_filename_new.str()<< endl;

        error += rename(ss_filename_old.str().c_str(),ss_filename_new.str().c_str());
    }
    delete [] counts;
    delete [] displ;
    return;
}
#endif // __MPI_ON__


void BundleData(const char* filenamebegin, int nrparticles) {
    //variables 
    int myid =0;
#ifdef __MPI_ON__
    // MPI
    myid = MPI::COMM_WORLD.Get_rank();
    //  	int numprocs = MPI::COMM_WORLD.Get_size();  
    Rename_part_file(filenamebegin, nrparticles);
    MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__
    char line[516];
    //char * tmp;
    char c[nrparticles];
    /*read variables*/
    int index;
    string name;
    double x,y,z; //in (mm) 
    double vx,vy,vz; 
    double R, rplus,rmin,Energy,Temp;  //mm
    double t; //(ms)
    ifstream infiles[nrparticles];
    ofstream output;
    stringstream op;

    char filename[50];  
    bool isopen;
    int nropen=0,nrclosed=0;
    string filename_temp = filenamebegin;
#ifdef __MPI_ON__
    filename_temp.erase(filename_temp.size() - 2, filename_temp.size()-1); // remove the "-0"
    //cout << filename_temp << endl;
#endif // __MPI_ON__    
    if(myid==0)
    {    

        op<<filename_temp<<"_bundled.txt";
        op>>filename;
        output.open(filename);
        output<<"alive\t index\t name\t x\t y\t z\t (mm) vx \t vy \t vz (m/s)\t R \t t\n";

        /*create streams, set stream at first line.*/
        for(int i=0; i<nrparticles; i++){
            stringstream ss;// (stringstream::in | stringstream::out);
            ss<<filename_temp<<"_p"<<i+1<<".txt"; 
            char filenam[50];//different from _filenam (parameter) he!!
            ss >>filenam;

            infiles[i].open(filenam);
            if(!infiles[i]){
                cout << "can't open file " << filenam << endl;infiles[i].close();
                exit(EXIT_FAILURE);
            }else{ 
                //cout << filenam <<" opened.\n";

            }           
        }   


        /*read all the data...*/
        for(int i=0; i<nrparticles; i++){ 
            //read all the comments		 
            while(infiles[i].peek() == '#'){
                infiles[i].getline(line, 516);
                //cout<<"comment found: "<<line<<" \n";
            }	
            int streamposition; 
            //continue reading the stream untill at the end of the file            	
            while(infiles[i].peek() != EOF && infiles[i].peek() != '#'){ 	
                streamposition=infiles[i].tellg();
                infiles[i].getline(line,256);
            }
            //set stream back to the last line and read this line
            infiles[i].seekg(streamposition); 
            //index name x y z (mm) vx vy vz (m/s)     r+(mm)          r-(mm)          R(mm)   Energy(eV)  Temp(K)   t(ms)	
            infiles[i]>>index>>name>>x>>y>>z>>vx>>vy>>vz>>rplus>>rmin>>R>>Energy>>Temp>>t;
            //cout<<"particle read: "<<"\t"<<i<<"\t"<<name<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<R<<"\t"<<t<<endl;
            if(!infiles[i].eof()){		
                //check if the particle is still alive
                infiles[i].getline(line,256); //this is the ******************** line
                infiles[i].getline(line,256);infiles[i].getline(line,256);
                char str2[] = "lost";
                if(strstr(line, str2) == NULL){
                    //particle alive when program ended           
                    output<<1;
                }else{
                    //particle lost when program ended 
                    output<<0;                                           
                }
                output<<"\t"<<i<<"\t"<<name<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t";
                output<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<R<<"\t"<<t<<endl;
                infiles[i].close();       
            }else{
                output<<1<<"\t"<<i<<"\t"<<name<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<R<<"\t"<<t<<endl;
            }
        }   


        /*close streams*/
        for(int i=0; i<nrparticles; i++){
            if(infiles[i]) infiles[i].close(); 
        } 
        output.close();
        cout<<"files succesfully bundled in "<<filename<<endl;
    } //myid==0
#ifdef  __MPI_ON__ 	
    MPI::COMM_WORLD.Barrier();
#endif // __MPI_ON__  
    //printf("%c[%d;%dmProgam Ended%c[%dm\n",27,4,31,27,0);    
}

void PrintIonCloudinfo(const char *beginstream, int nrparticles, double diaphragm_radius_mm) {
    /*read variables*/
    ifstream infiles[nrparticles];
    int index;
    string name; 
    double x,y,z; //in (mm) 	
    double vx,vy,vz,vxy,vxy2; 
    char star; //file: ** 	 
    double rplus,rmin;
    double R;  //mm
    double t; //(ms)
    double Energy;
    double Temp;
    bool CheckFilesOpen = true;
    //other variables
    char line[256]; 
    char c[nrparticles];   
    for(int i=0; i<nrparticles;i++){c[i]='0';}

    //IonCloud variables
    double x_min, x_max, x_mean = 0.0;
    double y_min, y_max, y_mean = 0.0;
    double z_min, z_max, z_mean = 0.0;
    double r_min, r_max, r_mean = 0.0;
    double rp_min, rp_max, rp_mean = 0.0;
    double rm_min, rm_max, rm_mean = 0.0;
    double Energy_min, Energy_max, Energy_mean = 0.0;
    double Temp_mean = 0.0;
    int nr_ions_throughpumping_dia = 0; //Radius of the pumping diaphraghm.
    bool FirstRead = true;
    int nr_particles_alive = 0;


    stringstream ssioncl;// (stringstream::in | stringstream::out);
    ssioncl<<beginstream<<"_cloudEvolution.txt"; 
    char filenameioncl[556];//different from _filename (parameter) he!!
    ssioncl >>filenameioncl;
    ofstream ioncloudinfo;

    ioncloudinfo.open(filenameioncl);
    ioncloudinfo<<"#->1.time(ms) \t 2.#particles \t 3.x_min x_max x_mean (mm) \t 6.y_min, y_max, y_mean (mm)\t 9.z_min, z_max, z_mean (mm)\t 12.r_min, r_max, r_mean (mm)\t 15.rp_min, rp_max, rp_mean (mm)\t 18.rm_min, rm_max, rm_mean (mm)\t 21.Energy_min Energy_max Energy_mean (eV)\t 24.Temp_mean (K)\t 25. %ions with R<PUMPING_DIA \n";




    /*create streams, set stream at first line.*/
    for(int i=0; i<nrparticles; i++){
        stringstream ss;// (stringstream::in | stringstream::out);
        ss<<beginstream<<"_p"<<i+1<<".txt"; 
        char filename[556];//different from _filename (parameter) he!!
        ss >>filename;

        infiles[i].open(filename);
        if(!infiles[i]){
            cout << "can't open file " << filename << endl;
            cin.get();
        }else{ 
            //cout << filenam <<" opened.\n";
            while(infiles[i].peek() == '#'){
                infiles[i].getline(line, 516);
                //cout<<"comment found: "<<line<<" \n";
            }		
        }           
    }


    /*read all the data...*/
    //Read the file he 
    int i;
    int nropen=0,nrclosed=0;
    for(int i=0; i<nrparticles; i++){
        c[i]='a';
        if(infiles[i].is_open()) {nropen++;}
        else nrclosed++;
    } 
    if(nropen == 0) CheckFilesOpen = false;
    else CheckFilesOpen = true;
    while(nrclosed < nrparticles){
        FirstRead= true;
        for(i=0; i<nrparticles; i++){  
            if(!infiles[i].eof() && infiles[i].peek() != '#' && infiles[i].good()){ 
                //index name x y z (mm) vx vy vz (m/s)     r+(mm)          r-(mm)          R(mm)   Energy(eV)     t(ms)	
                infiles[i]>>index>>name>>x>>y>>z>>vx>>vy>>vz>>rplus>>rmin>>R>>Energy>>Temp>>t;infiles[i].getline(line,256);

                if(FirstRead){				    
                    x_min= x_max= x_mean = x;
                    y_min= y_max= y_mean = y;
                    z_min= z_max= z_mean = z;
                    r_min= r_max= r_mean = R;
                    rp_min= rp_max= rp_mean = rplus;
                    rm_min= rm_max= rm_mean = rmin;
                    Energy_min=Energy_max=Energy_mean=Energy; Temp_mean = Temp;
                    nr_particles_alive=1;
                    nr_ions_throughpumping_dia=0;
                    FirstRead = false;
                }else{
                    if(x < x_min) x_min = x;
                    if(x > x_max) x_max = x;
                    x_mean += x;
                    if(y < y_min) y_min = y;
                    if(y > y_max) y_max = y;
                    y_mean += y;
                    if(z < z_min) z_min = z;
                    if(z > z_max) z_max = z;
                    z_mean += z;
                    if(R < r_min) r_min = R;
                    if(R > r_max) r_max = R;
                    r_mean += R;
                    if(rplus < rp_min) rp_min = rplus;
                    if(rplus > rp_max) rp_max = rplus;
                    rp_mean += rplus;
                    if(rmin < rm_min) rm_min = rmin;
                    if(rmin > rm_max) rm_max = rmin;
                    rm_mean += rmin;
                    if(Energy < Energy_min) Energy_min = Energy;
                    if(Energy > Energy_max) Energy_max = Energy;
                    Energy_mean += Energy;
                    Temp_mean += Temp;
                    nr_particles_alive++;
                }//end if(Firstread)					
                if(R < diaphragm_radius_mm){nr_ions_throughpumping_dia++;}
            }else{
                if(infiles[i].is_open()){ 
                    infiles[i].close();nrclosed++;
                    //cout<<"file number "<<i<<" closed with index "<<index<<" nrclosed:"<<nrclosed<<" close because peek is: "<<infiles[i].peek()<< endl;
                }
            }
        }//End of loop over particles


        /*all particles read for one timestep. */
        //calculate means
        double nr_particles_alive_inv=1.0/nr_particles_alive;
        x_mean = x_mean * nr_particles_alive_inv;
        y_mean = y_mean * nr_particles_alive_inv;
        z_mean = z_mean * nr_particles_alive_inv;
        r_mean = r_mean * nr_particles_alive_inv;
        rp_mean = rp_mean * nr_particles_alive_inv;
        rm_mean = rm_mean * nr_particles_alive_inv;
        Energy_mean = Energy_mean * nr_particles_alive_inv;
        Temp_mean = Temp_mean * nr_particles_alive_inv;

        //nr_ions_throughpumping_dia = nr_ions_throughpumping_dia;


        // "1.time(ms) \t 2.#particles \t 3.x_min x_max x_mean \t 6.y_min, y_max, y_mean \t 9.z_min, z_max, z_mean \t 12.r_min, r_max, r_mean \t 15.rp_min, rp_max, rp_mean 
        //	    \t 18.rm_min, rm_max, rm_mean \t 21.Energy_min Energy_max Energy_mean \t 24.Temp_mean \n";

        ioncloudinfo<<t<<"\t"<<nr_particles_alive<<"\t"<<x_min<<" "<<x_max<<" "<<x_mean<<"\t"<<y_min<<" "<< y_max<<" "<< y_mean<<" \t "<<z_min<<" "<< z_max<<" "<< z_mean<<" \t "
            <<r_min<<" "<< r_max<<" "<< r_mean<<" \t "<<rp_min<<" "<< rp_max<<" "<< rp_mean<<" \t "
            <<rm_min<<" "<< rm_max<<" "<< rm_mean<<" \t "<<Energy_min<<" "<<Energy_max<<" "<<Energy_mean<<" \t "<<Temp_mean<<" \t "<<nr_ions_throughpumping_dia<<" \n";


    }//End of loop over while(!Fileopen)



    /*close streams*/
    for(int i=0; i<nrparticles; i++){
        if(infiles[i].is_open()) infiles[i].close(); 
    } 

    ioncloudinfo.close();

    cout<<"cloud evolution information stored in "<<beginstream<<"_cloudEvolution.txt \n";
}

void PrintIonCloudGaussEvo(const char *beginstream, int nrparticles, double diaphragm_radius_mm, bool skipfirstparticle) {
    // function prints out the gaussian distr evolution of x,y,z vx,vy,vz ,rp, rmin(and E, T (why not..)
    /*read variables*/
    unsigned int Noffset = 0; 
    if(skipfirstparticle) nrparticles--;
    if(skipfirstparticle) {Noffset = 1;cout<<"WARNING "<<"skipping the first particle in calculating the cloud evolution\n";}
    else Noffset = 0;


    ifstream infiles[nrparticles];
    int index;
    string name; 
    double x[nrparticles],y[nrparticles],z[nrparticles]; //in (mm) 	
    double vx[nrparticles],vy[nrparticles],vz[nrparticles];
    char star; //file: ** 	 
    double rplus[nrparticles],rmin[nrparticles];
    double R[nrparticles];  //mm
    double t; //(ms)
    double Energy[nrparticles];
    double Temp[nrparticles];
    bool CheckFilesOpen = true;
    //other variables
    char line[256]; 
    char c[nrparticles];   
    for(int i=Noffset; i<nrparticles;i++){
        c[i]='0';
        x[i]=y[i]=z[i]=vx[i]=vy[i]=vz[i]=0;}

        //IonCloud vars
        double x_min, x_max = 0.0;
        double y_min, y_max = 0.0;
        double z_min, z_max = 0.0;
        double r_min, r_max = 0.0;
        double rp_min, rp_max = 0.0;
        double rm_min, rm_max = 0.0;
        double Energy_min, Energy_max = 0.0;
        int nr_ions_throughpumping_dia = 0; //Radius of the pumping diaphraghm.
        bool FirstRead = true;
        int nr_particles_alive = 0;
        //Gaussian vars
        double xu,yu,zu,vxu,vyu,vzu,rpu,rmu,ru,Eu,Tu;
        long double xs,ys,zs,vxs,vys,vzs,rps,rms,rs,Es,Ts;
        xu=yu=zu=vxu=vyu=vzu=rpu=rmu=ru=Eu=Tu=0.0;
        xs=ys=zs=vxs=vys=vzs=rps=rms=rs=Es=Ts=0.0;

        //file vars
        stringstream ssioncl;// (stringstream::in | stringstream::out);
        ssioncl<<beginstream<<"_cloud_GaussEvo.txt"; 
        char filenameioncl[556];//different from _filename (parameter) he!!
        ssioncl >>filenameioncl;
        ofstream ioncloudinfo;

        ioncloudinfo.open(filenameioncl);
        ioncloudinfo<<"->1.time(ms) \t 2.#particles \t xu \t yu \t zu \t 6.vxu \t vyu \t vzu \t \t 9.xs \t ys \t zs \t 12.vxs \t vys \t vzs"<<"\t 15.ru \t  rmu \t rpu \t 18.rs \t rms \t rps \t 21.Eu \t Tu \t Es \t Ts \n";

        ofstream tmp;
        tmp.open("tmp.txt");


        /*create streams, set stream at first line.*/
        for(int i=Noffset; i<nrparticles; i++){
            stringstream ss;// (stringstream::in | stringstream::out);
            ss<<beginstream<<"_p"<<i+1<<".txt"; 
            char filename[556];//different from _filename (parameter) he!!
            ss >>filename;

            infiles[i].open(filename);
            if(!infiles[i]){
                cout << "can't open file " << filename << endl;
                cin.get();
            }else{ 
                //cout << filename <<" opened.\n";
                while(infiles[i].peek() == '#'){
                    infiles[i].getline(line, 516);
                    //cout<<"comment found: "<<line<<" \n";
                }		
            }           
        }


        /*read all the data...*/
        //Read the file he 
        int i;
        int nropen=0,nrclosed=0;
        for(int i=Noffset; i<nrparticles; i++){
            c[i]='a';
            if(infiles[i].is_open()) {nropen++;}
            else nrclosed++;
        } 
        if(nropen == 0) CheckFilesOpen = false;
        else CheckFilesOpen = true;

        while(nrclosed < nrparticles){
            FirstRead= true;
            for(i=Noffset; i<nrparticles; i++){  
                if(!infiles[i].eof() && infiles[i].peek() != '#' && infiles[i].good()){ 	
                    //index name x y z (mm) vx vy vz (m/s)     r+(mm)          r-(mm)          R(mm)   Energy(eV)     t(ms)	
                    infiles[i]>>index>>name>>x[i]>>y[i]>>z[i]>>vx[i]>>vy[i]>>vz[i]>>rplus[i]>>rmin[i]>>R[i]>>Energy[i]>>Temp[i]>>t;infiles[i].getline(line,256);
                    if(FirstRead){				    
                        x_min= x_max=x[i];
                        y_min= y_max=y[i];
                        z_min= z_max=z[i];
                        r_min= r_max=R[i];
                        rp_min= rp_max= rpu = rplus[i];
                        rm_min= rm_max= rmu = rmin[i];
                        Energy_min=Energy_max=Energy[i]; 
                        nr_particles_alive=1;
                        nr_ions_throughpumping_dia=0;
                        FirstRead = false;
                    }else{
                        if(x[i] < x_min) x_min = x[i];
                        if(x[i] > x_max) x_max = x[i];
                        if(y[i] < y_min) y_min = y[i];
                        if(y[i] > y_max) y_max = y[i];
                        if(z[i] < z_min) z_min = z[i];
                        if(z[i] > z_max) z_max = z[i];
                        if(R[i] < r_min) r_min = R[i];
                        if(R[i] > r_max) r_max = R[i];
                        if(rplus[i] < rp_min) rp_min = rplus[i];
                        if(rplus[i] > rp_max) rp_max = rplus[i];
                        if(rmin[i] < rm_min) rm_min = rmin[i];
                        if(rmin[i] > rm_max) rm_max = rmin[i];
                        if(Energy[i] < Energy_min) Energy_min = Energy[i];
                        if(Energy[i] > Energy_max) Energy_max = Energy[i];
                        nr_particles_alive++;
                    }//end if(Firstread)					
                    if(R[i] < diaphragm_radius_mm){nr_ions_throughpumping_dia++;}
                }else{
                    if(infiles[i].is_open()){ 
                        infiles[i].close();nrclosed++;
                        //cout<<"file number "<<i<<" closed with index "<<index<<" nrclosed:"<<nrclosed<<" close because peek is: "<<infiles[i].peek()<< endl;
                    }

                }
            }//End of loop over particles


            /*all particles read for one timestep. */
            //calculate means
            double nr_particles_alive_inv=1.0/nr_particles_alive;
            //init values again
            xu=yu=zu=vxu=vyu=vzu=rpu=rmu=ru=Eu=Tu=0.0;
            xs=ys=zs=vxs=vys=vzs=rps=rms=rs=Es=Ts=0.0;

            double tt= 0 ;
            //calculate all averages: DONE FINE!
            for(int i=Noffset ; i<nrparticles; i++){
                //multiply with inv of particles before otherwise it runs away...
                xu +=x[i];yu +=y[i];zu +=z[i];
                vxu += vx[i];vyu +=vy[i]; vzu +=vz[i];
                ru +=R[i];
                rpu +=rplus[i];
                rmu +=rmin[i];
                Eu +=Energy[i];
                Tu +=Temp[i];
            }
            xu = xu * nr_particles_alive_inv;
            yu = yu * nr_particles_alive_inv;
            zu = zu * nr_particles_alive_inv;
            vxu = vxu * nr_particles_alive_inv;
            vyu = vyu * nr_particles_alive_inv;
            vzu = vzu * nr_particles_alive_inv;
            ru = ru * nr_particles_alive_inv;
            rpu = rpu * nr_particles_alive_inv;
            rmu = rmu * nr_particles_alive_inv;
            Eu = Eu * nr_particles_alive_inv;
            Tu = Tu * nr_particles_alive_inv;

            /*finds the sample variance
             *=(sum[(x-mean(x))^2]/n)
             */
            for(int i=Noffset ; i<nrparticles -1 ; i++){
                //multiply with inv of particles before otherwise it runs away...
                xs +=pow((x[i]-xu),2);
                ys +=pow((y[i]-yu),2);
                zs +=pow((z[i]-zu),2);
                vxs +=pow((vx[i]-vxu),2);
                vys +=pow((vy[i]-vyu),2);
                vzs +=pow((vz[i]-vzu),2);
                rs +=pow((R[i]-ru),2);
                rps +=pow((rplus[i]-rpu),2);
                rms +=pow((rmin[i]-rmu),2);
                Es +=pow((Energy[i]-Eu),2);
                Ts +=pow((Temp[i]-Tu),2);
            }
            xs =sqrt(xs*nr_particles_alive_inv);
            ys =sqrt(ys*nr_particles_alive_inv);
            zs =sqrt(zs*nr_particles_alive_inv);
            vxs =sqrt(vxs*nr_particles_alive_inv);
            vys =sqrt(vys*nr_particles_alive_inv);
            vzs =sqrt(vzs*nr_particles_alive_inv);
            rs =sqrt(rs*nr_particles_alive_inv);
            rps =sqrt(rps*nr_particles_alive_inv);
            rms =sqrt(rms*nr_particles_alive_inv);
            Es =sqrt(Es*nr_particles_alive_inv);
            Ts =sqrt(Ts*nr_particles_alive_inv);

            //cout<<"time: "<<t <<"\t average stdev z:"<<zu<<"\t"<<zs<<" min  max:" <<z_min << "\t" <<z_max<<endl;
            //cout<<"time: "<<t <<"avg "<<zu<<"\t stdev "<<zs<<endl;
            //print like incloudevolution with the u calculated differently here for comparison
            // "1.time(ms) \t 2.#particles \t 3.x_min x_max x_mean \t 6.y_min, y_max, y_mean \t 9.z_min, z_max, z_mean \t 12.r_min, r_max, r_mean \t 15.rp_min, rp_max, rp_mean 
            //	    \t 18.rm_min, rm_max, rm_mean \t 21.Energy_min Energy_max Energy_mean \t 24.Temp_mean \n";
            //ioncloudinfo<<t<<"\t"<<nr_particles_alive<<"\t"<<x_min<<" "<<x_max<<" "<<xu<<"\t"<<y_min<<" "<< y_max<<" "<< yu<<" \t "<<z_min<<" "<< z_max<<" "<< zu<<" \t " 		
            // <<r_min<<" "<< r_max<<" "<< ru<<" \t "<<rp_min<<" "<< rp_max<<" "<< rpu<<" \t "<<rm_min<<" "
            // << rm_max<<" "<< rmu<<" \t "<<Energy_min<<" "<<Energy_max<<" "<<Eu<<" \t "<<Tu<<" \t "<<nr_ions_throughpumping_dia<<" \n";
            //ioncloudinfo<<"->1.time(ms) \t 2.#particles \t xu \t yu \t zu \t 6.vxu \t vyu \t vzu \t \t 9.xs \t ys \t zs \t 12.vxs \t vys \t vzs \n";
            ioncloudinfo<<t<<"\t"<<nr_particles_alive<<"\t"<<xu<<"\t"<< yu<<"\t"<< zu<<"\t"<< vxu<<"\t"<< vyu<<"\t"<< vzu<<"\t"<< xs<<"\t"<< ys<<"\t"<<zs <<"\t"<< vxs <<"\t"<< vys<<"\t"<< vzs
                <<"\t"<<ru<<"\t"<<rmu<<"\t"<<rpu<<"\t"<<rs<<"\t"<<rms<<"\t"<<rps<<"\t"<<Eu<<"\t"<<Tu<<"\t"<<Es<<"\t"<<Ts
                <<"\n";

        }//while



        /*close streams*/
        for(int i=Noffset ; i<nrparticles; i++){
            if(infiles[i].is_open()) infiles[i].close(); 
        } 

        ioncloudinfo.close();

        cout<<"Gaussian cloud evolution information stored in "<<filenameioncl<<endl;
}
