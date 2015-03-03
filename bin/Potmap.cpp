#include "Potmap.h"
#include "MPI_simbuca.h"

using namespace std;

double abs(double _x)
{
    if( _x >0 ) return _x;
    else return -_x;
}
//______________________________________________________________________________
PotMap::PotMap()
{
    rmin = 0;
    rmax = 0;
    zmin = 1e8;
    zmax = -1e8;
    verbose = false;
    //interpolator
    ap.resize(9);
    fp.resize(9);
}
//______________________________________________________________________________

//______________________________________________________________________________
PotMap::~PotMap()
{

}
//______________________________________________________________________________

//______________________________________________________________________________
void PotMap::ReadEPotential(string Potmap_name)
{
    ifstream Potmap_file;

    Potmap_file.open(Potmap_name.c_str(),ios::in);
    if (!Potmap_file){
        SLogger slogger("trapparameters");
        slogger << ERROR << "Can't open "<<Potmap_name << SLogger::endmsg;exit(1);
    }

    double r,z,pot;


    bool compute_dr = true;
    bool second_line = true;
    int compt_line = 0;
    char firstline[256];
    char comment= '%';
    while (Potmap_file.peek() == comment){
        Potmap_file.getline(firstline, 256);
    }



    //first line
    Potmap_file>>r>>z>>pot;
    compt_line++;
    rp.push_back(r);
    zp.push_back(z);
    potential.push_back(pot);
    if(r<rmin){rmin=r;}if(r>rmax){rmax=r;}if(z<zmin){zmin=z;} if(z>zmax){zmax=z;}
    while ( !Potmap_file.eof())
    {
        // keep reading until end-of-file
        Potmap_file>>r>>z>>pot;
        compt_line++;


        if(second_line)
        {
            dz =  z -zp.back();
            second_line = false;
        }
        if(compute_dr&&(rp.back()!=r))
        {
            dr = r- rp.back();
            numz =compt_line -1;
            compute_dr = false;
        }
        rp.push_back(r);
        zp.push_back(z);
        potential.push_back(pot);
        if(r<rmin){rmin=r;}if(r>rmax){rmax=r;}if(z<zmin){zmin=z;} if(z>zmax){zmax=z;}

    }
    rp.pop_back();
    zp.pop_back();
    potential.pop_back();
    Potmap_file.close();
    compt_line--;

    numr = compt_line/numz;



    if(verbose)
    {
        cout << "number of line " << compt_line << endl;
        cout << "dr " << dr << endl;
        cout << "dz " << dz << endl;
        cout << "numz " << numz << endl;
        cout << "dz check " << (zmax-zmin)/(numz-1) << endl;
        cout << "numr " << numr << endl;
        cout << "dr check " << (rmax-rmin)/(numr-1) << endl;
    }

    if(abs(dr - (rmax-rmin)/(numr-1)) > 1e-9)
    {
        cout << "Potential map not regular: check dr step " <<endl;
        exit(-1);
    }

    if(abs(dz - ((zmax-zmin)/(numz-1))) > 1e-9)
    {

        cout << dz*(numz-1) << endl;
        cout << (zmax-zmin) << endl;
        cout << "Potential map not regular: check dz step " <<endl;
        exit(-1);
    }
    int myid=0;
#ifdef __MPI_ON__
    myid = MPI::COMM_WORLD.Get_rank();
#endif // _MPI_ON__
    if(myid==0)
    {
        cout<<"\t R: "<<rmin<<" "<<rmax<<" "<<dr<<endl;
        cout<<"\t z: "<<zmin<<" "<<zmax<<" "<<dz<<endl;

        cout<<"file: "<<Potmap_name<<" read\n";
    }
    //cout << "end" << endl;
    maxIndex = (numz-1)*(numr-1);
    return;
}
//______________________________________________________________________________

//______________________________________________________________________________
int PotMap::GetIndex(const double &r,const double &z)
{
    if(r<rmin||r>rmax||z<zmin||z>zmax)
    {
        cout << "Particle out of the mesh " << endl;
        cout << 0 << " " << r << " " << rmax << endl;
        cout << zmin << " " << z << " " << zmax<< endl;
        exit(-1);
    }
    int i,j;
    double r_index= r-rmin;
    double z_index= z-zmin;


    i=fabs(floor(r_index/dr));
    j=floor(z_index/dz);
    //cout << "i , j " << i << " " << j << endl;
    while(i>numr-3)
    {
        i--;
    }
    while(j>numz-3)
    {
        j--;
    }

    //cout << "i , j " << i << " " << j << endl;
    return i*numz+j;


}
//______________________________________________________________________________


//______________________________________________________________________________
double PotMap::GetPotential(const double &x,const double & y, const double & z)
{
    //z ,r
    double r = sqrt(x*x+y*y);

    double res;
    int index00 = GetIndex(r,z);
    int index10 = index00 + 1;
    int index20 = index00 + 2;
    int index01 = index00 + numz;
    int index11 = index00 + numz + 1;
    int index21 = index00 + numz + 2;
    int index02 = index00 + 2*numz;
    int index12 = index00 + 2*numz + 1;
    int index22 = index00 + 2*numz + 2;

    /*
       cout << index00 << endl;
       cout << index10 << endl;
       cout << index20 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       */



    fp[0] = potential[index00];//f11
    fp[1] = potential[index10];//f21
    fp[2] = potential[index20];//f31
    fp[3] = potential[index01];//f12
    fp[4] = potential[index11];//f22
    fp[5] = potential[index21];//f32
    fp[6] = potential[index02];//f13
    fp[7] = potential[index12];//f23
    fp[8] = potential[index22];//f33


    /*
       for(int k=0;k<9;k++)
       cout << k<<" " << fp[k] << endl;

*/


    double x0x1 = dz;
    double x0x2 = 2*dz;
    double x1x2 = dz;
    double y0y1 = dr;
    double y0y2 = 2*dr;
    double y1y2 = dr;

    double facx1 = x0x1*x0x2;
    double facx2 = -x0x1*x1x2;
    double facx3 = x0x2*x1x2;
    double facy1 = y0y1*y0y2;
    double facy2 = -y0y1*y1y2;
    double facy3 = y0y2*y1y2;



    ap[0] = fp[0] / (facx1*facy1);
    ap[1] = fp[1] / (facx2*facy1);
    ap[2] = fp[2] / (facx3*facy1);
    ap[3] = fp[3] / (facx1*facy2);
    ap[4] = fp[4] / (facx2*facy2);
    ap[5] = fp[5] / (facx3*facy2);
    ap[6] = fp[6] / (facx1*facy3);
    ap[7] = fp[7] / (facx2*facy3);
    ap[8] = fp[8] / (facx3*facy3);

    double xx1 = (z - zp[index00]);
    double xx2 = (z - zp[index11]);
    double xx3 = (z - zp[index22]);
    double yy1 = (r - rp[index00]);
    double yy2 = (r - rp[index11]);
    double yy3 = (r - rp[index22]);

    facx1 = xx2*xx3;
    facx2 = xx1*xx3;
    facx3 = xx1*xx2;
    facy1 = yy2*yy3;
    facy2 = yy1*yy3;
    facy3 = yy1*yy2;

    res =     ap[0]*facx1*facy1
        +     ap[1]*facx2*facy1
        +     ap[2]*facx3*facy1
        +     ap[3]*facx1*facy2
        +     ap[4]*facx2*facy2
        +     ap[5]*facx3*facy2
        +     ap[6]*facx1*facy3
        +     ap[7]*facx2*facy3
        +     ap[8]*facx3*facy3;

    return res;
}
//______________________________________________________________________________

//______________________________________________________________________________
double PotMap::GetEr(const double &x,const double & y, const double & z)
{
    //z ,r
    double r = sqrt(x*x+y*y);
    if(r==0)
        return 0;
    double res;
    int index00 = GetIndex(r,z);
    int index10 = index00 + 1;
    int index20 = index00 + 2;
    int index01 = index00 + numz;
    int index11 = index00 + numz + 1;
    int index21 = index00 + numz + 2;
    int index02 = index00 + 2*numz;
    int index12 = index00 + 2*numz + 1;
    int index22 = index00 + 2*numz + 2;

    /*
       cout << index00 << endl;
       cout << index10 << endl;
       cout << index20 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       */



    fp[0] = potential[index00];//f11
    fp[1] = potential[index10];//f21
    fp[2] = potential[index20];//f31
    fp[3] = potential[index01];//f12
    fp[4] = potential[index11];//f22
    fp[5] = potential[index21];//f32
    fp[6] = potential[index02];//f13
    fp[7] = potential[index12];//f23
    fp[8] = potential[index22];//f33


    /*
       for(int k=0;k<9;k++)
       cout << k<<" " << fp[k] << endl;

*/


    double x0x1 = dz;
    double x0x2 = 2*dz;
    double x1x2 = dz;
    double y0y1 = dr;
    double y0y2 = 2*dr;
    double y1y2 = dr;

    double facx1 = x0x1*x0x2;
    double facx2 = -x0x1*x1x2;
    double facx3 = x0x2*x1x2;
    double facy1 = y0y1*y0y2;
    double facy2 = -y0y1*y1y2;
    double facy3 = y0y2*y1y2;



    ap[0] = fp[0] / (facx1*facy1);
    ap[1] = fp[1] / (facx2*facy1);
    ap[2] = fp[2] / (facx3*facy1);
    ap[3] = fp[3] / (facx1*facy2);
    ap[4] = fp[4] / (facx2*facy2);
    ap[5] = fp[5] / (facx3*facy2);
    ap[6] = fp[6] / (facx1*facy3);
    ap[7] = fp[7] / (facx2*facy3);
    ap[8] = fp[8] / (facx3*facy3);

    /*
       double xx1 = (z - zp[index00]);
       double xx2 = (z - zp[index11]);
       double xx3 = (z - zp[index22]);
       double yy1 = (r - rp[index00]);
       double yy2 = (r - rp[index11]);
       double yy3 = (r - rp[index22]);

       facx1 = xx2*xx3;
       facx2 = xx1*xx3;
       facx3 = xx1*xx2;
       facy1 = yy2*yy3;
       facy2 = yy1*yy3;
       facy3 = yy1*yy2;
       */
    double xx1 = (z - zp[index00]);
    double xx2 = (z - zp[index11]);
    double xx3 = (z - zp[index22]);

    double dyy1 = (2.*r - rp[index11]  - rp[index22]);
    double dyy2 = (2.*r - rp[index00]  - rp[index22]);
    double dyy3 = (2.*r - rp[index00]  - rp[index11]);




    facx1 = xx2*xx3;
    facx2 = xx1*xx3;
    facx3 = xx1*xx2;

    res =     ap[0]*facx1*dyy1
        +     ap[1]*facx2*dyy1
        +     ap[2]*facx3*dyy1
        +     ap[3]*facx1*dyy2
        +     ap[4]*facx2*dyy2
        +     ap[5]*facx3*dyy2
        +     ap[6]*facx1*dyy3
        +     ap[7]*facx2*dyy3
        +     ap[8]*facx3*dyy3;

    return -res;
}
//______________________________________________________________________________

//______________________________________________________________________________
double PotMap::GetEz(const double &x,const double & y, const double & z)
{
    //z ,r
    double r = sqrt(x*x+y*y);

    double res;
    int index00 = GetIndex(r,z);
    int index10 = index00 + 1;
    int index20 = index00 + 2;
    int index01 = index00 + numz;
    int index11 = index00 + numz + 1;
    int index21 = index00 + numz + 2;
    int index02 = index00 + 2*numz;
    int index12 = index00 + 2*numz + 1;
    int index22 = index00 + 2*numz + 2;

    /*
       cout << index00 << endl;
       cout << index10 << endl;
       cout << index20 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       */



    fp[0] = potential[index00];//f11
    fp[1] = potential[index10];//f21
    fp[2] = potential[index20];//f31
    fp[3] = potential[index01];//f12
    fp[4] = potential[index11];//f22
    fp[5] = potential[index21];//f32
    fp[6] = potential[index02];//f13
    fp[7] = potential[index12];//f23
    fp[8] = potential[index22];//f33


    /*
       for(int k=0;k<9;k++)
       cout << k<<" " << fp[k] << endl;

*/


    double x0x1 = dz;
    double x0x2 = 2*dz;
    double x1x2 = dz;
    double y0y1 = dr;
    double y0y2 = 2*dr;
    double y1y2 = dr;

    double facx1 = x0x1*x0x2;
    double facx2 = -x0x1*x1x2;
    double facx3 = x0x2*x1x2;
    double facy1 = y0y1*y0y2;
    double facy2 = -y0y1*y1y2;
    double facy3 = y0y2*y1y2;



    ap[0] = fp[0] / (facx1*facy1);
    ap[1] = fp[1] / (facx2*facy1);
    ap[2] = fp[2] / (facx3*facy1);
    ap[3] = fp[3] / (facx1*facy2);
    ap[4] = fp[4] / (facx2*facy2);
    ap[5] = fp[5] / (facx3*facy2);
    ap[6] = fp[6] / (facx1*facy3);
    ap[7] = fp[7] / (facx2*facy3);
    ap[8] = fp[8] / (facx3*facy3);

    /*
       double xx1 = (z - zp[index00]);
       double xx2 = (z - zp[index11]);
       double xx3 = (z - zp[index22]);
       double yy1 = (r - rp[index00]);
       double yy2 = (r - rp[index11]);
       double yy3 = (r - rp[index22]);

       facx1 = xx2*xx3;
       facx2 = xx1*xx3;
       facx3 = xx1*xx2;
       facy1 = yy2*yy3;
       facy2 = yy1*yy3;
       facy3 = yy1*yy2;
       */
    double dxx1 = (2.*z - zp[index11]  - zp[index22]);
    double dxx2 = (2.*z - zp[index00]  - zp[index22]);
    double dxx3 = (2.*z - zp[index00]  - zp[index11]);

    double yy1 = (r - rp[index00]);
    double yy2 = (r - rp[index11]);
    double yy3 = (r - rp[index22]);



    facy1 = yy2*yy3;
    facy2 = yy1*yy3;
    facy3 = yy1*yy2;

    res =     ap[0]*dxx1*facy1
        +     ap[1]*dxx2*facy1
        +     ap[2]*dxx3*facy1
        +     ap[3]*dxx1*facy2
        +     ap[4]*dxx2*facy2
        +     ap[5]*dxx3*facy2
        +     ap[6]*dxx1*facy3
        +     ap[7]*dxx2*facy3
        +     ap[8]*dxx3*facy3;

    return -res;
}
//______________________________________________________________________________


//______________________________________________________________________________
pair<double,double > PotMap::GetE(const double &x,const double & y, const double & z)
{
    //z ,r
    double r = sqrt(x*x+y*y);
    pair<double,double > res;
    double Er,Ez;
    int index00 = GetIndex(r,z);
    int index10 = index00 + 1;
    int index20 = index00 + 2;
    int index01 = index00 + numz;
    int index11 = index00 + numz + 1;
    int index21 = index00 + numz + 2;
    int index02 = index00 + 2*numz;
    int index12 = index00 + 2*numz + 1;
    int index22 = index00 + 2*numz + 2;

    /*
       cout << index00 << endl;
       cout << index10 << endl;
       cout << index20 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       */



    fp[0] = potential[index00];//f11
    fp[1] = potential[index10];//f21
    fp[2] = potential[index20];//f31
    fp[3] = potential[index01];//f12
    fp[4] = potential[index11];//f22
    fp[5] = potential[index21];//f32
    fp[6] = potential[index02];//f13
    fp[7] = potential[index12];//f23
    fp[8] = potential[index22];//f33


    /*
       for(int k=0;k<9;k++)
       cout << k<<" " << fp[k] << endl;

*/


    double x0x1 = dz;
    double x0x2 = 2*dz;
    double x1x2 = dz;
    double y0y1 = dr;
    double y0y2 = 2*dr;
    double y1y2 = dr;

    double facx1 = x0x1*x0x2;
    double facx2 = -x0x1*x1x2;
    double facx3 = x0x2*x1x2;
    double facy1 = y0y1*y0y2;
    double facy2 = -y0y1*y1y2;
    double facy3 = y0y2*y1y2;



    ap[0] = fp[0] / (facx1*facy1);
    ap[1] = fp[1] / (facx2*facy1);
    ap[2] = fp[2] / (facx3*facy1);
    ap[3] = fp[3] / (facx1*facy2);
    ap[4] = fp[4] / (facx2*facy2);
    ap[5] = fp[5] / (facx3*facy2);
    ap[6] = fp[6] / (facx1*facy3);
    ap[7] = fp[7] / (facx2*facy3);
    ap[8] = fp[8] / (facx3*facy3);

    /*
       double xx1 = (z - zp[index00]);
       double xx2 = (z - zp[index11]);
       double xx3 = (z - zp[index22]);
       double yy1 = (r - rp[index00]);
       double yy2 = (r - rp[index11]);
       double yy3 = (r - rp[index22]);

       facx1 = xx2*xx3;
       facx2 = xx1*xx3;
       facx3 = xx1*xx2;
       facy1 = yy2*yy3;
       facy2 = yy1*yy3;
       facy3 = yy1*yy2;
       */
    double dxx1 = (2.*z - zp[index11]  - zp[index22]);
    double dxx2 = (2.*z - zp[index00]  - zp[index22]);
    double dxx3 = (2.*z - zp[index00]  - zp[index11]);

    double yy1 = (r - rp[index00]);
    double yy2 = (r - rp[index11]);
    double yy3 = (r - rp[index22]);

    double xx1 = (z - zp[index00]);
    double xx2 = (z - zp[index11]);
    double xx3 = (z - zp[index22]);

    double dyy1 = (2.*r - rp[index11]  - rp[index22]);
    double dyy2 = (2.*r - rp[index00]  - rp[index22]);
    double dyy3 = (2.*r - rp[index00]  - rp[index11]);

    facy1 = yy2*yy3;
    facy2 = yy1*yy3;
    facy3 = yy1*yy2;

    facx1 = xx2*xx3;
    facx2 = xx1*xx3;
    facx3 = xx1*xx2;

    Ez =     ap[0]*dxx1*facy1
        +     ap[1]*dxx2*facy1
        +     ap[2]*dxx3*facy1
        +     ap[3]*dxx1*facy2
        +     ap[4]*dxx2*facy2
        +     ap[5]*dxx3*facy2
        +     ap[6]*dxx1*facy3
        +     ap[7]*dxx2*facy3
        +     ap[8]*dxx3*facy3;



    Er =     ap[0]*facx1*dyy1
        +     ap[1]*facx2*dyy1
        +     ap[2]*facx3*dyy1
        +     ap[3]*facx1*dyy2
        +     ap[4]*facx2*dyy2
        +     ap[5]*facx3*dyy2
        +     ap[6]*facx1*dyy3
        +     ap[7]*facx2*dyy3
        +     ap[8]*facx3*dyy3;

    if(r==0)
        Er = 0;

    res.first = -Er;
    res.second =-Ez;
    return res;


}


//______________________________________________________________________________
vector<double > PotMap::GetEField_from_Pot(const double &x,const double & y, const double & z)
{
    //z ,r
    vector<double> electric_field(3);
    double r = sqrt(x*x+y*y);

    double Er,Ez;
    int index00 = GetIndex(r,z);
    int index10 = index00 + 1;
    int index20 = index00 + 2;
    int index01 = index00 + numz;
    int index11 = index00 + numz + 1;
    int index21 = index00 + numz + 2;
    int index02 = index00 + 2*numz;
    int index12 = index00 + 2*numz + 1;
    int index22 = index00 + 2*numz + 2;

    /*
       cout << index00 << endl;
       cout << index10 << endl;
       cout << index20 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       cout << index01 << endl;
       cout << index11 << endl;
       cout << index21 << endl;
       */



    fp[0] = potential[index00];//f11
    fp[1] = potential[index10];//f21
    fp[2] = potential[index20];//f31
    fp[3] = potential[index01];//f12
    fp[4] = potential[index11];//f22
    fp[5] = potential[index21];//f32
    fp[6] = potential[index02];//f13
    fp[7] = potential[index12];//f23
    fp[8] = potential[index22];//f33


    /*
       for(int k=0;k<9;k++)
       cout << k<<" " << fp[k] << endl;

*/


    double x0x1 = dz;
    double x0x2 = 2*dz;
    double x1x2 = dz;
    double y0y1 = dr;
    double y0y2 = 2*dr;
    double y1y2 = dr;

    double facx1 = x0x1*x0x2;
    double facx2 = -x0x1*x1x2;
    double facx3 = x0x2*x1x2;
    double facy1 = y0y1*y0y2;
    double facy2 = -y0y1*y1y2;
    double facy3 = y0y2*y1y2;



    ap[0] = fp[0] / (facx1*facy1);
    ap[1] = fp[1] / (facx2*facy1);
    ap[2] = fp[2] / (facx3*facy1);
    ap[3] = fp[3] / (facx1*facy2);
    ap[4] = fp[4] / (facx2*facy2);
    ap[5] = fp[5] / (facx3*facy2);
    ap[6] = fp[6] / (facx1*facy3);
    ap[7] = fp[7] / (facx2*facy3);
    ap[8] = fp[8] / (facx3*facy3);

    /*
       double xx1 = (z - zp[index00]);
       double xx2 = (z - zp[index11]);
       double xx3 = (z - zp[index22]);
       double yy1 = (r - rp[index00]);
       double yy2 = (r - rp[index11]);
       double yy3 = (r - rp[index22]);

       facx1 = xx2*xx3;
       facx2 = xx1*xx3;
       facx3 = xx1*xx2;
       facy1 = yy2*yy3;
       facy2 = yy1*yy3;
       facy3 = yy1*yy2;
       */
    double dxx1 = (2.*z - zp[index11]  - zp[index22]);
    double dxx2 = (2.*z - zp[index00]  - zp[index22]);
    double dxx3 = (2.*z - zp[index00]  - zp[index11]);

    double yy1 = (r - rp[index00]);
    double yy2 = (r - rp[index11]);
    double yy3 = (r - rp[index22]);

    double xx1 = (z - zp[index00]);
    double xx2 = (z - zp[index11]);
    double xx3 = (z - zp[index22]);

    double dyy1 = (2.*r - rp[index11]  - rp[index22]);
    double dyy2 = (2.*r - rp[index00]  - rp[index22]);
    double dyy3 = (2.*r - rp[index00]  - rp[index11]);

    facy1 = yy2*yy3;
    facy2 = yy1*yy3;
    facy3 = yy1*yy2;

    facx1 = xx2*xx3;
    facx2 = xx1*xx3;
    facx3 = xx1*xx2;

    Ez =     ap[0]*dxx1*facy1
        +     ap[1]*dxx2*facy1
        +     ap[2]*dxx3*facy1
        +     ap[3]*dxx1*facy2
        +     ap[4]*dxx2*facy2
        +     ap[5]*dxx3*facy2
        +     ap[6]*dxx1*facy3
        +     ap[7]*dxx2*facy3
        +     ap[8]*dxx3*facy3;



    Er =     ap[0]*facx1*dyy1
        +     ap[1]*facx2*dyy1
        +     ap[2]*facx3*dyy1
        +     ap[3]*facx1*dyy2
        +     ap[4]*facx2*dyy2
        +     ap[5]*facx3*dyy2
        +     ap[6]*facx1*dyy3
        +     ap[7]*facx2*dyy3
        +     ap[8]*facx3*dyy3;

    if(r==0)
    {
        electric_field[0] = 0;
        electric_field[1] = 0;
    }
    else
    {
        electric_field[0] = -Er*x/r;
        electric_field[1] = -Er*y/r;
    }
    electric_field[2] = -Ez;
    return electric_field;


}
//______________________________________________________________________________

