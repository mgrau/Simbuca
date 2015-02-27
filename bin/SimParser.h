#ifndef SIMPARSER_H
#define SIMPARSER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "SLogger.h"


#ifdef __GUI_ON__
#include <QTextBrowser>
#endif

using namespace std;
//----------------------------------
struct Operation_sim{
    string name;
    double time;
    double amp;
    double amp1;
    double amp2;
    int order;
    double frequency;
    double freq1;
    double freq2;
    double freq_bias;
    string EigenLett;
    string Element;
    int size;
    string name_file;
};
//----------------------------------
class SimParser {

 public:
  SimParser(const char * simfile);
  ~SimParser();
  void PrintParams();
  void PrintFile();
  void PrintFile(char * filename);

  void SetTrapConfigOff();
#ifdef __GUI_ON__
 void PrintParser(QTextEdit * textwindow);
#endif



  vector < Operation_sim > operation_vec;
 private:
  void ProcessFile();

  

 private:
  ifstream infile;
  vector< pair < string, std::vector< string > > > config_map;
  vector< vector < string > > operation_map;

  std::vector < string > comments;
  std::map<string, const int> label_map;
  std::map<string, const int> label_operation_map;
    enum labels{LCREATECLOUD = 1, LCLOUDPARTS = 2, LCLOUDCOORD = 3, LTEMP = 4, LCREATEPARTICLES = 5, LPARTICLES = 6, LBUFFER = 7, LODE = 8, LCOULOMB = 9, LOUTPUT = 10, LREALTRAP = 11, LNE = 12, LIMPORTDATA = 13, LDEW = 14, LDEB = 15, LQEW = 16, LQEB = 17, LOEW = 18, LOEB = 19, LRWB = 20, LRWW = 21, LAWB = 22, LAWW = 23, LIDEALTRAP = 24, LNONIDEALTRAP = 25, LCALCTRAP = 26, LAC = 27, LSC = 28, LPOTTRAP = 29, LAR=30, LFB=31, LSWEEP=32,LSWIFT=33,LEXC_EMAP=34,LPI_PULSE=35,LASYM_ARW=36};


 public:
  struct CREATECLOUD{
    bool flag;
    int cloudsize;
    int constituents;
    void printme(){cout << "cloudsize: " << cloudsize << ", Nconstituents: " <<  constituents << endl;};
  };
  CREATECLOUD mycreatecloud;
  //----------------------------------
  struct CLOUDPARTS{
    bool flag;
     std::vector<std::pair<string, double> > cloudfracs;
    double totalfrac;
    void printme(){
      cout << "cloudfracs: " << endl;
      for(int i=0; i < cloudfracs.size(); i++)
	cout << "\t" << cloudfracs[i].first << ", " << cloudfracs[i].second << endl;
    } 
  };
  CLOUDPARTS mycloudparts;
  //----------------------------------
  struct CLOUDCOORD{
    bool flag;
    double semiaxis_cloud[3];
    double offset_cloud[3];
    void printme(){
      cout << "cloud coordinates: " << endl;
      cout << "\t " << flush;
      for(int i = 0; i < 3; i++)
	cout << semiaxis_cloud[i] << " " << offset_cloud[i] << " " << flush;
      cout << endl;
    }
  };
  CLOUDCOORD mycloudcoord;
  //----------------------------------
  struct TEMP{
    bool flag;
    vector<double> temp;
    void printme(){
        for(unsigned int i=0;i<temp.size();i++)
            cout << "temp: " << temp[i] << endl;
    }
  };
  TEMP mytemp;
  //----------------------------------
  struct CREATEPARTICLES{
    bool flag;
    int nparts;
    void printme(){
      cout << "create particles: " << nparts << endl;
    }
  };
  CREATEPARTICLES mycreateparticles;
  //----------------------------------
  struct PARTICLES{
    bool flag;
    std::vector<std::pair<string, std::vector<double> > > particles; //a vector with each element being a map of a string with a vector. The string holds the particle name, the vector holds the particle positions.
    void printme(){
      cout << "single particles: "<<endl;
      for(int i=0; i < particles.size(); i++){
	cout << "\t" << particles[i].first << " " << flush;
	vector<double> myvec = particles[i].second;
	for(int i = 0 ; i < myvec.size(); i++)
	  cout << myvec[i] << " " << flush;
	cout<<endl;
      }
      cout << endl;
    }
  };
  PARTICLES myparticles;
  //----------------------------------
  struct BUFFER{
    bool flag;
    int index;
    double pressure;
    void printme(){
      cout << "buffer: " << index << " " << pressure << endl;
    }
  };
  BUFFER mybuffer;
  //----------------------------------
  struct ODE{
    bool flag;
    int order;
    int adaptive;// 0 - no, 1 - yes
    double stepsize;
    void printme(){
      cout << "ODE: " << order << " " << adaptive << " " << stepsize << endl;
    }
  };
  ODE myode;
  //----------------------------------
  struct COULOMB{
    bool flag;
    int index;// 0 - no, 1 - yes
    double weight; 
    void printme(){
      cout << "coulomb: " << index << " " << weight << endl;
    }
  };
  COULOMB mycoulomb;
  //----------------------------------
  struct OUTPUT{
    bool flag;
    string outfile;
    double sample_time;
    int separate; // 0 - no, 1 -yes
    void printme(){
      cout << "output: " << outfile << " " << sample_time << " " << separate << endl;
    }
  };
  OUTPUT myoutput;
  //----------------------------------
  struct REALTRAP{
    bool flag;
    int nfiles;
    std::vector<string> filenames;
    void printme(){
      cout << "realtrap: " << endl;
      cout << "\t" << flush;
      for(int i = 0 ; i < filenames.size(); i++)
	cout << filenames[i] << " " << flush;
      cout << endl;
    }
	
  };
  REALTRAP myrealtrap;
    //----------------------------------
    struct POTTRAP{
        bool flag;
        string potmap_filename;
        double Bz;
        
        void printme(){
            cout << "pottrap: " << endl;
            cout << "\t" << flush;
            cout << "Map : " << potmap_filename << ", Bz = " << Bz <<flush;
            cout << endl;
        }
        
    };
    POTTRAP mypottrap;
    //----------------------------------
  struct IDEALTRAP{
    bool flag;
    double r_electrode;
    double Ud2;
    double B;
    void printme(){
      cout << "idealtrap:" << endl;
      cout << "\t " << r_electrode << " " << Ud2 << " " << B << endl;
    }
  };
  IDEALTRAP myidealtrap;
  //----------------------------------
  struct NONIDEALTRAP{
    bool flag;
    string trapconfigfile;
    void printme(){
      cout << "non-idealtrap: " << trapconfigfile << endl;
    }

  };
  NONIDEALTRAP mynonitrap;
  //----------------------------------
  struct IMPORTDATA{
    bool flag;
    string prev_simu_file;
    void printme(){
      cout << "import data: " << prev_simu_file << endl;
    }

  };
  IMPORTDATA myimportdata;
  //----------------------------------
  struct CALCTRAP{
    bool flag;
    double param1;
    double param2;
    void printme(){
      cout << "calctrap: " << param1 << " " << param2 << endl;
    }
  };
  CALCTRAP mycalctrap;
  //----------------------------------
  struct NE{
    bool flag;
    double time;
    void printme(){
      cout << "No excitation: " << time << endl;
    }
  };
  NE myne;
  //----------------------------------
  struct DEW{
    bool flag;
    int size;
    double time;
    string EigenLett;
    string Element;
    double freq_bias;
    double amp;
    void printme(){
      cout << "Dipole excitation without buffer gas: " << endl;
      if(size == 6)cout << "\t " << time << " " << EigenLett << " " << Element << " " << freq_bias << " " << amp << endl;
      if(size == 5)cout << "\t " << time << " " << EigenLett << " " << freq_bias << " " << amp << endl;
    }
  };
  DEW mydew;
  //----------------------------------
  struct DEB{
    bool flag;
    int size;
    double time;
    string EigenLett;
    string Element;
    double freq_bias;
    double amp;
    void printme(){
      cout << "Dipole excitation with buffer gas: " << endl;
      if(size == 6)cout << "\t " << time << " " << EigenLett << " " << Element << " " << freq_bias << " " << amp << endl;
      if(size == 5)cout << "\t " << time << " " << EigenLett << " " << freq_bias << " " << amp << endl;
    }
  };
  DEB mydeb;
  //----------------------------------
  struct QEW{
    bool flag;
    int size;
    double time;
    string EigenLett;
    string Element;
    double freq_bias;
    double amp;
    void printme(){
      cout << "Quadrupole excitation without buffer gas: " << endl;
      if(size == 6)cout << "\t " << time << " " << EigenLett << " " << Element << " " << freq_bias << " " << amp << endl;
      if(size == 5)cout << "\t " << time << " " << EigenLett << " " << freq_bias << " " << amp << endl;
    }
  };
  DEW myqew;
  //----------------------------------
  struct QEB{
    bool flag;
    int size;
    double time;
    string EigenLett;
    string Element;
    double freq_bias;
    double amp;
    void printme(){
      cout << "Quadrupole excitation with buffer gas: " << endl;
      if(size == 6)cout << "\t " << time << " " << EigenLett << " " << Element << " " << freq_bias << " " << amp << endl;
      if(size == 5)cout << "\t " << time << " " << EigenLett << " " << freq_bias << " " << amp << endl;
    }
  };
  QEB myqeb;
  //----------------------------------
  struct OEW{
    bool flag;
    int size;
    double time;
    string EigenLett;
    string Element;
    double freq_bias;
    double amp;
    void printme(){
      cout << "Octupole excitation without buffer gas: " << endl;
      if(size == 6)cout << "\t " << time << " " << EigenLett << " " << Element << " " << freq_bias << " " << amp << endl;
      if(size == 5)cout << "\t " << time << " " << EigenLett << " " << freq_bias << " " << amp << endl;
    }
  };
  OEW myoew;
  //----------------------------------
  struct OEB{
    bool flag;
    int size;
    double time;
    string EigenLett;
    string Element;
    double freq_bias;
    double amp;
    void printme(){
      cout << "Octupole excitation with buffer gas: " << endl;
      if(size == 6)cout << "\t " << time << " " << EigenLett << " " << Element << " " << freq_bias << " " << amp << endl;
      if(size == 5)cout << "\t " << time << " " << EigenLett << " " << freq_bias << " " << amp << endl;
    }
  };
  OEB myoeb;
  //----------------------------------
  struct RWW{
    bool flag;
    double time;
    int order;
    double frequency;
    double amp;
    void printme(){
      cout << "Rotating wall excitation without buffer gas: " << endl;
      cout << "\t " << time << " " << order << " " << frequency << " " << amp << endl;
    }
  };
  RWW myrww;
  //----------------------------------
  struct RWB{
    bool flag;
    double time;
    int order;
    double frequency;
    double amp;
    void printme(){
      cout << "Rotating wall excitation with buffer gas: " << endl;
      cout << "\t " << time << " " << order << " " << frequency << " " << amp << endl;
    }
  };
  RWB myrwb;  
  //----------------------------------
  struct AWW{
    bool flag;
    double time;
    int order;
    double frequency;
    double amp;
    void printme(){
      cout << "Anti-Rotating wall excitation without buffer gas: " << endl;
      cout << "\t " << time << " " << order << " " << frequency << " " << amp << endl;
    }
  };
  AWW myaww;
  //----------------------------------
  struct AWB{
    bool flag;
    double time;
    int order;
    double frequency;
    double amp;
    void printme(){
      cout << "Anti-Rotating wall excitation with buffer gas: " << endl;
      cout << "\t " << time << " " << order << " " << frequency << " " << amp << endl;
    }
  };
  AWB myawb;  
  //----------------------------------
  struct AC{
    bool flag;
    double time;
    int order;
    string Element;
    double amp;
    void printme(){
      cout << "Axial coupling excitation: " << endl;
      cout << "\t " << time << " " << order << " " << Element << " " << amp << endl;
    }
  };
  AC myac;  
  //----------------------------------
  struct SC{
    bool flag;
    double time;
    int order;
    string Element;
    double freq1;
    double freq2;
    double amp1;
    double amp2;
    void printme(){
      cout << "Simco excitation: " << endl;
      cout << "\t " << time << " " << order << " " << freq1 << " " << freq2 << " " << amp1 << " " << amp2 << endl;
    }
  };
  SC mysc;  
    //----------------------------------
  struct AR{
    bool flag;
    double time;
    double amp1;
    double freq1;  
    double freq2;
    void printme(){
      cout << "AR excitation: " << endl;
      cout << "\t " << time << " " << amp1 << " " << freq1 << " " << freq2 << endl;
    }
  };
  AR myar; 
   //----------------------------------
  struct FB{
    bool flag;
    double time;
    double amp1;
    double freq1;  
    double freq2;
    void printme(){
      cout << "AR excitation: " << endl;
      cout << "\t " << time << " " << amp1 << " " << freq1 << " " << freq2 << endl;
    }
  };
  FB myfb; 
   //----------------------------------
    struct SWEEP{
        bool flag;
        double wi;
        double wf;
        double par;
        void printme(){
            cout << "SWEEP pulsation from " << wi << " to " << wf << endl;
        }
    };
    SWEEP mysweep;
    //----------------------------------

};

#endif
