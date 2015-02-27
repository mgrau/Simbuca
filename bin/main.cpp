#include "MPI_simbuca.h"


#include <iostream>
//using namespace std;
#include <fstream>
#include "Initialization.h"
#include "SLogger.h"
#include <sstream>


#ifdef __NBODY_ON__
#include "nbody.h"
#endif //__NBODY_ON__

#ifdef __GUI_ON__
#include "mainwindow.h"
#include <QApplication>
#endif

int main(int argc,  char* argv[] ){
string simFile;

#ifdef __GUI_ON__
      QApplication a(argc, argv);
      MainWindow w;
      w.show();
      w.setWindowTitle("Simbuca GUI");
#else //__GUI_ON__

  if(argc==2){
      simFile = argv[1];
      // SIMBUCA FILE INPUT
      string extension_simbuca = ".sim";
      size_t found_ext;
      found_ext=simFile.rfind(extension_simbuca);
      if(found_ext!=string::npos){
          SimParser sparser(simFile.c_str());
          Run(argc,argv,sparser);
      }
  }
  else{ // argc==2
      SLogger slogger("trapparameters");
      slogger << ERROR << "Simbuca requires an inputfile\n" << SLogger::endmsg;exit(1);
  }  // argc!=2 .sim file

#endif //__GUI_ON__

#ifdef __GUI_ON__
    return a.exec();
#else
     return 0;
#endif //__GUI_ON__
}




