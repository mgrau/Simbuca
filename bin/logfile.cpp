#include "logfile.h"

//definition of member functions
LogFile::LogFile(){
    //Do Nothing.
} 
 
LogFile::~LogFile(){
       logfile.close();
} 

void LogFile::open(char* filename){ 
       logfile.open(filename);
       logfile<<"logfile created on: \t";
       time_t rawtime;struct tm * timeinfo;time ( &rawtime );timeinfo = localtime ( &rawtime );
       logfile<<asctime (timeinfo);
}
 

void LogFile::operator <<(char* is){
      // XXX throw error if file is not opened
      logfile<<is;
      logfile<<flush;
}

void LogFile::operator <<(double d){
     logfile<<d;  
     logfile<<flush;   
}

void LogFile::operator <<(string s){
     logfile<<s;  
     logfile<<flush;   
}

/*
void LogFile::operator <<(int i){
     logfile<<i;     
}*/

//define static variable
ofstream LogFile::logfile; 


