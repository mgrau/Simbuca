//try to implement Matrix class, based on 
// http://www.parashift.com/c++-faq-lite/freestore-mgmt.html#faq-16.16

#ifndef LOGFILE_H
#define LOGFILE_H
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;



 class LogFile {
 public:
   LogFile();
  ~LogFile();
 
   //Access methods:
   void open(char* filename);
   void operator <<(char* is);
   void operator <<(double d);
   void operator <<(string s);
  // void operator <<(int i);
 private:
   static ofstream logfile; 
 };
 
#endif


