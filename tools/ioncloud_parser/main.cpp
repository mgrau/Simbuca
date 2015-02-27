#include <iostream>
using namespace std;
#include <vector>
#include <fstream>
#include <math.h>
#include "../../bin/globals.h"
#include "../../bin/parser.h"
#include <sstream>
#include <assert.h>


int main (int argc,char *argv[]){	
   // argc is the number of items on the command line.
   // The first (argv[0]) is the command name (simple)
   // and then there is one real argument.  So argc 
   // should be 2.
   
   if(argc != 3) cout<<"first input is the filename-prefix, second input the number of files. \n";
   assert(argc == 3);

   // the first command-line input (argv[1]) is a ASCII string
   //   containing the digits of an integer.  atoi translates
   //   this string into an integer.
   if(argc == 3) {
   	int nrparticles = atoi(argv[2]);
   	char * filename_prefix = argv[1];

    	/* PostProces the data */
       	//PrintIonCloudinfo(filename_prefix, nrparticles, 2.0);
    	PrintIonCloudGaussEvo(filename_prefix, nrparticles, 1.0, true);
	//BundleData(filename_prefix,nrparticles);
    }
    
    //email_it("G100_bundled.txt" , "blabla@blablab.com");
    //see http://cc.byexamples.com/20070120/print-color-string-without-ncurses/
    //printf("%c[%d;%dmProgam Ended%c[%dm\n",27,4,31,27,0);    
	return 0;
}





