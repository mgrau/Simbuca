// Simon Van Gorp
// Created: 18 April 2013
// Modified: 7 May 2013
//
//
// USAGE: command line inputs:
//		1) filename (lines starting with '#' will be skipped)
//		2) column 1
//		3) column 2

#include <iostream>
using namespace std;
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <stdlib.h>     /* atoi */
#include <sstream>
#include <algorithm>    // std::find_if
#include "../../bin/matrix.h"


int NCol(char *filename){
    //returns the number of columns in a file
    string comment;
    fstream fptr(filename,ios::in|ios::out);
    char ch = ' ';
    char line[512];
    int col=1;
    //clrscr();
    fptr.seekg(0,ios::beg);
    while(ch!='\n'){
	if( fptr.peek() != '#' && ch != '#'){
            fptr.get(ch);
            if(ch == '\t' || ch == ' '){
		col++;
	    //read all the following whitespaces and remove them
           while(fptr.peek() == ' ')fptr.get(ch);
	   }
	}else{     
          getline(fptr,comment);
          comment.erase(comment.begin(), find_if(comment.begin(), comment.end(), not1(ptr_fun<int, int>(isspace))));
          //cout<<"comment found: "<<comment<<" \n";
}
   }
    cout<<"ncols = "<<col<<endl;
    fptr.clear();
    return col;
}

int NRow(char *filename){
    //returns the number of columns in a file
    string comment;
    fstream fptr(filename,ios::in|ios::out);
    char line[512];
    int row=1;
    //clrscr();
    fptr.seekg(0,ios::beg);
    while(!fptr.eof()){
	 if( fptr.peek() != '#'){
         	fptr.getline(line,512);
     	 	row++;
	}else{
          getline(fptr,comment);
          comment.erase(comment.begin(), find_if(comment.begin(), comment.end(), not1(ptr_fun<int, int>(isspace))));
          //cout<<"comment found: "<<comment<<" \n";
	}
	
     }
    cout<<"nrows = "<<row<<endl;
    fptr.clear();
    return row;
}

void columnselect(char * filename, int col1, int col2){ 
    string comment;
    //read the columns in a matrix and print out the matrix!
    ifstream myfile;
    myfile.open(filename);
    //check if the filename exists then open
    if(!myfile.good())cout<<"Can't open file: "<<filename<<endl;
    assert(myfile.good());

    
    //check if the requested column is lower than the #columns of the file.
    int ncol = NCol(filename);
    if(col1 > ncol || col2 > ncol) cout<<"requested column is higher than "<<ncol<<". Which is the number of columns of the file"<<endl;
  assert(col1 <= ncol && col2 <=ncol );

    /*declare Matrix as a buffer to put the read values in */
    int nrow = NRow(filename);
    Matrix<double> m(nrow,2);

    /* Start reading the file */

    double ix,iy;
    double n = 1;
    char line[512];
    char input[100];
    //read the file and when at the correct column put the number in a file.
    int row=0;
    while(!myfile.eof()){
         myfile>>input;
	 if( input[0] != '#'){
		//cout<<n<<"\n"<<input<<endl;cin.get();
        	if(n == col1) ix = atof(input);
        	if(n == col2) iy = atof(input);
        	n++;
        	if( n > ncol) {
			//insert values in the matrix
			m(row,0)=ix;
			m(row,1)=iy;
			n = 1;  
			row++;
		}

	}
	else {
          getline(myfile,comment);
          comment.erase(comment.begin(), find_if(comment.begin(), comment.end(), not1(ptr_fun<int, int>(isspace))));
          //cout<<"comment found: "<<comment<<" \n";
	}

    }


	/* cout the information of the matrix  */
	for(int i=0;i<row;i++){
		//cout<<i<<" "<<m(i,0)<<"\t"<<m(i,1)<<endl;
		cout<<m(i,0)<<"\t"<<m(i,1)<<endl;
		
	}
}


int main(int argc, char * argv[])
{
  if(argc != 4){ 
	cout<<"The column selector requires 3 inputs:\n";
	cout<<"\t 1: filename\n \t 2: column #1 \n \t 3: column #2\n";
   }
  assert(argc == 4);
  char * filename = argv[1];
  int col1 = atoi(argv[2]);
  int col2 = atoi(argv[3]);

  columnselect(filename,col1,col2);
  
  return 0;
}


