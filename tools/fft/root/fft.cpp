// Simon Van Gorp
// Created: 18 April 2013
// Modified: 18 April 2013
//
// 
// This program uses the root FFT function
//
//
// USAGE: command line inputs:
//		1) filename (the first line in the beginning will be skipped during reading)
//		2) the column with the time (in ms)
//		3) the column with the value of which you want the take the FFT
//
//	Code can be used command line or in root itself
//	command line: ./fft ../../CUSPsims/1pbarplasma_cloud_GaussEvo.txt 1 5
//      root: fft("../../CUSPsims/1pbarplasma_cloud_GaussEvo.txt",1,5)
//

#include <iostream>
using namespace std;
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <stdlib.h>     /* atoi */

//#include <iterator> 
//#include <math.h>
#include <sstream>

#include <TApplication.h>
#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "../../../bin/matrix.h"

bool verbose = true;

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
    int row=0;
    //clrscr();
    fptr.seekg(0,ios::beg);
    while(!fptr.eof()){
	 if( fptr.peek() != '#'){
         	fptr.getline(line,512);
     	 	row++;
		//cout<<row<<"\t"<<line<<endl;
	}else{
          getline(fptr,comment);
          comment.erase(comment.begin(), find_if(comment.begin(), comment.end(), not1(ptr_fun<int, int>(isspace))));
          //cout<<"comment found: "<<comment<<" \n";
	}
	
     }  
    row--; //since skip the endl on the end
    cout<<"nrows = "<<row<<endl;
    fptr.clear();
    return row;
}

void fft(char * filename, int col1, int col2, char * fftfilename){
    string comment;
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
    Matrix<double> m(nrow+1,2);

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
			m(row,0)=ix; //input in seconds
			m(row,1)=iy;
			n = 1;  
			row++;
			//cout<<row<<endl;
		}
	}
	else {
          getline(myfile,comment);
          comment.erase(comment.begin(), find_if(comment.begin(), comment.end(), not1(ptr_fun<int, int>(isspace))));
          //cout<<"comment found: "<<comment<<" \n";
	}
     }
	/* Define histogram with the information from the Matrix */	
	Double_t endtime = m(nrow-1,0); //in mseconds
	cout<<"endtime = "<<endtime<<" [ms]"<<endl;
	Int_t nbins = nrow;
	TH1D * hist = new TH1D("hist", "col1 vs col2", nbins, 0. , endtime);

	/* copy information from matrix to histogram */
	for(int i=0;i<nbins;i++){
		if(!verbose) cout<<i<<" "<<m(i,0)<<"\t"<<m(i,1)<<endl;
		hist->Fill(m(i,0),m(i,1));
		
	}


        //do the actual FFT:
 	TH1 *hm = 0;
	TVirtualFFT::SetTransform(0);
	hm = hist->FFT(hm, "MAG");
	hm->Draw();
	
	TCanvas *c = new TCanvas("c", "c", 1);
	c->Divide(2,1);
	c->cd(1);
	hist->Draw();
	ofstream fftfile;
	fftfile.open(fftfilename);
	//scale the xmax value of the FFt specytrum to get the correct frequency size;
	TH1D * hres = new TH1D("hres", "FFT", nbins, 0, nbins/endtime);
	for(int i= 1; i < hm->GetNbinsX(); i++){
		hres->SetBinContent(i, hm->GetBinContent(i));
		fftfile<<hm->GetXaxis()->GetBinCenter(i)*1000.<<"\t"<<hm->GetBinContent(i)<<endl;
	}
	hres->SetBinContent(1,0);
	fftfile.close();
	c->cd(2);
	hres->Draw();
	
	
}




# ifndef __CINT__
int main(int argc, char * argv[])
{
  if(argc != 4 && argc !=5){ 
	cout<<"The FFT requires 3 inputs:\n";
	cout<<"\t 1: filename\n \t 2: column #1(=time) is x1000 in the program \n \t 3: column #2(=amplitude)\n \t (4: v for verbose - optional)\n";
   }


  assert(argc == 5 || argc == 6);
  char * filename = argv[1];
  int col1 = atoi(argv[2]);
  int col2 = atoi(argv[3]);
  char * fftfilename = argv[4];

  if(argc == 6){
  	assert(argv[5][0] == 'v');
	verbose = false;
  }
  TApplication* theApp = new TApplication("App", &argc, argv);


  fft(filename,col1,col2,fftfilename);

  theApp->Run();
  
  return 0;
}
# endif

