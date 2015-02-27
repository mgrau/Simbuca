#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <time.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TRegexp.h>
#include <TROOT.h>

using namespace std;


std::vector<string> * list_files(const char *dirname="C:/root/folder/", TRegexp expression=".root")
{
  std::vector<string> * ret = new std::vector<string>;
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
	if (!file->IsDirectory() && fname.Contains(expression)) {
	//cout << fname.Data() << endl;
	ret->push_back(std::string(fname.Data()));
      }
    }
  }
  return ret;
}


void draw(){

  clock_t start = clock();
  /*Do something*/

  float labelsize = 0.05; 
  TRegexp expression("p1.txt");
  std::vector<string> * myfiles = list_files(".", expression);
  if(myfiles->size() == 0){
    cerr << "No files found - exiting..." << endl;
    exit(-1);
  }

  // grand map of coordinates
  map<int , map<int, double > > mt_x;
  map<int , map<int, double > > mt_y;
  map<int , map<int, double > > mt_z;
  map<int , map<int, double > > mt_vx;
  map<int , map<int, double > > mt_vy;
  map<int , map<int, double > > mt_vz;
  map<int , map<int, double > > mt_E;

  std::vector< double > v_imq;
  std::vector< double > v_imi;
  std::vector< double > v_EED;
  std::vector< double > v_C_Q;
  std::vector< double > v_I_Qs;
  std::vector< double > v_C;
  std::vector< double > v_I_A;

  std::vector<double> master_time;
  int i = 0;
  int Nmax = 1000;
  for(std::vector<string>::iterator it = myfiles->begin(); it!=myfiles->end(); ++it){
    i++;
    if(i == Nmax)break;

    string mfile = *it;
    std::cout << "opening particle file : " << mfile << std::endl;
    const char * filename = mfile.c_str();

    ifstream ifile;
    ifile.open(filename);

    // skip comments
    char comment= '#';
    int ncols = 14;
    int index = 0;
    double x,y,z, vx, vy, vz, rp, rm, R, E, Temp, t;
    char input[100];
    char name[100];
    
    int n = 0; 
    int c = 0;
    std::map< int, double > pmt_x, pmt_y, pmt_z, pmt_R, pmt_vx, pmt_vy, pmt_vz, pmt_E;
    while(!ifile.eof()){
      std::string str;
      std::getline(ifile, str);
      if(str[0] == '#' || str.size() == 0)continue;
      std::stringstream ss(str);
      ss >> index >> name >> x >> y >> z >> vx >> vy >> vz >> rp >> rm >> R >> E >> Temp >> t ;
            //std::cout << "line " << c << ": " <<  index <<  ", " << x << ", " << y << ", " << z << std::endl;cin.get();
      pmt_x.insert(std::pair<int, double>(c, x));
      pmt_y.insert(std::pair<int, double>(c, y));
      pmt_z.insert(std::pair<int, double>(c, z));
      pmt_vx.insert(std::pair<int, double>(c, vx));
      pmt_vy.insert(std::pair<int, double>(c, vy));
      pmt_vz.insert(std::pair<int, double>(c, vz));
      pmt_E.insert(std::pair<int, double>(c, E));
      c++;      
      if(i == 1){
	master_time.push_back(t);
      }
    }

    mt_x.insert(std::pair< int, std::map<int, double > >(index, pmt_x));
    mt_y.insert(std::pair< int, std::map<int, double > >(index, pmt_y));
    mt_z.insert(std::pair< int, std::map<int, double > >(index, pmt_z));
    mt_vx.insert(std::pair< int, std::map<int, double > >(index, pmt_vx));
    mt_vy.insert(std::pair< int, std::map<int, double > >(index, pmt_vy));
    mt_vz.insert(std::pair< int, std::map<int, double > >(index, pmt_vz));
    mt_E.insert(std::pair< int, std::map<int, double > >(index, pmt_E));

    ifile.close();
  }
// read image charge file
   expression = "imagecharges.txt";	  
   std::vector<string> * myfiles_images = list_files(".", "imagecharges.txt");
   i = 0;

   for(std::vector<string>::iterator it = myfiles_images->begin(); it!=myfiles_images->end(); ++it){
    i++;
    if(i == Nmax)break;

    string mfile = *it;
    std::cout << "opening imagecharges file : " << mfile << std::endl;
    const char * filename = mfile.c_str();

    ifstream ifile;
    ifile.open(filename);

    // skip comments
    char comment= '#';
    int ncols = 10;
    int index = 0;
    double t,z,vz,imq,imi,EED,C_Q,I_Qs,C,I_A;
    int n = 0; 
    int c = 0;

    while(!ifile.eof()){
      std::string str;
      std::getline(ifile, str);
      if(str[0] == '#' || str.size() == 0)continue;
      std::stringstream ss(str);
      //# (1)t (2)z(mm)   (3)vz    (4)im_q   (5)im_i (6)EED (7)BesselCharge[q]  (8)BesselCurrent[q/s] (9)BesselCharge[C] (10)BesselCurrent[A]
      ss >> t >> z >> vz >> imq >> imi >> EED >> C_Q >> I_Qs >> C >> I_A;
      //      std::cout << "line " << c << ": " <<  t <<  ", " << C_Q << ", " << I_Qs << ", " << C << ", " << I_A << std::endl;cin.get();
      v_imq.push_back(imq);
      v_imi.push_back(imi);
      v_EED.push_back(EED);
      v_C_Q.push_back(C_Q);
      v_I_Qs.push_back(I_Qs);
      v_C.push_back(C);
      v_I_A.push_back(I_A*1e15);
      c++;      
      if(i == 1){
	master_time.push_back(t);
      }
    }
    ifile.close();
  }
// end of read image charges file



  clock_t end1 = clock();
  float seconds = (float)(end1 - start) / CLOCKS_PER_SEC;
  std::cout << "reading in master data: " << seconds  << " sec" << std::endl;


  TCanvas * c = new TCanvas("c", "c", 10, 10, 1600, 1000 );
  c->Divide(2,2);
  TGraph * graphxy;
  TGraph * graphxy_u;
  TGraph * graphyz;
  TGraph * graphyz_u;
  TGraph * grapht_az;
  TGraph * graphzvz;
  TH1F * h_Er;
  TH1F * h_Ez;

  h_Ez = (TH1F*)gROOT->FindObject("h_Ez");
  if(h_Ez)delete h_Ez;
  h_Ez = new TH1F("h_Ez", "Axial Energy", 40.0, 0.0, 40.0);
  h_Ez->SetFillColor(kBlue);
  c->cd(3);
  h_Ez->Draw("hist");



  // h_Er = (TH1F*)gROOT->FindObject("h_Er");
  // if(h_Er)delete h_Er;
  // h_Er = new TH1F("h_Er", "Radial Energy", 100, 0.0, 100.);
  // h_Er->SetFillColor(kBlue);
  // c->cd(4);
  // h_Er->Draw("hist");

  std::vector< double > cxv, cyv, czv;
  std::vector< double > uxv, uyv, uzv;
  std::vector< double > tv, vaz, vvz;
  std::vector< double > ctI_A;


  graphxy = (TGraph*)(gROOT->FindObject("graphxy"));
  if(graphxy)delete graphxy;
  graphxy = new TGraph();
  // graphxy = new TGraph(cxv.size(), &cxv[0], &cyv[0]);
  graphxy->SetName("graphxy");
  graphxy->SetMarkerStyle(7);
  graphxy->SetMarkerColor(kBlue);
  
  graphxy_u = (TGraph*)(gROOT->FindObject("graphxy_u"));
  if(graphxy_u)delete graphxy_u;
  graphxy_u = new TGraph();
  // graphxy_u = new TGraph(uxv.size(), &uxv[0], &uyv[0]);
  graphxy_u->SetName("graphxy_u");
  graphxy_u->SetMarkerStyle(8);
  graphxy_u->SetMarkerColor(kRed);
  
  graphyz = (TGraph*)(gROOT->FindObject("graphyz"));
  if(graphyz)delete graphyz;
  graphyz = new TGraph();
  // graphyz = new TGraph(cyv.size(), &czv[0], &cyv[0]);
  graphyz->SetName("graphyz");
  graphyz->SetMarkerStyle(7);
  graphyz->SetMarkerColor(kBlue);
  
  graphyz_u = (TGraph*)(gROOT->FindObject("graphyz_u"));
  if(graphyz_u)delete graphyz_u;
  graphyz_u = new TGraph();
  //  graphyz_u = new TGraph(uyv.size(), &uzv[0], &uyv[0]);
  graphyz_u->SetName("graphyz_u");
  graphyz_u->SetMarkerStyle(8);
  graphyz_u->SetMarkerColor(kRed);
  
  grapht_az = (TGraph*)(gROOT->FindObject("grapht_az"));
  if(grapht_az)delete grapht_az;
  grapht_az = new TGraph(tv.size(), &tv[0], &vaz[0]);
  grapht_az->SetName("grapht_az");
  grapht_az->SetMarkerStyle(7);
  grapht_az->SetMarkerColor(kRed);
  
  graphzvz = (TGraph*)(gROOT->FindObject("graphzvz"));
  if(graphzvz)delete graphzvz;
  graphzvz = new TGraph();
  // graphyz = new TGraph(cyv.size(), &czv[0], &cyv[0]);
  graphzvz->SetName("graphzvz");
  graphzvz->SetMarkerStyle(7);
  graphzvz->SetMarkerColor(kBlue);

  map< int, map<int, double > >::iterator itx;
  map< int, map<int, double > >::iterator ity;
  map< int, map<int, double > >::iterator itz;
  map< int, map<int, double > >::iterator itvx;
  map< int, map<int, double > >::iterator itvy;
  map< int, map<int, double > >::iterator itvz;
  map< int, map<int, double > >::iterator itE;

  //std::map< int, double >::iterator iimq;
  //std::map< int, double >::iterator iimi;
  //std::map< int, double >::iterator iEED;
  //std::map< int, double >::iterator iC_Q;
  //std::map< int, double >::iterator iI_Qs;
  //std::map< int, double >::iterator iC;
  //std::map< int, double >::iterator iI_A;

  
  double cx, cy, cz, cvx, cvy, cvz, cEr, cEz, cE, az;
  int pid = 100;
  int count = 0;
  int pindex = 0;
  for(vector< double >::iterator itt = master_time.begin(); itt != master_time.end(); ++itt){
    count++;// step over master time
    double time = *itt;    
    std::cout << "current time: " << time << std::endl;
    tv.push_back(time);

    cxv.clear();
    cyv.clear();
    czv.clear();
    uxv.clear();
    uyv.clear();
    uzv.clear();
    //    vaz.clear();

    itx = mt_x.begin();
    ity = mt_y.begin();
    itz = mt_z.begin();
    itvx = mt_vx.begin();
    itvy = mt_vy.begin();
    itvz = mt_vz.begin();
    itE = mt_E.begin();

    // h_Er->Reset();
    h_Ez->Reset();
    clock_t start2 = clock();
    int cc =  0;// count files
    az = 0;
    map< int, double >::iterator mapx;
    map< int, double >::iterator mapy;
    map< int, double >::iterator mapz;
    map< int, double >::iterator mapvx;
    map< int, double >::iterator mapvy;
    map< int, double >::iterator mapvz;
    map< int, double >::iterator mapE;
    //loop over time
    for(; itx != mt_x.end() && ity != mt_y.end() && itz != mt_z.end() && itvx != mt_vx.end() && itvy != mt_vy.end() && itE != mt_E.end(); ++itx, ++ity, ++itz, ++itvx, ++itvy, ++itE){

      pindex = itx->first;
      mapx = itx->second.begin();
      mapy = ity->second.begin();
      mapz = itz->second.begin();
      mapvx = itvx->second.begin();
      mapvy = itvy->second.begin();
      mapvz = itvz->second.begin();
      mapE = itE->second.begin();
      int tc = 0;
      //loop over files
      for(; mapx!= itx->second.end() && mapy!=ity->second.end() && mapz != itz->second.end() &&
	    mapvx!= itvx->second.end() && mapvy!=itvy->second.end() && mapvz != itvz->second.end() &&
	    mapE!= itE->second.end(); ++mapx, ++mapy, ++mapz, ++mapvx, ++mapvy, ++mapvz, ++mapE){
	tc++; // step over proper time
	
	
	if(tc == count){
	  cx = mapx->second;
	  cy = mapy->second;
	  cz = mapz->second;
	  cvx = mapvx->second;
	  cvy = mapvy->second;
	  cvz = mapvz->second;
	  cE = mapE->second;
	  
	  cEr = cE*(cvx*cvx + cvy*cvy)/(cvx*cvx + cvy * cvy + cvz*cvz);
	  cEr = cEr*10000;
	  //      cEz = cE - cEr;
	  cEz = cE *(cvz*cvz)/(cvx*cvx + cvy * cvy + cvz*cvz);
	  //cout << cx << ", " << cy << ", " << cz << endl;
	  if (pindex != pid){
	    cxv.push_back(cx);
	    cyv.push_back(cy);
	    czv.push_back(cz);
	    //	cout << cE << endl;
	    //	    h_Er->Fill(fabs(cEr));
	    h_Ez->Fill(fabs(cEz));
	  }else if(pindex == pid){
	    uxv.push_back(cx);
	    uyv.push_back(cy);
	    uzv.push_back(cz);
	  }

	  az += cz;
	  break;
	}
      }
      cc++; // step over files
    }
    az = az/cc;
    vaz.push_back(az);
    vvz.push_back(cvz);
   
    clock_t end2 = clock();
    seconds = (float)(end2 - start2) / CLOCKS_PER_SEC;
    //std::cout << "finding data took: " << seconds << " sec " << endl;
    

    //    cout << "cxv size: " << cxv.size() << endl;
  
    clock_t start3 = clock();


 /*   c->cd(1);
    TH1F* h1 = c->cd(1)->DrawFrame(-1, -1, 1, 1);
    h1->SetTitle("Y vs. X");
    h1->GetXaxis()->SetTitle("X coordinate [mm]");
    h1->GetYaxis()->SetTitle("Y coordinate [mm]");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(labelsize);
    h1->GetYaxis()->SetLabelSize(labelsize);
    h1->GetXaxis()->SetTitleSize(labelsize);
    h1->GetYaxis()->SetTitleSize(labelsize);
    graphxy->DrawGraph(cxv.size(), &cxv[0], &cyv[0], "P");
    graphxy_u->DrawGraph(uxv.size(), &uxv[0], &uyv[0], "Psame");
    //graphxy->Draw("P");
    //  graphxy_u->Draw("Psame");
    c->cd(1)->SetGrid();
    c->cd(1)->Update();
*/
/*    c->cd(2);
    TH1F* h2 = c->cd(2)->DrawFrame(-240, -1, -60, 1);
    h2->GetXaxis()->SetTitle("Z coordinate [mm]");
    h2->GetYaxis()->SetTitle("Y coordinate [mm]");
    h2->SetTitle("Y vs. Z");
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();
    h2->GetXaxis()->SetLabelSize(labelsize);
    h2->GetYaxis()->SetLabelSize(labelsize);
    h2->GetXaxis()->SetTitleSize(labelsize);
    h2->GetYaxis()->SetTitleSize(labelsize);
    graphyz->DrawGraph(cyv.size(), &czv[0], &cyv[0], "P");
    graphyz_u->DrawGraph(uyv.size(), &uzv[0], &uyv[0], "Psame");
    c->cd(2)->SetGrid();
    c->cd(2)->Update();
*/
    c->cd(1);
    TH1F* h1 = c->cd(1)->DrawFrame(0., -150.,1., 150.);
    h1->SetTitle("Image current vs. time");
    h1->GetXaxis()->SetTitle("time [ms]");
    h1->GetYaxis()->SetTitle("Induced current on U5 [fA]");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->SetNdivisions(5, 0, 0);
    h1->GetXaxis()->SetLabelSize(labelsize);
    h1->GetYaxis()->SetLabelSize(labelsize);
    h1->GetXaxis()->SetTitleSize(labelsize);
    h1->GetYaxis()->SetTitleSize(labelsize);
    graphxy->DrawGraph(tv.size(), &tv[0], &v_I_A[0], "P");
    c->cd(1)->SetGrid();
    c->cd(1)->Update();

    c->cd(2);
    TH1F* h2 = c->cd(2)->DrawFrame(-190., -100000., -100., 100000.);
    h2->GetXaxis()->SetTitle("Z coordinate [mm]");
    h2->GetYaxis()->SetTitle("axial velocity [m/s]");
    h2->SetTitle("Z vs. V_z");
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();
    h2->GetXaxis()->SetLabelSize(labelsize);
    h2->GetYaxis()->SetLabelSize(labelsize);
    h2->GetXaxis()->SetTitleSize(labelsize);
    h2->GetYaxis()->SetTitleSize(labelsize);
    graphzvz->DrawGraph(tv.size(), &vaz[0], &vvz[0], "P");
    //graphyz_u->DrawGraph(uyv.size(), &uzv[0], &uyv[0], "Psame");
    c->cd(2)->SetGrid();
    c->cd(2)->Update();


    c->cd(3);
    h_Ez->Draw("hist");
    h_Ez->GetXaxis()->SetTitle("Energy [eV]");
    h_Ez->GetYaxis()->SetTitle("Entries");    
    h_Ez->GetXaxis()->CenterTitle();
    h_Ez->GetYaxis()->CenterTitle();
    h_Ez->GetXaxis()->SetLabelSize(labelsize);
    h_Ez->GetYaxis()->SetLabelSize(labelsize);
    h_Ez->GetXaxis()->SetTitleSize(labelsize);
    h_Ez->GetYaxis()->SetTitleSize(labelsize);
    // h_Er->Draw("hist");
    // h_Er->GetXaxis()->SetTitle("Energy [eV*1e+04]");
    // h_Er->GetYaxis()->SetTitle("Entries");    
    // h_Er->GetXaxis()->CenterTitle();
    // h_Er->GetYaxis()->CenterTitle();
    // h_Er->GetXaxis()->SetLabelSize(labelsize);
    // h_Er->GetYaxis()->SetLabelSize(labelsize);
    // h_Er->GetXaxis()->SetTitleSize(labelsize);
    // h_Er->GetYaxis()->SetTitleSize(labelsize);
    c->cd(3)->SetGrid();
    c->cd(3)->Update();

    c->cd(4);
    TH1F* h4 = c->cd(4)->DrawFrame(0, -190, 1., -100);
    h4->GetXaxis()->SetTitle("time [ms]");
    h4->GetYaxis()->SetTitle("Average Z coordinate [mm]");
    h4->SetTitle("Average Z vs. time");
    h4->GetXaxis()->CenterTitle();
    h4->GetYaxis()->CenterTitle();
    h4->GetXaxis()->SetNdivisions(5, 0, 0);
    h4->GetXaxis()->SetLabelSize(labelsize);
    h4->GetYaxis()->SetLabelSize(labelsize);
    h4->GetXaxis()->SetTitleSize(labelsize);
    h4->GetYaxis()->SetTitleSize(labelsize);
    grapht_az->DrawGraph(tv.size(), &tv[0], &vaz[0], "P");
    c->cd(4)->SetGrid();
    c->cd(4)->Update();

    clock_t end3 = clock();
    seconds = (float)(end3 - start3) / CLOCKS_PER_SEC;    
    //cout << "creating plots took: " << seconds << " sec " << endl;

    clock_t start4 = clock();

    char str[15];
    sprintf(str, "%d", count);
    string files = "test_long" + string(str) + ".png";
    c->Print(files.c_str());
    cout<<"made file: "<<files.c_str()<<endl;
    //    c->Print("test_long3.gif+");

    clock_t end4 = clock();
    seconds = (float)(end4 - start4) / CLOCKS_PER_SEC;    
    //cout << "saving plot took: " << seconds << " sec " << endl;

  }
  

}

# ifndef __CINT__
int main(int argc, char * argv[])
{

  gROOT->SetBatch( kTRUE );

  TApplication* theApp = new TApplication("App", &argc, argv);

  draw();

  theApp->Run();
  return 0;
}
# endif
