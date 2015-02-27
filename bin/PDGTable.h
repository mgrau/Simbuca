//
//  PDGTable.h
//
//
//  Created by balint radics on 08/02/13.
//
//

////////////////////////////////////////////////////////////////////////
//
//  Particle database manager class
//
//  This manager creates a list of particles which by default is
//  initialised from with the constants used by PYTHIA6 (plus some
//  other particles added). See definition and the format of the default
//  particle list in libraries/pdg/pdg_table.txt
//
////////////////////////////////////////////////////////////////////////


#ifndef _PDGTable_h
#define _PDGTable_h

#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "SLogger.h"

using namespace std;
class PDGTable {
  
 public:
  PDGTable();
  PDGTable(string filename_table_pdg);
  ~PDGTable();

  void ReadPDGTable(string filename_table_pdg);
  
  string GetPDGName(int pdgid);
  int GetPDGId(string name = "");
  double GetPDGMass(int pdgid);
  double GetPDGCharge(int pdgid);
  
  bool TestPDGName(string name_);
  bool TestPDGNames(vector<string > names_);
  
  void Print(); // print full table
  void Print(int pdgid); // print properties of particle with id
  void Print(vector<int> pdgids);  

  void MakeVerbose(bool verbose = false);
  
 private:
  vector<string> pName;
  vector<int> pId;
  vector<double> pMass;
  vector<double> pCharge;

  bool pVerbose;
  
  // to do: manage decay branching ratios and daughter particles

};

#endif
