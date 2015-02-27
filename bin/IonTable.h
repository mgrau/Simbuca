//
//  IonTable.h
//
//
//  Created by pierre dupré on 27/09/12.
//
//

#ifndef _IonTable_h
#define _IonTable_h

#include "SLogger.h"
#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;
class IonTable {
public:
    IonTable();
    IonTable(string filename_table_mass);
    ~IonTable();
    inline int GetionN(int i_){return ionN[i_];};
    inline int GetionZ(int i_){return ionZ[i_];};
    inline int GetionA(int i_){return ionA[i_];};
    inline string GetionEl(int i_){return ionEl[i_];};
    inline string GetionName(int i_){return ionName[i_];};
    inline double GetionMass(int i_){return ionMass[i_];};
    inline int GetNtot(){return N_ions;};
    void ReadTable();
    bool TestIonName(string name_);
    bool TestIonNames(vector<string > names_);
    double ResearchMass(string name_);
private:
    int N_ions;// number of ion in the table
    vector<int > ionN; // atomic N
    vector<int > ionZ; // atomic Z
    vector<int > ionA; // atomic A
    vector<string > ionEl; // name of the element ex: Li for Lithium
    vector<string > ionName; // full name of the element ex 7Li
    vector<double> ionMass; // mass (amu)
    vector<double> ionMass_error; // mass error (micro amu)
};

#endif
