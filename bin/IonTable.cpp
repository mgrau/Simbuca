//  IonTable.cpp
#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "IonTable.h"

using namespace std;

IonTable::IonTable(string filename_table_mass)
{   
    ifstream file;
    file.open(filename_table_mass.c_str(),ios::in);
    string dump;
    string line;
    int N_tmp;
    int A_tmp;
    int Z_tmp;
    int int_tmp;
    double mass_tmp;
    stringstream ss;
    if(!file)
    {   SLogger slogger("trapparameters");
        slogger << ERROR << "file " << filename_table_mass << " not found, mass table cannot be loaded!" << SLogger::endmsg;
	exit (EXIT_FAILURE);
    }
    else
    {
        while(getline(file,line))
        {
            if(line.size()>115 && line.size()!=124)
            {
                if(line.find("#") != string::npos)
                {
                    do{line.replace(line.find("#"),1," ");}
                    while(line.find("#") != string::npos);
                }
                
                unsigned short num;
                char c[14];
                
                num = line.copy(c,3,16);
                c[num] = '\0';
                A_tmp = atoi(c);
                
                num = line.copy(c,3,11);
                c[num] = '\0';
                Z_tmp = atoi(c);
                
                N_tmp = A_tmp - Z_tmp;
                
                num = line.copy(c,3,96);
                c[num] = '\0';
                int_tmp = atoi(c);
                
                num = line.copy(c,13,100);
                c[num] = '\0';
                mass_tmp = atof(c);
                
                
                dump = line.substr(20,2);
                if(dump.find(" ") != string::npos)
                {
                    dump.erase(dump.find(" "),1);
                }
                
                ionN.push_back(N_tmp);
                ionZ.push_back(Z_tmp);
                ionA.push_back(A_tmp);
                ionEl.push_back(dump);
                
                ss << A_tmp << dump;
                ionName.push_back(ss.str());
                ss.str("");
                
                mass_tmp = int_tmp + mass_tmp*1e-6;
                ionMass.push_back(mass_tmp);
                cout.precision(9);
                //    cout << dump << "  " << A_tmp << "  " << Z_tmp << "  " <<N_tmp << "  " << int_tmp << " " <<mass_tmp <<   endl;
            }
            
        }
        N_ions = ionN.size();
        
    }
    file.close();
    return;
}

IonTable::~IonTable()
{
    
}

void IonTable::ReadTable()
{
    cout.precision(9);
    for(int i=0;i<N_ions;i++)
    {
        cout << ionName[i] << "  " << ionMass[i] << endl;
    }
    return;
}

bool IonTable::TestIonName(string name_)
{
	bool res=false;
	for(int i=0;i<N_ions;i++)
    	{
        	if(ionName[i]==name_)
			res=true;
    	}
 return res;
}

bool IonTable::TestIonNames(vector<string > names_)
{
        bool res= true;
	for(unsigned j=0;j<names_.size();j++)
	{
		res &= TestIonName(names_[j]);
	}
	return !res;
}

double IonTable::ResearchMass(string name_)
{ 
    int i =-1;
    
    while(i<N_ions)
    {
        i++;
	
        if(name_==ionName[i])
        {
           return  ionMass[i];			  
        }
    }
  
    cout << "Ion not in the Table" << endl; 
    exit(EXIT_FAILURE);
   


}
