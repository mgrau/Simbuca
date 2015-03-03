#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "ioncloud.h"
#include "ionFly.h"
#include "Operation.h"

using namespace std;

Operation::Operation()
{
  name="NE";
  time_operation=1e-9;
  p_buff_mbar = 0.;
  frequency = 0.;
  amplitude = 0.;
}

Operation::Operation(string _name)
{
  name=_name;
  time_operation=1e-9;
  p_buff_mbar = 0.;
  frequency = 0.;
  amplitude = 0.;
  buffergas = true;
}

Operation::~Operation(){};

void Operation::Launch(IonCloud &_cloud, _ode_vars & odev)
{
    stringstream ss_tmp;
    
  if(name=="NE")
  {
  	DoNoExcitation(time_operation,buffergas,p_buff_mbar/1000.0,_cloud,odev);
  }
  {
      DoSIMCOWithoutBuffergas(time_operation,frequency,amplitude,frequency2,amplitude2,_cloud,odev);
  }
    
  return;
}

void Operation::Write()
{
    cout.precision(9);

	if(name=="NE")
	{
		cout << "No excitation during " << time_operation << " s\n" << endl;
	}
    if(name=="DE"||name=="QE"||name=="OE")
    {
        if(name=="DE")
            cout << "Dipole excitation ";
        if(name=="QE")
            cout << "Quadrupole excitation ";
        if(name=="OE")
            cout << "Octupole excitation" ;
        
		if(p_buff_mbar==0.)
		{
			cout << "without buffer gas ";
        }
        else
        {
            cout << "with buffer gas ";
        }
        cout << "during " << time_operation << " s" << endl;
        cout << "\twith frequency: " << frequency << "Hz" << endl;
        cout << "\twith amplitude: "<<  amplitude << " V\n" << endl;
        
	}
	
	string str_order;
	
	if(name=="RW")
	{
		cout.precision(9);
		switch ( order )
      		{ 
			case 1:         
           		cout <<"Dipole RW ";
            		break;
         		case 2:
			cout<<"Quadrupole RW ";
            		break;
         		case 3:
			cout<<"Octupole RW ";
            		break;
			default:
			cout <<"Dipole RW ";order=1;
			break;
		}
		if(p_buff_mbar==0) 
			cout << "without buffer gas ";
		else
			cout << "with buffer gas ";	
		cout << "during " << time_operation << " s" << endl; 
		cout << "\twith frequency : "<< frequency << " Hz " << endl;
		cout << "\twith amplitude: "<<  amplitude << " V\n" << endl;
	}
	if(name=="AW")
	{
		cout.precision(9);
		switch ( order )
      		{ 
			case 1:         
           		cout <<"Dipole Anti-RW ";
            		break;
         		case 2:
			cout<<"Quadrupole Anti-RW ";
            		break;
         		case 3:
			cout<<"Octupole Anti-RW ";
            		break;
			default:
			cout <<"Dipole Anti-RW ";order=1;
			break;
		}
		if(p_buff_mbar==0) 
			cout << "without buffer gas ";
		else
			cout << "with buffer gas ";	
		cout << "during " << time_operation << " s" << endl; 
		cout << "\twith frequency : "<< frequency << " Hz " << endl;
		cout << "\twith amplitude: "<<  amplitude << " V\n" << endl;
	}
	if(name=="AC")
	{
		cout.precision(9);
		cout << "Axial Coupling excitation without buffer gas during " << time_operation << " s" << endl; 
		cout << "\twith frequency : " << frequency << "Hz" ;
		cout << "\twith amplitude: "<<  amplitude << " V\n" << endl;	
	}
    if(name=="SC")
    {
		cout << "SIMCO  excitation without buffer gas during " << time_operation << " s" << endl;
		cout << "\twith dipolar frequency : " << frequency << "Hz" ;
		cout << "\twith amplitude: "<<  amplitude << " V" << endl;
		cout << "\twith quadrupolar frequency : " << frequency2 << "Hz" ;
		cout << "\twith amplitude: "<<  amplitude2 << " V\n" << endl;
    }
    if(name=="AR")
    {
		cout << "Auto Resonance drive for " << time_operation << " s" << endl;
		cout << "\twith amplitude: "<<  amplitude << " V" << endl;
 		cout << "\twith start frequency : " << frequency << " Hz" <<endl ;
 		cout << "\twith stop frequency : " << frequency2 << " Hz"<<endl ;
    }
    if(name=="FB")
    {
		cout << "FB excitation drive for " << time_operation << " s" << endl;
		cout << "\twith amplitude: "<<  amplitude << " V" << endl;
 		cout << "\twith start frequency : " << frequency << " Hz" <<endl ;
 		cout << "\twith stop frequency : " << frequency2 << " Hz"<<endl ;
    }
    if(name=="SWIFT")
    {
        cout << "SWIFT dipole excitation  for " << time_operation << " s" << endl;
    }
    if(name=="RW_SWIFT")
    {
        cout << "Rotating Wall SWIFT dipole excitation  for " << time_operation << " s" << endl;
    }
    if(name=="EXC_EMAP")
    {
        cout << "Excitation with E map : " << file_name <<" during " << time_operation << " s" << endl;
    }
    if(name=="PI_PULSE")
    {
        cout << "Pi PULSE excitation during " << time_operation << " with amplitude : " << amplitude<< " V" << endl;
        
    }
	return;
}

Operations::Operations(){ope_N=0;ope_total_time=0;};

Operations::~Operations(){};

void Operations::AddOperation(Operation _Ope)
{
	Ope.push_back(_Ope);
    ope_total_time+= _Ope.GetTime();
	ope_N +=1;
	return;
}

void Operations::Launch(IonCloud &_cloud, _ode_vars & odev)
{
 if(ope_N!=0) 
 {
    SetTotalTime_of_Simu(ope_total_time,odev);
 	for(int i=0;i<ope_N;i++)
	{	
		Ope[i].Launch(_cloud, odev);
	}
 }
 else
 {
     cout << "No operation loaded" << endl;
     exit(-1);
 }
#ifdef __MPI_ON__
     odev.delete_mpiv();
#endif // _MPI_ON__
     
 
 return;
}

void Operations::Write()
{
 if(ope_N!=0) 
 {
 	for(int i=0;i<ope_N;i++)
	{	
		cout << "Operation " << i << " : ";
		Ope[i].Write();
	}
 }
 else
 {
 	cout << "No Operation loaded " << endl;
 }
 return;

}
