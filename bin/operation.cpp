#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "ioncloud.h"
#include "ionFly.h"
#include "operation.h"

using namespace std;

operation::operation()
{
    name="NE";
    time=1e-9;
}

operation::operation(string _name)
{
    name=_name;
    time=1e-9;
}

operation::~operation(){};

void operation::launch(IonCloud &_cloud, _ode_vars & odev)
{
    stringstream ss_tmp;
    if(name=="NE")
        DoNoExcitation(time,false,0.0,_cloud,odev);
}

void operation::write()
{
    cout.precision(9);
    if(name=="NE")
    {
        cout << "No excitation during " << time << " s\n" << endl;
    }
    if(name=="DE"||name=="QE"||name=="OE")
    {
        if(name=="DE")
            cout << "Dipole excitation ";
        if(name=="QE")
            cout << "Quadrupole excitation ";
        if(name=="OE")
            cout << "Octupole excitation" ;
        cout << "during " << time << " s" << endl;
    }
    string str_order;
    if(name=="RW")
    {
        cout.precision(9);
        cout << "during " << time << " s" << endl; 
    }
    if(name=="AW")
    {
        cout.precision(9);
        cout << "during " << time << " s" << endl; 
    }
    if(name=="AC")
    {
        cout.precision(9);
        cout << "Axial Coupling excitation without buffer gas during " << time << " s" << endl; 
    }
    if(name=="SC")
    {
        cout << "SIMCO  excitation without buffer gas during " << time << " s" << endl;
    }
    if(name=="AR")
    {
        cout << "Auto Resonance drive for " << time << " s" << endl;
    }
    if(name=="FB")
    {
        cout << "FB excitation drive for " << time << " s" << endl;
    }
    if(name=="SWIFT")
    {
        cout << "SWIFT dipole excitation  for " << time << " s" << endl;
    }
    if(name=="RW_SWIFT")
    {
        cout << "Rotating Wall SWIFT dipole excitation  for " << time << " s" << endl;
    }
}

operations::operations() {};

operations::~operations(){};

void operations::add(operation op)
{
    ops.push_back(op);
}

double operations::total_time() {
    double time = 0.0;
    for (std::vector<operation>::iterator op = ops.begin(); op != ops.end(); ++op)
        time += op->time;
    return time;
}

void operations::launch(IonCloud &_cloud, _ode_vars & odev)
{
    if(ops.size()!=0) 
    {
        SetTotalTime_of_Simu(total_time(),odev);
        for (std::vector<operation>::iterator op = ops.begin(); op != ops.end(); ++op)
            op->launch(_cloud, odev);
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

void operations::write()
{
    if(ops.size()!=0) 
    {
        int i = 0;
        for (std::vector<operation>::iterator op = ops.begin(); op != ops.end(); ++op)
        {
            cout << "Operation " << i++ << " : ";
            op->write();
        }
    }
    else
        cout << "No Operation loaded " << endl;
}
