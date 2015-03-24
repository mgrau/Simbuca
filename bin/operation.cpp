#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "ioncloud.h"
#include "ion_fly.h"
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

void operation::launch(IonCloud &cloud, _ode_vars & odev)
{
    if(name=="normal")
        normal_operation(time,cloud,odev);
    else
        normal_operation(time,cloud,odev);
}

void operation::write()
{
    cout.precision(9);
    if (name == "normal")
        cout << "normal operation for " << time << " s" << endl;
    else
        cout << "No operation defined for " << time << " s!" << endl;
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
        odev.total_time = total_time();
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
