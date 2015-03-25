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
    name="normal";
    time=1e-9;
    print_time = 0.0;
}

operation::operation(string _name)
{
    name=_name;
    time=1e-9;
    print_time= 0.0;
}

operation::~operation(){};

void operation::launch(IonCloud &cloud, _ode_vars & odev)
{
    if (print_time) {
        SetPrintInterval(print_time);
        if (odev.h > print_time) {
            odev.hnext = print_time;
            odev.h = print_time;
        }
    }
    odev.forcev.reset_ops();
    if (name=="normal") {
        cout << "normal operation launching" << endl;
        odev.forcev.trap = true;
        normal_operation(time,cloud,odev);
    }
    else if (name == "tof") {
        cout << "time of flight operation launching" << endl;
        odev.forcev.trap = false;
        odev.forcev.tof = true;
        normal_operation(time,cloud,odev);
    }
    else {
        cout << "free flight operation launching" << endl;
        normal_operation(time,cloud,odev);
    }
}

void operation::write()
{
    cout.precision(9);
    if (name == "normal")
        cout << "normal operation for " << time << " s" << endl;
    else if (name == "tof")
        cout << "time of flight kick for " << time << " s" << endl;
    else
        cout << "No operation defined for " << time << " s!" << endl;
}

operations::operations() {};

operations::~operations(){};

void operations::add(operation op) {
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
