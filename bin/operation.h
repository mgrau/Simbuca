#ifndef _OPERATION_h
#define _OPERATION_h

#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "ioncloud.h"
#include "ion.h"

using namespace std;

struct operation {
    operation();
    operation(string _name);
    ~operation();
    void launch(IonCloud &_cloud, _ode_vars & odev); 
    void write();
    string name;
    double time;
    double print_time;
};

class operations {
    public:
        operations();
        ~operations();
        void launch(IonCloud &_cloud, _ode_vars & odev);
        void add(operation op);
        double total_time();
        void write();
    private:
        vector<operation> ops;
};

#endif
