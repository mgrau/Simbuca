//
//  IonTable.h
//
//
//  Created by pierre dupre on 27/09/12.
//
//

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

class Operation {
public:
  Operation();
  Operation(string _name);
  ~Operation();
  void Launch(IonCloud &_cloud, _ode_vars & odev); 
  void Write();
  inline void SetTime(double t_){time_operation = t_;};
  inline void SetFrequency(double f_){frequency = f_;};
  inline void SetFrequency2(double f_){frequency2 = f_;};
  inline void SetFrequency3(double f_){f3 = f_;};
  inline void SetFrequency4(double f_){f4 = f_;};
  inline void SetFrequencyBias(double f_){f_bias = f_;};
  inline void SetAmplitude(double a_){amplitude = a_;};
  inline void SetAmplitude2(double a_){amplitude2 = a_;};
  inline void SetAmplitude3(double a_){a3 = a_;};
  inline void SetAmplitude4(double a_){a4 = a_;};
  inline void SetBuffBool(bool _buffergas){buffergas = _buffergas;};
  inline void SetBuff(double p_){p_buff_mbar = p_;};
  inline void SetOrder(double o_){order = o_;};
  inline void SetName(string name_){name=name_;};
  inline void SetIonName(string str_){ion_name=str_;};
  inline void SetFreqType(string str_){freqtype=str_;};
    inline double GetTime(){return time_operation;};
  inline double GetFrequency(){return frequency;};
  inline double GetAmplitude(){return amplitude;};
  inline double GetFrequency2(){return frequency2;};
  inline double GetAmplitude2(){return amplitude2;};
  //inline void SetIon(Ion i_){ion=i_;};
  inline string GetIonName(){return ion_name;};
  inline void SetFileName(string f_){file_name = f_;};
  
private:
  string name; // "NE" "DE" "QE" "OE" "RW" "DT"
  double time_operation;
  double p_buff_mbar;
  double frequency;
  double frequency2;
  double f_bias;
  double amplitude;
  double amplitude2;
  bool buffergas;
  int order;
  string ion_name;
  string freqtype;
    double f3;
    double f4;
    double a3;
    double a4;
    string file_name;
  //Ion ion;
};

class Operations {
public:
    Operations();
    ~Operations();
     void   Launch(IonCloud &_cloud, _ode_vars & odev);
     void AddOperation(Operation _Ope);
     inline int GetNumber_of_Ope(){return ope_N;};
     inline double GetTotal_time_of_Ope(){return ope_total_time;};
     void Write();
private:
    int ope_N; // Number of operations;
    double ope_total_time; // Time of all the operations
    vector<Operation > Ope;
};

#endif
