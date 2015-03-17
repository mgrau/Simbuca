#include "trapparameters.h"

// int countword2(string line)
// {
//     int wordcount;
//     if (line.at(0)==' ')
//         wordcount=0;
//     else
//         wordcount=1;
//     for(int i=1; i<line.length() ;i++)
//         if (line[i] != ' ' && line[i-1] == ' ')
//             wordcount++;
//     //cout<<"The string has "<<wordcount<<" words"<<endl;
//     return wordcount ;
// }
_trap_param::_trap_param() {
    Vrf = 50;
    Vdc = 5;
    kappa = 0.001;
    r0 = 0.04;
    z0 = 0.04;
    freq_rf = 50e5;
    trap_type = 0;
}

//______________________________________________________________________________
_trap_param::~_trap_param() {
}

int _trap_param::ReadFile(string filename) {
    // ifstream file;
    // file.open(filename.c_str(),ios::in);
    // string line,comment;
    // char line2[256];
    // stringstream streamline;
    // stringstream streamline2;
    // int _countword;
    // B2 = 0.;
    // c4 =0.;
    // c6 =0.;
    // if(!file)
    // {   SLogger slogger("trapparameters");
    //     slogger << ERROR << "can`t find config file <" << filename << ">" << SLogger::endmsg;
    //     return 1;
    // }
    // else
    // {
    //     //check for comments in the beginning of the file
    //     while(file.peek() == '#'){
    //         getline(file,comment);
    //         comment.erase(comment.begin(), find_if(comment.begin(), comment.end(), not1(ptr_fun<int, int>(isspace))));
    //         //cout<<"comment found: "<<comment<<" \n";
    //     }
    //
    //     // 1st line
    //     getline(file,line); // MAGNETIC FIELD
    //     // 2nd line
    //     getline(file,line); // B0 (B2)
    //     _countword = countword2(line);
    //     switch(_countword)
    //     {
    //         case 1:
    //             streamline.str(line);
    //             streamline >> B0;
    //             break;
    //         case 2:
    //             streamline.str(line);
    //             streamline >> B0;
    //             streamline >> B2;
    //             break;
    //         default:
    //             SLogger slogger("trapparameters");
    //             slogger << ERROR << "wrong magnetic field parameters" << SLogger::endmsg;
    //             return 1;
    //     }
    //     getline(file,line); //ELECTRIC POTENTIAL
    //
    //     getline(file,line); //c2 (c4) (c6)
    //     _countword = countword2(line);
    //
    //     switch(_countword)
    //     {	
    //         case 1:
    //             streamline2.str(line.c_str());
    //             streamline2 >> Ud2;
    //
    //             break;
    //         case 2:
    //             streamline2.str(line);
    //             streamline2 >> U0;
    //             streamline2 >> c2;
    //             break;
    //         case 3:
    //             streamline2.str(line);
    //             streamline2 >> U0;
    //             streamline2 >> c2;
    //             streamline2 >> c4;
    //         case 4:
    //             streamline2.str(line);
    //             streamline2 >> U0;
    //             streamline2 >> c2;
    //             streamline2 >> c4;
    //             streamline2 >> c6;
    //             break;
    //         default:
    //             SLogger slogger("trapparameters");
    //             slogger << ERROR << "Error in electric potential  parameters" << SLogger::endmsg;
    //             return 1;
    //     }
    //     getline(file,line);//ELECTRODE RADIUS
    //     file >> r_0; //r0
    //
    //     getline(file,line);//CORRECTION FACTOR
    //     getline(file,line);
    //     file >> a;//a
    //
    //     getline(file,line);//PUMPING DIAPHRAGM
    //     getline(file,line);
    //     file >> PUMPING_DIAPRHAGM_RADIUS;//Rd
    //
    //     if (_countword != 1){
    //         Ud2 = U0*c2;
    //     }        
    //
    //
    //     trap_config = 1;
    //     file.close();
    // }
    //
    // return 0;
}

void _trap_param::Print() {
}
