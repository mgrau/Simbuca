
#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdlib.h>

#include "SimParser.h"


SimParser::SimParser(const char * simfile){

  // map configuration labels
  label_map.insert(pair<string, const int>("CREATECLOUD", LCREATECLOUD));
  label_map.insert(pair<string, const int>("CLOUDPARTS", LCLOUDPARTS));
  label_map.insert(pair<string, const int>("CLOUDCOORD", LCLOUDCOORD));
  label_map.insert(pair<string, const int>("TEMP", LTEMP));
  label_map.insert(pair<string, const int>("CREATEPARTICLES", LCREATEPARTICLES));
  label_map.insert(pair<string, const int>("PARTICLES", LPARTICLES));
  label_map.insert(pair<string, const int>("BUFFER", LBUFFER));
  label_map.insert(pair<string, const int>("ODE", LODE));
  label_map.insert(pair<string, const int>("COULOMB", LCOULOMB));
  label_map.insert(pair<string, const int>("OUTPUT", LOUTPUT));
  label_map.insert(pair<string, const int>("IMPORTDATA", LIMPORTDATA));
  label_map.insert(pair<string, const int>("REALTRAP", LREALTRAP));
  label_map.insert(pair<string, const int>("POTTRAP", LPOTTRAP));
  label_map.insert(pair<string, const int>("IDEALTRAP", LIDEALTRAP));
  label_map.insert(pair<string, const int>("NONIDEALTRAP", LNONIDEALTRAP));
  label_map.insert(pair<string, const int>("CALCTRAP", LCALCTRAP));
  label_map.insert(pair<string, const int>("NE", LNE));
  label_map.insert(pair<string, const int>("DEW", LDEW));
  label_map.insert(pair<string, const int>("DEB", LDEB));
  label_map.insert(pair<string, const int>("QEW", LQEW));
  label_map.insert(pair<string, const int>("QEB", LQEB));
  label_map.insert(pair<string, const int>("OEW", LOEW));
  label_map.insert(pair<string, const int>("OEB", LOEB));
  label_map.insert(pair<string, const int>("RWW", LRWW));
  label_map.insert(pair<string, const int>("RWB", LRWB));
  label_map.insert(pair<string, const int>("AWW", LAWW));
  label_map.insert(pair<string, const int>("AWB", LAWB));
  label_map.insert(pair<string, const int>("AC", LAC));
  label_map.insert(pair<string, const int>("SC", LSC));
  label_map.insert(pair<string, const int>("AR", LAR));
  label_map.insert(pair<string, const int>("FB", LFB));
  label_map.insert(pair<string, const int>("SWEEP", LSWEEP));
  label_map.insert(pair<string, const int>("SWIFT", LSWIFT));
  label_map.insert(pair<string, const int>("EXC_EMAP", LEXC_EMAP));
  label_map.insert(pair<string, const int>("PI_PULSE", LPI_PULSE));
  label_map.insert(pair<string, const int>("ASYM_ARW", LASYM_ARW));
    label_operation_map.insert(pair<string, const int>("NE", LNE));
    label_operation_map.insert(pair<string, const int>("DEW", LDEW));
    label_operation_map.insert(pair<string, const int>("DEB", LDEB));
    label_operation_map.insert(pair<string, const int>("QEW", LQEW));
    label_operation_map.insert(pair<string, const int>("QEB", LQEB));
    label_operation_map.insert(pair<string, const int>("OEW", LOEW));
    label_operation_map.insert(pair<string, const int>("OEB", LOEB));
    label_operation_map.insert(pair<string, const int>("RWW", LRWW));
    label_operation_map.insert(pair<string, const int>("RWB", LRWB));
    label_operation_map.insert(pair<string, const int>("AWW", LAWW));
    label_operation_map.insert(pair<string, const int>("AWB", LAWB));
    label_operation_map.insert(pair<string, const int>("AC", LAC));
    label_operation_map.insert(pair<string, const int>("SC", LSC));
    label_operation_map.insert(pair<string, const int>("AR", LAR));
    label_operation_map.insert(pair<string, const int>("FB", LFB));
    label_operation_map.insert(pair<string, const int>("SWIFT", LSWIFT));
    label_operation_map.insert(pair<string, const int>("EXC_EMAP", LEXC_EMAP));
    label_operation_map.insert(pair<string, const int>("PI_PULSE", LPI_PULSE));
    label_operation_map.insert(pair<string, const int>("ASYM_ARW", LASYM_ARW));
  // initialize flags to false
  myimportdata.flag = mycreatecloud.flag = mycloudparts.flag = mycloudcoord.flag = mytemp.flag = mycreateparticles.flag = myparticles.flag = mybuffer.flag = myode.flag = mycoulomb.flag = myoutput.flag = myidealtrap.flag = myrealtrap.flag = myne.flag = mydew.flag = mydeb.flag = myqew.flag = myqeb.flag = myoew.flag = myoeb.flag = myrww.flag = myrwb.flag = myaww.flag = myawb.flag = myac.flag = mysc.flag = mypottrap.flag = mynonitrap.flag =mycalctrap.flag = myar.flag=myfb.flag=mysweep.flag=false;

  infile.open(simfile);
  if(infile){

    string line;
    while(!infile.eof()){
       std::getline(infile, line);
      if(line[0] == '#')comments.push_back(line);
      string word;
      std::vector<string> fields;
      stringstream mystream(line);
      while(getline(mystream, word, ' ')){
	fields.push_back(word);
      }
      if(fields.size() > 0 && fields[0] != "#" )
      {
          if( label_operation_map[fields[0]]!=0)
          {
              //cout << fields[0] << " is an operation" << endl;
              operation_map.push_back(fields);
          }
          else
          {
             config_map.push_back(std::pair<std::string, std::vector<string> >(fields[0], fields));
          }
      }
    }

    infile.close();
    //  PrintFile();
  
    ProcessFile();
      
  }else{
    SLogger slogger("simparser");
    slogger << ERROR << "Simbuca inputfile <" << simfile << "> does not exist!" << SLogger::endmsg;
    exit(-1);
  }

}

SimParser::~SimParser(){}

void SimParser::ProcessFile(){


  // Loop through all the label strings and see if it matches any
  for(vector< pair < string, std::vector<string> >  >::iterator it = config_map.begin(); it != config_map.end(); ++it){
    string label = it->first;
    std::vector<string> myvec = it->second;
    int index = label_map[label];
    //    cout << "processing key " << label << endl;
    if(myvec.size() == 0){
      cout << "No parameters given for <" << label << ">" << endl;
      exit(-1);
    }

    //----------------------------------
    if(index == LCREATECLOUD){
      if(myvec.size() != 3){
        SLogger slogger("simparser");
        slogger << ERROR << "wrong label in createcloud " << label << ": Wrong number("<< myvec.size() << ") of parameter!" << SLogger::endmsg;
	exit(-1);
      }else{
	mycreatecloud.flag = true;
	mycreatecloud.cloudsize = atoi(myvec[1].c_str());
	mycreatecloud.constituents = atoi(myvec[2].c_str());
		//mycreatecloud.printme();
      }

      //----------------------------------
    }else if(index == LCLOUDPARTS){
      mycloudparts.flag = true;
      for(int i = 1; i < myvec.size()-1; i = i + 2)
	(mycloudparts.cloudfracs).push_back(pair<string,double>(myvec[i], atof(myvec[i+1].c_str())));
      //      mycloudparts.printme();
      //----------------------------------
    }else if(index == LCLOUDCOORD){
      mycloudcoord.flag = true;
      mycloudcoord.semiaxis_cloud[0] = atof(myvec[1].c_str());
      mycloudcoord.semiaxis_cloud[1] = atof(myvec[2].c_str());
      mycloudcoord.semiaxis_cloud[2] = atof(myvec[3].c_str());
      mycloudcoord.offset_cloud[0] = atof(myvec[4].c_str());
      mycloudcoord.offset_cloud[1] = atof(myvec[5].c_str());
      mycloudcoord.offset_cloud[2] = atof(myvec[6].c_str());
      //mycloudcoord.printme();
      //----------------------------------
    }else if(index == LTEMP){
      mytemp.flag = true;
        for(unsigned int i = 1; i < myvec.size(); i++)
            mytemp.temp.push_back( atof(myvec[i].c_str()));
      // mytemp.printme();
      //----------------------------------
    }else if(index == LCREATEPARTICLES){
      mycreateparticles.flag = true;
      mycreateparticles.nparts = atoi(myvec[1].c_str());
      //mycreateparticles.printme();
    //----------------------------------
    }else if(index == LPARTICLES){
      myparticles.flag = true;
      for(int i = 1; i < myvec.size(); i = i+7){
	std::vector<double> params;
	string pname;
	pname = myvec[i];
	params.push_back(atof(myvec[i+1].c_str()));
	params.push_back(atof(myvec[i+2].c_str()));
	params.push_back(atof(myvec[i+3].c_str()));
	params.push_back(atof(myvec[i+4].c_str()));
	params.push_back(atof(myvec[i+5].c_str()));
	params.push_back(atof(myvec[i+6].c_str()));
	pair<string, std::vector<double> > tmppair; 
	tmppair.first = pname; tmppair.second = params;
	(myparticles.particles).push_back(tmppair);
         //cout << i << " " <<  (myparticles.particles).size() << endl;
      }
       //myparticles.printme();
      //----------------------------------
    }else if(index == LBUFFER){
      mybuffer.flag = true;
      mybuffer.index = atoi(myvec[1].c_str());
      mybuffer.pressure = atof(myvec[2].c_str());
      //       mybuffer.printme();
      //----------------------------------
    }else if(index == LODE){
      myode.flag = true;
      myode.order = atoi(myvec[1].c_str());
      myode.adaptive = atoi(myvec[2].c_str());
      myode.stepsize = atof(myvec[3].c_str());
      //       myode.printme();
      //----------------------------------
    }else if(index == LCOULOMB){
      mycoulomb.flag = true;
      mycoulomb.index = atoi(myvec[1].c_str());
      mycoulomb.weight = atof(myvec[2].c_str());
      //       mycoulomb.printme();
      //----------------------------------
    }else if(index == LOUTPUT){
      myoutput.flag = true;
      myoutput.outfile = myvec[1];
      myoutput.sample_time = atof(myvec[2].c_str());
      myoutput.separate = atoi(myvec[3].c_str());
      //       myoutput.printme();
      //----------------------------------
    }else if(index == LREALTRAP){
      myrealtrap.flag = true;
      myrealtrap.nfiles = atoi(myvec[1].c_str());
      for(int i = 2 ; i < myvec.size(); i++)
	(myrealtrap.filenames).push_back(myvec[i]);
      // myrealtrap.printme();
      //----------------------------------
    }else if(index == LPOTTRAP){
        mypottrap.flag = true;
        mypottrap.potmap_filename = myvec[1];
        mypottrap.Bz = atof(myvec[2].c_str());
        // myrealtrap.printme();
        //----------------------------------
    }else if(index == LIDEALTRAP){
      myidealtrap.flag = true;
      myidealtrap.r_electrode = atof(myvec[1].c_str());
      myidealtrap.Ud2 = atof(myvec[2].c_str());
      myidealtrap.B = atof(myvec[3].c_str());
      //myidealtrap.printme();
      //----------------------------------
    }else if(index == LNONIDEALTRAP){
      mynonitrap.flag = true;
      mynonitrap.trapconfigfile = myvec[1];
      // mynonitrap.printme();
      //----------------------------------
    }else if(index == LCALCTRAP){
      mycalctrap.flag = true;
      mycalctrap.param1 = atof(myvec[1].c_str());
      mycalctrap.param2 = atof(myvec[2].c_str());
      // mycalctrap.printme();
      //----------------------------------
    }else if(index == LIMPORTDATA){
      myimportdata.flag = true;
      myimportdata.prev_simu_file = myvec[1];
      //       myimportdata.printme();
      //----------------------------------
    }
    /*
    else if(index == LNE){
      myne.flag = true;
      myne.time = atof(myvec[1].c_str());
      // myne.printme();
      //----------------------------------
    }else if(index == LDEW){
      mydew.flag = true;
      mydew.size = myvec.size();
      if(myvec.size() == 6){
	mydew.time = atof(myvec[1].c_str());
	mydew.EigenLett = myvec[2];
	mydew.Element = myvec[3];
	mydew.freq_bias = atof(myvec[4].c_str());
	mydew.amp = atof(myvec[5].c_str());
      }else if(myvec.size() == 5){
	mydew.time = atof(myvec[1].c_str());
	mydew.EigenLett = myvec[2];
	mydew.freq_bias = atof(myvec[3].c_str());
	mydew.amp = atof(myvec[4].c_str());
      }
      // mydew.printme();
      //----------------------------------
    }else if(index == LDEB){
      mydeb.flag = true;
      mydeb.size = myvec.size();
      if(myvec.size() == 6){
	mydeb.time = atof(myvec[1].c_str());
	mydeb.EigenLett = myvec[2];
	mydeb.Element = myvec[3];
	mydeb.freq_bias = atof(myvec[4].c_str());
	mydeb.amp = atof(myvec[5].c_str());
      }else if(myvec.size() == 5){
	mydeb.time = atof(myvec[1].c_str());
	mydeb.EigenLett = myvec[2];
	mydeb.freq_bias = atof(myvec[3].c_str());
	mydeb.amp = atof(myvec[4].c_str());
      }
      // mydeb.printme();
      //----------------------------------
    }else if(index == LQEW){
      myqew.flag = true;
      myqew.size = myvec.size();
      if(myvec.size() == 6){
	myqew.time = atof(myvec[1].c_str());
	myqew.EigenLett = myvec[2];
	myqew.Element = myvec[3];
	myqew.freq_bias = atof(myvec[4].c_str());
	myqew.amp = atof(myvec[5].c_str());
      }else if(myvec.size() == 5){
	myqew.time = atof(myvec[1].c_str());
	myqew.EigenLett = myvec[2];
	myqew.freq_bias = atof(myvec[3].c_str());
	myqew.amp = atof(myvec[4].c_str());
      }
      // myqew.printme();
      //----------------------------------
    }else if(index == LQEB){
      myqeb.flag = true;
      myqeb.size = myvec.size();
      if(myvec.size() == 6){
	myqeb.time = atof(myvec[1].c_str());
	myqeb.EigenLett = myvec[2];
	myqeb.Element = myvec[3];
	myqeb.freq_bias = atof(myvec[4].c_str());
	myqeb.amp = atof(myvec[5].c_str());
      }else if(myvec.size() == 5){
	myqeb.time = atof(myvec[1].c_str());
	myqeb.EigenLett = myvec[2];
	myqeb.freq_bias = atof(myvec[3].c_str());
	myqeb.amp = atof(myvec[4].c_str());
      }
      // myqeb.printme();
      //----------------------------------
    }else if(index == LOEW){
      myoew.flag = true;
      myoew.size = myvec.size();
      if(myvec.size() == 6){
	myoew.time = atof(myvec[1].c_str());
	myoew.EigenLett = myvec[2];
	myoew.Element = myvec[3];
	myoew.freq_bias = atof(myvec[4].c_str());
	myoew.amp = atof(myvec[5].c_str());
      }else if(myvec.size() == 5){
	myoew.time = atof(myvec[1].c_str());
	myoew.EigenLett = myvec[2];
	myoew.freq_bias = atof(myvec[3].c_str());
	myoew.amp = atof(myvec[4].c_str());
      }
      //      myoew.printme();
      //----------------------------------
    }else if(index == LOEB){
      myoeb.flag = true;
      myoeb.size = myvec.size();
      if(myvec.size() == 6){
	myoeb.time = atof(myvec[1].c_str());
	myoeb.EigenLett = myvec[2];
	myoeb.Element = myvec[3];
	myoeb.freq_bias = atof(myvec[4].c_str());
	myoeb.amp = atof(myvec[5].c_str());
      }else if(myvec.size() == 5){
	myoeb.time = atof(myvec[1].c_str());
	myoeb.EigenLett = myvec[2];
	myoeb.freq_bias = atof(myvec[3].c_str());
	myoeb.amp = atof(myvec[4].c_str());
      }
      // myoeb.printme();
      //----------------------------------
    }else if(index == LRWW){
      myrww.flag = true;
      myrww.time = atof(myvec[1].c_str());
      myrww.order = atoi(myvec[2].c_str());
      myrww.frequency = atof(myvec[3].c_str());
      myrww.amp = atof(myvec[4].c_str());
      //myrww.printme();
      //----------------------------------
    }else if(index == LRWB){
      myrwb.flag = true;
      myrwb.time = atof(myvec[1].c_str());
      myrwb.order = atoi(myvec[2].c_str());
      myrwb.frequency = atof(myvec[3].c_str());
      myrwb.amp = atof(myvec[4].c_str());
      // myrwb.printme();
      //----------------------------------
    }else if(index == LAWW){
      myaww.flag = true;
      myaww.time = atof(myvec[1].c_str());
      myaww.order = atoi(myvec[2].c_str());
      myaww.frequency = atof(myvec[3].c_str());
      myaww.amp = atof(myvec[4].c_str());
      // myaww.printme();
      //----------------------------------
    }else if(index == LAWB){
      myawb.flag = true;
      myawb.time = atof(myvec[1].c_str());
      myawb.order = atoi(myvec[2].c_str());
      myawb.frequency = atof(myvec[3].c_str());
      myawb.amp = atof(myvec[4].c_str());
      // myawb.printme();
      //----------------------------------
    }else if(index == LAC){
      myac.flag = true;
      myac.time = atof(myvec[1].c_str());
      myac.order = atoi(myvec[2].c_str());
      myac.Element = myvec[3];
      myac.amp = atof(myvec[4].c_str());
      // myac.printme();
      //----------------------------------
    }else if(index == LSC){
      mysc.flag = true;
      if(myvec.size() == 7){
	mysc.time = atof(myvec[1].c_str());
	mysc.order = atoi(myvec[2].c_str());
	mysc.freq1 = atof(myvec[3].c_str());
	mysc.amp1 = atof(myvec[4].c_str());
	mysc.freq2 = atof(myvec[5].c_str());
	mysc.amp2 = atof(myvec[6].c_str());
      }else if(myvec.size() == 8){
	mysc.time = atof(myvec[1].c_str());
	mysc.order = atoi(myvec[2].c_str());
	mysc.Element = myvec[3];
	mysc.freq1 = atof(myvec[4].c_str());
	mysc.amp1 = atof(myvec[5].c_str());
	mysc.freq2 = atof(myvec[6].c_str());
	mysc.amp2 = atof(myvec[7].c_str());
      }
      // mysc.printme();
      //----------------------------------
    }else if(index == LAR){
      myar.flag = true;
      myar.time = atof(myvec[1].c_str());
      myar.amp1 = atof(myvec[2].c_str());
      myar.freq1 = atof(myvec[3].c_str());
      myar.freq2 = atof(myvec[4].c_str());
      //myar.printme();
      //----------------------------------
    }else if(index == LFB){
      myfb.flag = true;
      myfb.time = atof(myvec[1].c_str());
      myfb.amp1 = atof(myvec[2].c_str());
      myfb.freq1 = atof(myvec[3].c_str());
      myfb.freq2 = atof(myvec[4].c_str());
      //myfb.printme();
      //----------------------------------
    }*/
     else if(index == LSWEEP){
        mysweep.flag = true;
        mysweep.wi = atof(myvec[1].c_str());
        mysweep.wf = atof(myvec[2].c_str());
         mysweep.par=0;
         if(myvec.size()==4)
         {
          mysweep.par = atof(myvec[3].c_str());
         }
        //mysweep.printme();
      //----------------------------------  
    }else{
        SLogger slogger("simparser");
        slogger << ERROR << "unknown excitation, label: " << label <<SLogger::endmsg;
        exit(-1);
    }
      
    

  }
  // OPERATIONS
    for(int i=0;i<operation_map.size();i++)
    {
        int index = label_operation_map[operation_map[i][0]];
        std::vector<string> myvec = operation_map[i];
        Operation_sim Ope_tmp;
        Ope_tmp.name = myvec[0];
        if(index == LNE){
            Ope_tmp.time = atof(myvec[1].c_str());
            // myne.printme();
            //----------------------------------
        }else if(index == LDEW){
            Ope_tmp.size = myvec.size();
            if(myvec.size() == 6){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.Element = myvec[3];
                Ope_tmp.freq_bias = atof(myvec[4].c_str());
                Ope_tmp.amp = atof(myvec[5].c_str());
            }else if(myvec.size() == 5){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.freq_bias = atof(myvec[3].c_str());
                Ope_tmp.amp = atof(myvec[4].c_str());
            }
            // mydew.printme();
            //----------------------------------
        }else if(index == LDEB){
            Ope_tmp.size = myvec.size();
            if(myvec.size() == 6){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.Element = myvec[3];
                Ope_tmp.freq_bias = atof(myvec[4].c_str());
                Ope_tmp.amp = atof(myvec[5].c_str());
            }else if(myvec.size() == 5){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.freq_bias = atof(myvec[3].c_str());
                Ope_tmp.amp = atof(myvec[4].c_str());
            }
            // mydeb.printme();
            //----------------------------------
        }else if(index == LQEW){
            Ope_tmp.size = myvec.size();
            if(myvec.size() == 6){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.Element = myvec[3];
                Ope_tmp.freq_bias = atof(myvec[4].c_str());
                Ope_tmp.amp = atof(myvec[5].c_str());
            }else if(myvec.size() == 5){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.freq_bias = atof(myvec[3].c_str());
                Ope_tmp.amp = atof(myvec[4].c_str());
            }
            // myqew.printme();
            //----------------------------------
        }else if(index == LQEB){
            Ope_tmp.size = myvec.size();
            if(myvec.size() == 6){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.Element = myvec[3];
                Ope_tmp.freq_bias = atof(myvec[4].c_str());
                Ope_tmp.amp = atof(myvec[5].c_str());
            }else if(myvec.size() == 5){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.freq_bias = atof(myvec[3].c_str());
                Ope_tmp.amp = atof(myvec[4].c_str());
            }
            // myqeb.printme();
            //----------------------------------
        }else if(index == LOEW){
            Ope_tmp.size = myvec.size();
            if(myvec.size() == 6){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.Element = myvec[3];
                Ope_tmp.freq_bias = atof(myvec[4].c_str());
                Ope_tmp.amp = atof(myvec[5].c_str());
            }else if(myvec.size() == 5){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.freq_bias = atof(myvec[3].c_str());
                Ope_tmp.amp = atof(myvec[4].c_str());
            }
            //      myoew.printme();
            //----------------------------------
        }else if(index == LOEB){
            Ope_tmp.size = myvec.size();
            if(myvec.size() == 6){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.Element = myvec[3];
                Ope_tmp.freq_bias = atof(myvec[4].c_str());
                Ope_tmp.amp = atof(myvec[5].c_str());
            }else if(myvec.size() == 5){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.EigenLett = myvec[2];
                Ope_tmp.freq_bias = atof(myvec[3].c_str());
                Ope_tmp.amp = atof(myvec[4].c_str());
            }
            // myoeb.printme();
            //----------------------------------
        }else if(index == LRWW){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.order = atoi(myvec[2].c_str());
            Ope_tmp.frequency = atof(myvec[3].c_str());
            Ope_tmp.amp = atof(myvec[4].c_str());
            //myrww.printme();
            //----------------------------------
        }else if(index == LRWB){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.order = atoi(myvec[2].c_str());
            Ope_tmp.frequency = atof(myvec[3].c_str());
            Ope_tmp.amp = atof(myvec[4].c_str());
            // myrwb.printme();
            //----------------------------------
        }else if(index == LAWW){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.order = atoi(myvec[2].c_str());
            Ope_tmp.frequency = atof(myvec[3].c_str());
            Ope_tmp.amp = atof(myvec[4].c_str());
            // myaww.printme();
            //----------------------------------
        }else if(index == LAWB){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.order = atoi(myvec[2].c_str());
            Ope_tmp.frequency = atof(myvec[3].c_str());
            Ope_tmp.amp = atof(myvec[4].c_str());
            // myawb.printme();
            //----------------------------------
        }else if(index == LAC){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.order = atoi(myvec[2].c_str());
            Ope_tmp.Element = myvec[3];
            Ope_tmp.amp = atof(myvec[4].c_str());
            // myac.printme();
            //----------------------------------
        }else if(index == LSC){
            if(myvec.size() == 7){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.order = atoi(myvec[2].c_str());
                Ope_tmp.freq1 = atof(myvec[3].c_str());
                Ope_tmp.amp1 = atof(myvec[4].c_str());
                Ope_tmp.freq2 = atof(myvec[5].c_str());
                Ope_tmp.amp2 = atof(myvec[6].c_str());
            }else if(myvec.size() == 8){
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.order = atoi(myvec[2].c_str());
                Ope_tmp.Element = myvec[3];
                Ope_tmp.freq1 = atof(myvec[4].c_str());
                Ope_tmp.amp1 = atof(myvec[5].c_str());
                Ope_tmp.freq2 = atof(myvec[6].c_str());
                Ope_tmp.amp2 = atof(myvec[7].c_str());
            }
            // mysc.printme();
            //----------------------------------
        }else if(index == LAR){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.amp1 = atof(myvec[2].c_str());
            Ope_tmp.freq1 = atof(myvec[3].c_str());
            Ope_tmp.freq2 = atof(myvec[4].c_str());
            //myar.printme();
            //----------------------------------
        }else if(index == LFB){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.amp1 = atof(myvec[2].c_str());
            Ope_tmp.freq1 = atof(myvec[3].c_str());
            Ope_tmp.freq2 = atof(myvec[4].c_str());
            //myfb.printme();
            //----------------------------------
        }else if(index == LSWIFT){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.name_file = myvec[2].c_str();
        }else if(index == LEXC_EMAP){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.name_file = myvec[2].c_str();
            if(myvec.size()==4)
                Ope_tmp.amp = atof(myvec[3].c_str());
            else
                Ope_tmp.amp = 1;
        }else if(index == LPI_PULSE){
            Ope_tmp.time = atof(myvec[1].c_str());
            Ope_tmp.Element = myvec[2];
            
        }else if(index == LASYM_ARW){
            if(myvec.size()==5)
            {
                Ope_tmp.freq_bias = 0;
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.Element = myvec[2];
                Ope_tmp.frequency = atof(myvec[3].c_str());
                Ope_tmp.amp = atof(myvec[4].c_str());
            }
            if(myvec.size()==3)
            {
                Ope_tmp.freq_bias = 1;
                Ope_tmp.time = atof(myvec[1].c_str());
                Ope_tmp.name_file = myvec[2].c_str();
                
            }
        }
        operation_vec.push_back(Ope_tmp);
        
    }
  

}

void SimParser::PrintFile(){
  //print comments
  for(int s = 0; s <comments.size(); s++){
        cout<< comments[s] <<endl;
    }

  //print read fileparts
  for(std::vector< std::pair< string, std::vector< string> > >::iterator it = config_map.begin(); it != config_map.end(); ++it){

    std::vector<string> myvec = it->second;
    string lab = it->first;
    std::cout << lab << " "  << flush;
    for(int i = 1; i < myvec.size(); i++)
      cout << myvec[i] << " " << flush;
    cout << endl;
    
  }
  //print operations
  for(int s = 0; s < operation_map.size();s++){
      std::vector< std::string> myvec = operation_map[s];
     for(int i=0;i<myvec.size();i++){
         cout<<myvec[i]<<" ";
    }
     cout<<" "<<endl;
  }
}

void SimParser::PrintFile(char * filename){
  ofstream outfile;
  outfile.open(filename);
//print comments
for(int s = 0; s <comments.size(); s++){
      outfile<< comments[s] <<endl;
  }

//print read fileparts
for(std::vector< std::pair< string, std::vector< string> > >::iterator it = config_map.begin(); it != config_map.end(); ++it){

  std::vector<string> myvec = it->second;
  string lab = it->first;
  outfile << lab << " "  << flush;
  for(int i = 1; i < myvec.size(); i++)
    outfile << myvec[i] << " " << flush;
  outfile << endl;

}
//print operations
for(int s = 0; s < operation_map.size();s++){
    std::vector< std::string> myvec = operation_map[s];
   for(int i=0;i<myvec.size();i++){
       outfile<<myvec[i]<<" ";
  }
   outfile<<" "<<endl;
}
outfile.close();
}




#ifdef __GUI_ON__

void SimParser::PrintParser(QTextEdit * textwindow){

  textwindow->clear();
//print comments
  for(int s = 0; s <comments.size(); s++){
      textwindow->append(QString::fromStdString(comments[s]));
    }
  textwindow->insertPlainText("\n\n");


//print other information - operations

  if(mycreatecloud.flag){
      textwindow->insertPlainText("CREATECLOUD ");
      textwindow->insertPlainText(QString::number(mycreatecloud.cloudsize));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(mycreatecloud.constituents));
      textwindow->insertPlainText("\n");
  //}
  //if(mycloudparts.flag){
      textwindow->insertPlainText("CLOUDPARTS ");
    for(int i = 0; i < (mycloudparts.cloudfracs).size(); i++){
      textwindow->insertPlainText(QString::fromStdString((mycloudparts.cloudfracs)[i].first));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number((mycloudparts.cloudfracs)[i].second));
    }
    textwindow->insertPlainText("\n");
  //}
  //if(mycloudcoord.flag){
    textwindow->insertPlainText("CLOUDCOORD ");
    textwindow->insertPlainText(QString::number(mycloudcoord.semiaxis_cloud[0]));
    textwindow->insertPlainText(" ");
    textwindow->insertPlainText(QString::number(mycloudcoord.semiaxis_cloud[1]));
    textwindow->insertPlainText(" ");
    textwindow->insertPlainText(QString::number(mycloudcoord.semiaxis_cloud[2]));
    textwindow->insertPlainText(" ");
    textwindow->insertPlainText(QString::number(mycloudcoord.offset_cloud[0]));
    textwindow->insertPlainText(" ");
    textwindow->insertPlainText(QString::number(mycloudcoord.offset_cloud[1]));
    textwindow->insertPlainText(" ");
    textwindow->insertPlainText(QString::number(mycloudcoord.offset_cloud[2]));
    textwindow->insertPlainText("\n");
  //}
  //if(mytemp.flag){
      textwindow->insertPlainText("TEMP ");
    for(int i = 0; i < mytemp.temp.size(); i++){
      textwindow->insertPlainText(QString::number(mytemp.temp[i]));
    }textwindow->insertPlainText("\n");
  }
  if(mycreateparticles.flag){
    textwindow->insertPlainText("CREATEPARTICLES ");
    textwindow->insertPlainText(QString::number(mycreateparticles.nparts));
    textwindow->insertPlainText("\n");
  //}
  //if(myparticles.flag){
    textwindow->insertPlainText("PARTICLES ");
      for(int i = 0; i < (myparticles.particles).size(); i++){
        pair<string, std::vector<double> > tmppair = (myparticles.particles)[i];
        textwindow->insertPlainText(QString::fromStdString(tmppair.first));
        textwindow->insertPlainText(" ");
        for(int j=0; j < tmppair.second.size(); j++){
            textwindow->insertPlainText(QString::number(tmppair.second[j]));
            textwindow->insertPlainText(" ");
          }
    }
    textwindow->insertPlainText("\n");
    //----------------------------------
  }
  if(mybuffer.flag){
      textwindow->insertPlainText("BUFFER ");
      textwindow->insertPlainText(QString::number(mybuffer.index));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(mybuffer.pressure));
      textwindow->insertPlainText("\n");
  }
  if(myode.flag){
      textwindow->insertPlainText("ODE ");
      textwindow->insertPlainText(QString::number(myode.order));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(myode.adaptive));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(myode.stepsize));
      textwindow->insertPlainText("\n");
  }
  if(mycoulomb.flag){
      textwindow->insertPlainText("COULOMB ");
      textwindow->insertPlainText(QString::number(mycoulomb.index));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(mycoulomb.weight));
      textwindow->insertPlainText("\n");
   }
  if(myoutput.flag){
      textwindow->insertPlainText("OUTPUT ");
      textwindow->insertPlainText(QString::fromStdString(myoutput.outfile));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(myoutput.sample_time));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(myoutput.separate));
      textwindow->insertPlainText("\n");
  }
  if(myrealtrap.flag){
      textwindow->insertPlainText("REALTRAP ");
      textwindow->insertPlainText(QString::number(myrealtrap.nfiles));
      textwindow->insertPlainText(" ");
    for(int i = 0 ; i < (myrealtrap.filenames).size(); i++){
        textwindow->insertPlainText(QString::fromStdString((myrealtrap.filenames)[i]));
      textwindow->insertPlainText(" ");
    }
     textwindow->insertPlainText("\n");
  }
  if(mypottrap.flag){
      textwindow->insertPlainText("POTTRAP ");
      textwindow->insertPlainText(QString::fromStdString(mypottrap.potmap_filename));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(mypottrap.Bz));
      textwindow->insertPlainText("\n");
  }
  if(myidealtrap.flag){
      textwindow->insertPlainText("IDEALTRAP ");
      textwindow->insertPlainText(QString::number(myidealtrap.r_electrode));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(myidealtrap.Ud2));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(myidealtrap.B));
      textwindow->insertPlainText("\n");
  }
  if(mynonitrap.flag){
      textwindow->insertPlainText("NONIDEALTRAP ");
      textwindow->insertPlainText(QString::fromStdString(mynonitrap.trapconfigfile));
      textwindow->insertPlainText("\n");
  }
  if(mycalctrap.flag){
      textwindow->insertPlainText("CALCTRAP ");
      textwindow->insertPlainText(QString::number(mycalctrap.param1));
      textwindow->insertPlainText(" ");
      textwindow->insertPlainText(QString::number(mycalctrap.param2));
      textwindow->insertPlainText("\n");
  }
  if(myimportdata.flag){
      textwindow->insertPlainText("IMPORTDATA ");
      textwindow->insertPlainText(QString::fromStdString(myimportdata.prev_simu_file));
      textwindow->insertPlainText("\n");
  }

//print operations
  for(int s = 0; s < operation_map.size();s++){
      std::vector< std::string> myvec = operation_map[s];
     for(int i=0;i<myvec.size();i++){
         textwindow->insertPlainText(QString::fromStdString(myvec[i]+" "));
    }
    textwindow->insertPlainText("\n");
  }
}
#endif

void SimParser::PrintParams(){

  // Print first particle info
  cout  << "-----------------SIMBUCA CONFIGURATION-------------------" << endl;

  if(mycreatecloud.flag){
    cout << "------PARTICLES/CLOUDS: " << endl;
    if(mycloudparts.flag != true || mycloudcoord.flag != true){
      cout << "Cloud parameters or cloud coordinates not specified" << endl;
      exit(-1);
    }else{
      if(mytemp.flag != true){
	cout << "Cloud temperature not specified" << endl;
	exit(-1);
      }else{
	mycreatecloud.printme();
	mycloudparts.printme();
	mycloudcoord.printme();
	mytemp.printme();
      }
    }
  }
  
  if(mycreateparticles.flag){
    if(myparticles.flag != true){
      cout << "Particle info not specified" << endl;
      exit(-1);
    }else{
      myparticles.printme();
    }
  }

  // Print Buffer gas
  if(mybuffer.flag){
    cout << "------BUFFERGAS: " << endl;
    mybuffer.printme();
  }

  // Print ODE config
  if(myode.flag){
    cout << "------ODE: " << endl;
    myode.printme();
  }

  // Print Coulomb config
  if(mycoulomb.flag){
    cout << "------COULOMB: " << endl;
    mycoulomb.printme();
  }

  // Print Output config
  if(myoutput.flag){
    cout << "------OUTPUT: " << endl;
    myoutput.printme();
  } 

  // Print Realtrap config
  if(myrealtrap.flag){
    cout << "------REALTRAP: " << endl;
    myrealtrap.printme();
  } 

  // Print Idealtrap config
  if(myidealtrap.flag){
    cout << "------IDEALTRAP: " << endl;
    myidealtrap.printme();
  } 

  // Print Non-idealtrap config
  if(mynonitrap.flag){
    cout << "------NONIDEALTRAP: " << endl;
    mynonitrap.printme();
  } 

  // Print Importdata config
  if(myimportdata.flag){
    cout << "------IMPORTDATA: " << endl;
    myimportdata.printme();
  }

  // Print Calctrap config
  if(mycalctrap.flag){
    cout << "------CALCTRAP: " << endl;
    mycalctrap.printme();
  }

  // Print excitation configs
  if(myne.flag){
    cout << "------No excitation: " << endl;
    myne.printme();
  }

  if(mydew.flag){
    cout << "------DEW: " << endl;
    mydew.printme();
  }

  if(mydeb.flag){
    cout << "------DEB: " << endl;
    mydeb.printme();
  }

  if(myqew.flag){
    cout << "------QEW: " << endl;
    myqew.printme();
  }

  if(myqeb.flag){
    cout << "------QEB: " << endl;
    myqeb.printme();
  }


  if(myoew.flag){
    cout << "------OEW: " << endl;
    myoew.printme();
  }

  if(myoeb.flag){
    cout << "------OEB: " << endl;
    myoeb.printme();
  }

  if(myrww.flag){
    cout << "------RWW: " << endl;
    myrww.printme();
  }

  if(myrwb.flag){
    cout << "------RWB: " << endl;
    myrwb.printme();
  }


  if(myaww.flag){
    cout << "------AWW: " << endl;
    myaww.printme();
  }

  if(myawb.flag){
    cout << "------AWB: " << endl;
    myawb.printme();
  }

  if(myac.flag){
    cout << "------AC: " << endl;
    myac.printme();
  }

  if(mysc.flag){
    cout << "------SC: " << endl;
    mysc.printme();
  }
  
  if(myar.flag){
    cout << "------AR: " << endl;
    myar.printme();
  }
  
if(myfb.flag){
    cout << "------AR: " << endl;
    myfb.printme();
  }

if(mysweep.flag){
    cout << "------SWEEP: " << endl;
    mysweep.printme();
  }
  cout << "---------------------------------------------------------" << endl;

}

void SimParser::SetTrapConfigOff(){
    myrealtrap.flag = false;
    mypottrap.flag = false;
    myidealtrap.flag = false;
    mynonitrap.flag  = false;
    mycalctrap.flag = false;
}
