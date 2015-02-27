//
//  PDGTable.h
//
//
//  Created by balint radics on 08/02/13.
//
//

#include <stdlib.h>
#include "globals.h"
#include "PDGTable.h"

//______________________________________________________________________________
PDGTable::PDGTable()
{
  pVerbose = false;
}

//______________________________________________________________________________
PDGTable::PDGTable(string filename_table_pdg)
{
  pVerbose = false;
  this->ReadPDGTable(filename_table_pdg);
}

//______________________________________________________________________________
PDGTable::~PDGTable()
{

}


void PDGTable::MakeVerbose(bool verbose){
  pVerbose = verbose;
}

//______________________________________________________________________________
void PDGTable::Print() 
{

  vector<string>::iterator itN = pName.begin();
  vector<int>::iterator itId = pId.begin();
  vector<double>::iterator itMass = pMass.begin();
  vector<double>::iterator itCharge = pCharge.begin();
  for(; itN != pName.end() && itId != pId.end() && itMass != pMass.end() && itCharge != pCharge.end(); ++itN, ++itId, ++itMass, ++itCharge){
    cout << *itN << "\t" << *itId << "\t" << *itMass << " GeV \t" << *itCharge << endl;
  }

}

//______________________________________________________________________________
void PDGTable::Print(int pdgid) 
{

  vector<string>::iterator itN = pName.begin();
  vector<int>::iterator itId = pId.begin();
  vector<double>::iterator itMass = pMass.begin();
  vector<double>::iterator itCharge = pCharge.begin();
  for(; itN != pName.end() && itId != pId.end() && itMass != pMass.end() && itCharge != pCharge.end(); ++itN, ++itId, ++itMass, ++itCharge){
    if(*itId == pdgid){
      cout << *itN << "\t" << *itId << "\t" << *itMass << " GeV \t" << *itCharge << endl;
      break;
    }
  }


}

//______________________________________________________________________________
void PDGTable::Print(vector<int> pdgids) 
{

  for(vector<int>::iterator it = pdgids.begin(); it != pdgids.end(); ++it)
    this->Print(*it);

}

//______________________________________________________________________________
void PDGTable::ReadPDGTable(string filename_table_pdg)
{
   // read list of particles from a file
   // See libraries/pdg/pdg_table.txt to see the file format

   FILE* file = fopen(filename_table_pdg.c_str(),"r");
   if (file == 0) {
     SLogger slogger("trapparameters");
     slogger << ERROR << "PDGTable: Could not open PDG particle file " << filename_table_pdg.c_str() << SLogger::endmsg;
     exit(-1);
     return;
   }

   char      c[512];
   int       class_number, anti, isospin, i3, spin, tracking_code;
   int       ich, kf, nch, charge;
   char      name[30], class_name[30];
   double    mass, width, branching_ratio;
   int       dau[20];

   int       idecay, decay_type, flavor, ndau, stable;

   int       input;
   while ( (input=getc(file)) != EOF) {
      c[0] = input;
      if (c[0] != '#') {
         ungetc(c[0],file);
         // read channel number
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%i",&ich)) {;}
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%s",name  )) {;}
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%i",&kf   )) {;}
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%i",&anti )) {;}

         if (kf < 0) {
	   //            AddAntiParticle(name,kf);
	   // assuming particle has been already defined
	   int pdgid = abs(kf);
	   pName.push_back(name);
	   pMass.push_back(GetPDGMass(pdgid));
	   pCharge.push_back(-GetPDGCharge(pdgid));
	   pId.push_back(kf);

            // nothing more on this line
            if (fgets(c,200,file)) {;}
         } else {
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&class_number)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%s",class_name)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&charge)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%le",&mass)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%le",&width)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&isospin)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&i3)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&spin)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&flavor)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&tracking_code)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&nch)) {;}
            // nothing more on this line
            if (fgets(c,200,file)) {;}
            if (width > 1e-10) stable = 0;
            else               stable = 1;

            // add particle to vectors
	    pName.push_back(string(name));
	    pMass.push_back(mass);
	    pId.push_back(kf);
	    pCharge.push_back(charge/3.0);
	    

            if (nch) {
               // read in decay channels
               ich = 0;
               int c_input = 0;
               while ( ((c_input=getc(file)) != EOF) && (ich <nch)) {
                  c[0] = c_input;
                  if (c[0] != '#') {
                     ungetc(c[0],file);

                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%i",&idecay)) {;}
                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%i",&decay_type)) {;}
                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%le",&branching_ratio)) {;}
                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%i",&ndau)) {;}
                     for (int idau=0; idau<ndau; idau++) {
                        // coverity [secure_coding : FALSE]
                        if (fscanf(file,"%i",&dau[idau])) {;}
                     }
                     // add decay channel

		     //                     if (part) part->AddDecayChannel(decay_type,branching_ratio,ndau,dau);
                     ich++;
                  }
                  // skip end of line
                  if (fgets(c,200,file)) {;}
               }
            }
         }
      } else {
         // skip end of line
         if (fgets(c,200,file)) {;}
      }
   }
   // // in the end loop over the antiparticles and
   // // define their decay lists
   // TIter it(fParticleList);

   // Int_t code[20];
   // TParticlePDG  *ap, *p, *daughter;
   // TDecayChannel *dc;

   // while ((p = (TParticlePDG*) it.Next())) {

   //    // define decay channels for antiparticles
   //    if (p->PdgCode() < 0) {
   //       ap = GetParticle(-p->PdgCode());
   //       if (!ap) continue;
   //       nch = ap->NDecayChannels();
   //       for (ich=0; ich<nch; ich++) {
   //          dc = ap->DecayChannel(ich);
   //          if (!dc) continue;
   //          ndau = dc->NDaughters();
   //          for (int i=0; i<ndau; i++) {
   //             // conserve CPT

   //             code[i] = dc->DaughterPdgCode(i);
   //             daughter = GetParticle(code[i]);
   //             if (daughter && daughter->AntiParticle()) {
   //                // this particle does have an
   //                // antiparticle
   //                code[i] = -code[i];
   //             }
   //          }
   //          p->AddDecayChannel(dc->MatrixElementCode(),
   //                             dc->BranchingRatio(),
   //                             dc->NDaughters(),
   //                             code);
   //       }
   //       p->SetAntiParticle(ap);
   //       ap->SetAntiParticle(p);
   //    }
   // }

   fclose(file);
   return;
}

//______________________________________________________________________________
string PDGTable::GetPDGName(int pdgid){

  vector<string>::iterator itN = pName.begin();
  vector<int>::iterator itId = pId.begin();
  string ret = "";
  bool found = false;
  for(; itN != pName.end() && itId != pId.end(); ++itN, ++itId){
    if(*itId == pdgid){
      ret = *itN;found = true;
      break;
    }
  }  

  if(found){
    return ret;
  }else{
    if(pVerbose)cerr << "PDG Id [" << pdgid << "] not found!" << endl;
    return ret;
  }

  return ret;
}

//______________________________________________________________________________
int PDGTable::GetPDGId(string name){

  vector<string>::iterator itN = pName.begin();
  vector<int>::iterator itId = pId.begin();
  int ret = 123456789;
  bool found = false;
  for(; itN != pName.end() && itId != pId.end(); ++itN, ++itId){
    if(*itN == name){
      ret = *itId;found = true;
      break;
    }
  }  

  if(found){
    return ret;
  }else{
    if(pVerbose)cerr << "PDG Name [" << name << "] not found!" << endl;
    return ret;
  }

  return ret;
}

//______________________________________________________________________________
double PDGTable::GetPDGMass(int pdgid){
  vector<double>::iterator itMass = pMass.begin();
  vector<int>::iterator itId = pId.begin();
  double ret = -99999;
  bool found = false;
  for(; itMass != pMass.end() && itId != pId.end(); ++itMass, ++itId){
    if(*itId == pdgid){
      ret = *itMass;found = true;
      break;
    }
  }  

  if(found){
    return ret;
  }else{
    if(pVerbose)cerr << "PDG Id [" << pdgid << "] not found!" << endl;
    return ret;
  }

  return ret;

}

//______________________________________________________________________________
double PDGTable::GetPDGCharge(int pdgid){

  vector<double>::iterator itC = pCharge.begin();
  vector<int>::iterator itId = pId.begin();
  double ret = -99999999;
  bool found = false;
  for(; itC != pCharge.end() && itId != pId.end(); ++itC, ++itId){
    if(*itId == pdgid){
      ret = *itC;found = true;
      break;
    }
  }  

  if(found){
    return ret;
  }else{
    if(pVerbose)cerr << "PDG Id [" << pdgid << "] not found!" << endl;
    return ret;
  }

  return ret;

}

//______________________________________________________________________________
bool PDGTable::TestPDGName(string name_){

  int pdgid = this->GetPDGId(name_);
  if(pVerbose)cout << "Testing : " << name_ << flush;
  if(pdgid == 123456789){
    if(pVerbose)cout << ".. not exist" << endl;
    return false;
  }else{
    if(pVerbose)cout << ".. exist" << endl;
    return true;
  }
  return false;

}

//______________________________________________________________________________
bool PDGTable::TestPDGNames(vector<string > names_){
  
  bool res = true;
  for(vector<string>::iterator it = names_.begin(); it != names_.end(); ++it)
    res &= TestPDGName(*it);

  return res;
}
