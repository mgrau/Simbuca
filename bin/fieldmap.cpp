//fieldmap.cpp

#include "fieldmap.h"
#include <cstdlib> // for exit function
#include <sstream>
#include "math.h"
#include <stdexcept>

#include <fstream>
using std::ifstream;
using std::ofstream;


//int round(&double num) {
//  char sign = static_cast<char>(num/fabs(num));
//  return static_cast<int>(sign* (fabs(num)+ 0.5));
//}    
 
 
vector<double> field(3);

 
inline double fieldmap::round(double d){
  return floor(d + 0.5);
}


 vector<double> values;
 
 
bool fieldmap::IsElement(double _value){

  bool returnvalue = false;

  vector<int>::reverse_iterator rit;
  for ( int index= 0 ; index < values.size(); index++ ){
   	if ( values[index] == _value) returnvalue = true;
  }
  
  
  
return returnvalue;
 }  
    

//fieldmap::fieldmap(){
//}


fieldmap::fieldmap(int _nrows, int _ncols){
	interpolatepoints = true;
	nrows = _nrows;
	ncols = _ncols;
	matrix = new Matrix<std::pair<double,double> >(_nrows,_ncols);
    factor = 1.;
}

fieldmap::~fieldmap(){
	delete matrix;
}

void fieldmap::ReadField(char * _Frz_filename){
     //Reads from a file constructed as r z Er Ez separated by tabs or spaces.
     //For each r stepped, all z-stepped are being walked through.
     //steps over z-values first
     //should start with z-values from low to high values    
    ifstream Frz_str;
    
    Frz_str.open(_Frz_filename); 
    if (!Frz_str){
        SLogger slogger("fieldmap");
        slogger << ERROR << "opening "<<_Frz_filename<<". The program will be stopped." << SLogger::endmsg;exit(1);
      }
     
     double r,z,Fz,Fr,Fmod;
     double r2,z2;
     int r_stepper = 0, z_stepper = 0;
     char firstline[256]; 
     char comment= '%';
     while (Frz_str.peek() == comment){
     	Frz_str.getline(firstline, 256);
     }
     
     //read file minus first whiteline.
     bool valuesset = false;

     while ( !Frz_str.eof() && z_stepper != nrows && r_stepper != ncols ) { // keep reading until end-of-file 
        Frz_str>>r>>z>>Fr>>Fz;
        //Br_min,max en voor z ook bepalen he.
        if(!valuesset){valuesset=true;r_min=r_max=r;z_min=z_max=z;}
        if(r<r_min){r_min=r;}if(r>r_max){r_max=r;}if(z<z_min){z_min=z;} if(z>z_max){z_max=z;}        
        (*matrix)(z_stepper,r_stepper).first=Fz;
        (*matrix)(z_stepper,r_stepper).second=Fr;
        z_stepper++;
        if (z_stepper == nrows){z_stepper =0;r_stepper++;}                     
     } 
     
     
     
     //set the r and z_step
     r_step=(r_max-r_min)/(ncols-1);
     z_step=(z_max-z_min)/(nrows-1);
     
     cout<<"\t R: "<<r_min<<" "<<r_max<<" "<<r_step<<endl;
     cout<<"\t z: "<<z_min<<" "<<z_max<<" "<<z_step<<endl;
     /*
     for (unsigned i=0; i < Ematrix.nrows(); i++){
         for (unsigned j=0; j < Ematrix.ncols();j++){
             cout<<Ez_min+i*Ez_step<<" \t"<<Er_min+j*Er_step<<" \t"<<Ematrix(i,j).first<<endl;
             }
             cin.get();
     }    
    */
    cout<<"file: "<<_Frz_filename<<" read\n"; 
    Frz_str.close();               
}    
    
    
/*read function*/    
void fieldmap::ReadField(char *_Fr_filename,char *_Fz_filename){
     //Supposed it that the 2 files have the same number of indexes!!!
    ifstream Fr_str, Fz_str;
    
    Fr_str.open(_Fr_filename); 
    Fz_str.open(_Fz_filename);
    if (!Fr_str){
        SLogger slogger("fieldmap");
        slogger << ERROR << "opening "<<_Fr_filename<<". The program will be stopped." << SLogger::endmsg;exit(1);
      }
    if (!Fr_str){
        SLogger slogger("fieldmap");
        slogger << ERROR << "opening "<<_Fr_filename<<". The program will be stopped." << SLogger::endmsg;exit(1);
      }
     
     double r,z,Fz,Fr,Fmod;
     double r2,z2;
     int r_stepper = 0, z_stepper = 0;
     char firstline[256]; 
     char comment= '%';
     while (Fr_str.peek() == comment){
     	Fr_str.getline(firstline, 256);
     }
     while (Fz_str.peek() == comment){
     	Fz_str.getline(firstline, 256);
     }
     
     //read file minus first whiteline.
     bool valuesset = false;


     while ( !Fr_str.eof() && z_stepper != nrows && r_stepper != ncols ) { // keep reading until end-of-file 
        Fr_str>>r>>z>>Fr;
        Fz_str>>r>>z>>Fz;
        //Br_min,max en voor z ook bepalen he.
        if(!valuesset){valuesset=true;r_min=r_max=r;z_min=z_max=z;}
        if(r<r_min){r_min=r;}if(r>r_max){r_max=r;}if(z<z_min){z_min=z;} if(z>z_max){z_max=z;}        
        (*matrix)(z_stepper,r_stepper).first=Fz;
        (*matrix)(z_stepper,r_stepper).second=Fr;
        z_stepper++;
        if (z_stepper == nrows){z_stepper =0;r_stepper++;}                     
     }
     
     //set the r and z_step
     r_step=(r_max-r_min)/(ncols-1);
     z_step=(z_max-z_min)/(nrows-1);
     
     //cout<<"\t R: "<<r_min<<" "<<r_max<<" "<<r_step<<endl;
     //cout<<"\t z: "<<z_min<<" "<<z_max<<" "<<z_step<<endl;
     /*
     for (unsigned i=0; i < Ematrix.nrows(); i++){
         for (unsigned j=0; j < Ematrix.ncols();j++){
             cout<<Ez_min+i*Ez_step<<" \t"<<Er_min+j*Er_step<<" \t"<<Ematrix(i,j).first<<endl;
             }
             cin.get();
     }    
    */
    cout<<"file: "<<_Fr_filename<<" read\n";
    cout<<"file: "<<_Fz_filename<<" read\n";    
    Fr_str.close(); 
    Fz_str.close();               
}

  
/*Interpolate functions*/  
vector<double> fieldmap::getField(const double& _x,const double& _y,const double& _z){

     r=sqrt(_x*_x+_y*_y);if(r == 0)r=1e-10; //otherwise division by zero!     
     //check if particles are within the potential
     if(_x>r_max||_y>r_max||r>r_max ||_z>z_max||_z<(z_min+r_step)){
           cout<<"Particle position: "<<_x<<" "<<_y<<" "<<_z<<" not found in electric fieldmap"<<endl;}
       	   //scale both particles
      	   z_index= _z-z_min; 
      	   r_index= r-r_min;


     if(interpolatepoints){
      	   	r1_index=fabs(floor(r_index/r_step));r2_index=r1_index+1;//ceil(r_index/Er_step); //the next one he: (r+0.5*r_step)/r_step
      		z1_index=floor(z_index/z_step);z2_index=z1_index+1;//ceil(z_index/Ez_step);
        	//matrix:first row(z), then column(r)
      
        	m12=(*matrix)(z2_index,r1_index);m22=(*matrix)(z2_index,r2_index);
        	m11=(*matrix)(z1_index,r1_index);m21=(*matrix)(z1_index,r2_index);
      
      		r1=r1_index*r_step+r_min;r2=r2_index*r_step+r_min;
      		z1=z1_index*z_step+z_min;z2=z2_index*z_step+z_min;
      
   
      		//interpolate Potential!  http://en.wikipedia.org/wiki/Bilinear_interpolation
      		dzr=1.0/((r2-r1)*(z2-z1));
      		mzr.first = m11.first*(r2-r)*(z2-_z)*dzr
                 +m21.first*(r-r1)*(z2-_z)*dzr
                 +m12.first*(r2-r)*(_z-z1)*dzr
                 +m22.first*(r-r1)*(_z-z1)*dzr;
      		mzr.second = m11.second*(r2-r)*(z2-_z)*dzr
                 +m21.second*(r-r1)*(z2-_z)*dzr
                 +m12.second*(r2-r)*(_z-z1)*dzr
                 +m22.second*(r-r1)*(_z-z1)*dzr;
      		//if nan, solve it this way he!           
      		if(mzr.second != mzr.second){ 
            		mzr.second = 0.0;
                        SLogger slogger("fieldmap");
                        slogger << ERROR << "fieldmap error\t x: "<<_x<<" y:\t"<<_y<<" r:\t"<<r<<" z:\t"<<_z << SLogger::endmsg;
            		cout<<"fieldmap error: "<<dzr<<endl;
            		cout<<"x: "<<_x<<endl;
            		cout<<"y: "<<_y<<endl;
           		cout<<"r: "<<r<<endl;
            		cout<<"z: "<<_z<<endl; 
      			cout<<r1_index<<"\t"<<r2_index<<endl;
      			cout<<z1_index<<"\t"<<z2_index<<endl;            
            		field[0]=0.0;//E11.second;
            		field[1]=0.0;
      		}
      		else{
           		/*if (fabs(_x) < 0.0001)*/ field[0]=mzr.second*(_x/r); //else electric_field[0]=0.0;
           		/*if (fabs(_y) < 0.0001)*/ field[1]=mzr.second*(_y/r); //else electric_field[0]=0.0;
      		}
      
      		if(mzr.first != mzr.first){ 
           		mzr.first=0.0;
           		field[2]=0.0;                                  
      		}
      		else{
          		 field[2]=mzr.first;           
     		}
     }
     else{ // no interpolation = take neirest neighbour point
      	r1_index=fabs(round(r_index/r_step));
      	z1_index=(round(z_index/z_step));

	m12=(*matrix)(z1_index,r1_index);
        field[0]=m12.second*(_x/r);
        field[1]=m12.second*(_y/r);
        field[2]=m12.first;

    	//cout<<r1_index<<"\t"<<z1_index<<"\t"<<B12.first<<"\t"<<B12.second<<endl;
	//cin.get();
     }
     
       //cout<<"asked for E with z="<<_z<<" r="<<r<<"  ->Ez="<<Ezr.first<<" // Er="<<Ezr.second<<" V/m"<<endl;
     return field;                
}     

vector<double> fieldmap::getFieldXY(const double& _x,const double& _y){
    
    //r=sqrt(_x*_x+_y*_y);if(r == 0)r=1e-10; //otherwise division by zero!
    //check if particles are within the potential
    r = _x;
    if(_x>r_max||_y>z_max||_y<(r_min+r_step)||_y<(z_min+r_step)){
        cout<<"XY Particle position: "<<_x<<" "<<_y<<" not found in electric fieldmap"<<endl;
        field[0]=0.0;//E11.second;
        field[1]=0.0;
        return field;
    }
    //scale both particles
    z_index= _y-z_min;
    r_index= _x-r_min;
    
    
    if(interpolatepoints){
        r1_index=fabs(floor(r_index/r_step));r2_index=r1_index+1;//ceil(r_index/Er_step); //the next one he: (r+0.5*r_step)/r_step
        z1_index=floor(z_index/z_step);z2_index=z1_index+1;//ceil(z_index/Ez_step);
        //matrix:first row(z), then column(r)
        
        m12=(*matrix)(z2_index,r1_index);m22=(*matrix)(z2_index,r2_index);
        m11=(*matrix)(z1_index,r1_index);m21=(*matrix)(z1_index,r2_index);
        
        r1=r1_index*r_step+r_min;r2=r2_index*r_step+r_min;
        z1=z1_index*z_step+z_min;z2=z2_index*z_step+z_min;
        
        
        //interpolate Potential!  http://en.wikipedia.org/wiki/Bilinear_interpolation
        dzr=1.0/((r2-r1)*(z2-z1));
        mzr.first = m11.first*(r2-r)*(z2-_y)*dzr
        +m21.first*(r-r1)*(z2-_y)*dzr
        +m12.first*(r2-r)*(_y-z1)*dzr
        +m22.first*(r-r1)*(_y-z1)*dzr;
        mzr.second = m11.second*(r2-r)*(z2-_y)*dzr
        +m21.second*(r-r1)*(z2-_y)*dzr
        +m12.second*(r2-r)*(_y-z1)*dzr
        +m22.second*(r-r1)*(_y-z1)*dzr;
        //if nan, solve it this way he!
        if(mzr.second != mzr.second){
            mzr.second = 0.0;
            SLogger slogger("fieldmap");
            slogger << ERROR << "fieldmap error\t x: "<<_x<<" y:\t"<<_y<<" r:\t"<<r<< SLogger::endmsg;
            cout<<"fieldmap error: "<<dzr<<endl;
            cout<<"x: "<<_x<<endl;
            cout<<"y: "<<_y<<endl;
            //cout<<"r: "<<r<<endl;
            //cout<<"z: "<<_z<<endl;
            cout<<r1_index<<"\t"<<r2_index<<endl;
            cout<<z1_index<<"\t"<<z2_index<<endl;
            cin.get();
            field[0]=0.0;//E11.second;
            field[1]=0.0;
        }
        else{
            /*if (fabs(_x) < 0.0001)*/ field[0]=factor*mzr.second; //else electric_field[0]=0.0;
            /*if (fabs(_y) < 0.0001)*/ //field[1]=factor*mzr.second*(_y/r); //else electric_field[0]=0.0;
        }
        
        if(mzr.first != mzr.first){
            mzr.first=0.0;
            field[1]=0.0;
        }
        else{
            field[1]=factor*mzr.first;
        }
    }
    else{ // no interpolation = take neirest neighbour point
      	r1_index=fabs(round(r_index/r_step));
      	z1_index=(round(z_index/z_step));
        
        m12=(*matrix)(z1_index,r1_index);
        field[0]=factor*m12.second;
       // field[1]=m12.second*(_y/r);
        field[1]=factor*m12.first;
        
    	//cout<<r1_index<<"\t"<<z1_index<<"\t"<<B12.first<<"\t"<<B12.second<<endl;
        //cin.get();
    }
    
    //cout<<"asked for E with z="<<_z<<" r="<<r<<"  ->Ez="<<Ezr.first<<" // Er="<<Ezr.second<<" V/m"<<endl;
    return field;
}
void fieldmap::Print(){
	for(unsigned int i=0; i < nrows; i++){
		for(unsigned int j=0;j < ncols; j++){
			cout<<(*matrix)(i,j).first<<" ";
		}
		cout<<endl;
	}
}

void fieldmap::PrintOnAxis(){
	for(unsigned int i=0;i<nrows;i++){
		cout<<i<<"\t"<<(*matrix)(i,0).first<<"\t"<<(*matrix)(i,0).second<<endl;
	}
	
}

void fieldmap::SetInterpolate(double _interpolate){ 
	// interpolate the fielmap points yes or no...
       // standard the interpolation is ON
	interpolatepoints = _interpolate;
}

/*    
void ReadVoltages(char *_V_filename){
        
     ifstream V_str;
     V_str.open(_V_filename); 
     if (!V_str){cout<<"Error in opening "<<_V_filename<<". The program will be stopped."<<endl; exit(1);}
     
     double r,z,V;
     int r_stepper = 0, z_stepper = 0;

     char firstline[256]; V_str.getline(firstline, 256);
     //read file minus first whiteline.
     bool valuesset = false;
     while ( !V_str.eof() && z_stepper != Vmatrix.nrows() && r_stepper != Vmatrix.ncols() ) { // keep reading until end-of-file
        V_str>>r>>z>>V;
        //r_min,max en voor z ook bepalen he.
        if(!valuesset){valuesset=true;r_min=r_max=r;z_min=z_max=z;}
        if(r<r_min){r_min=r;}
        if(r>r_max){r_max=r;}
        if(z<z_min){z_min=z;}
        if(z>z_max){z_max=z;}        
        Vmatrix(z_stepper,r_stepper)=V;
        z_stepper++;
        if (z_stepper == Vmatrix.nrows()){z_stepper =0;r_stepper++;}                      
     }
     //set the r and z_step
     r_step=(r_max-r_min)/(Vmatrix.ncols()-1);
     z_step=(z_max-z_min)/(Vmatrix.nrows()-1);
    
    //  cout<<r_min<<" "<<r_max<<" "<<r_step<<endl;
    //  cout<<z_min<<" "<<z_max<<" "<<z_step<<endl;
     
     cout<<"file: "<<_V_filename<<" read\n";
    // cout<<"*-*"<<Vmatrix(55,55)<<"*-*\n";
     V_str.close(); 
} 
*/
  
  
  
/*
vector<double>& getElectr_via_voltage(const double& _x,const double& _y,const double& _z){
     double x,y,z;
     x = _x;y = _y;z= _z;
     r=sqrt(x*x+y*y);if(r == 0)r=1e-10; //otherwise division by zero!     
   
     
     //check of ze wel in of bound zitten he!
  //   if(x>r_max||y>r_max||z>z_max||x<r_min||y<r_min||z<z_min){
  //         cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;}
     //scale both particles he.
     //uppercorner is(z_min,r_min)
      z= z+(-z_min); 
      r= r+(-r_min);

     
     //find correct neighbours in array
      double electric_fieldr, electric_fieldz;
      double rs[4],zs[4];
      double volts[4];
      int coord[4];
      
      rs[0] = fabs(r - r_step*0.5);
      rs[1] = r + r_step*0.5;
      rs[2] = rs[3] = r;
      
      zs[2] = z + z_step*0.5;
      zs[3] = z - z_step*0.5;
      zs[0] = zs[1] = z;
      
      double u = 0, t = 0;
      
      for (int i = 0;i<4;i++){
          coord[0] = (int)(rs[i]/r_step);
          coord[1] = (int)(rs[i]/r_step)+1;
          coord[2] = (int)(zs[i]/z_step);
          coord[3] = (int)(zs[i]/z_step)+1;
      
          t = (rs[i]-coord[0]*r_step)/r_step;
          u = (zs[i]-coord[2]*z_step)/z_step;
          
         if((coord[1] >= Vmatrix.ncols())|| (coord[3]>=Vmatrix.nrows()) || coord[2]<0.0 || coord[0]<0.0){
                cout<<"Array at of bound error in 'void getElectr(double x,double y,double z)'..."<<endl;
                cout<<"rs[i] = "<<rs[i]<<endl;
                cout<<"coord[1] = "<<coord[1]<<endl;
                cout<<"zs[i] = "<<zs[i]<<endl;
                cout<<"coord[3] = "<<coord[3]<<endl;
                cout<<"Program terminated"<<endl;   
                exit(1); 
                break;  
          }
          volts[i] = (1-t)*(1-u)*Vmatrix(coord[2],coord[0])+
                t*(1-u)*Vmatrix(coord[2],coord[1])+
                t*u*Vmatrix(coord[3],coord[1])+
                (1-t)*u*Vmatrix(coord[3],coord[0]);
      }
      electric_fieldr = (volts[0]-volts[1])/r_step;
      electric_fieldz = (volts[3]-volts[2])/z_step;


      
      
      electric_field[0] = electric_fieldr*(x/r);
      electric_field[1] = electric_fieldr*(y/r);
      electric_field[2] = electric_fieldz;
      return electric_field;
 
}
*/


/*
vector<double>& getElectrInterpolation(const double& _x,const double& _y,const double& _z){
     //Interpolate potentials on r,z direction and calculate stuff
     //Currently works only if r>0 and z>0.
     
     double x,y,z,r;
     x = _x;y = _y;z= _z;
     r=sqrt(x*x+y*y);if(r == 0)r=1e-10; //otherwise division by zero!     
   
     
     //check of ze wel in of bound zitten he!
     if(x>r_max||y>r_max||z>z_max||z<z_min){
           cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;cin.get();}
     //scale both particles
      double z_index= z+(-z_min); 
      double r_index= r+(-r_min);

      double z1,z2,r1,r2;
      double z1_index,z2_index,r1_index,r2_index;

      r1_index=(int)(r_index/r_step-0.5);r2_index=(int)(r_index/r_step+0.5); //the next one he: (r+0.5*r_step)/r_step
      z1_index=(int)(z_index/z_step-0.5);z2_index=(int)(z_index/z_step+0.5);
      
      //matrix is first row, then column he.
      double V12=Vmatrix(z2_index,r1_index);double V22=Vmatrix(z2_index,r2_index);
      double V11=Vmatrix(z1_index,r1_index);double V21=Vmatrix(z1_index,r2_index);
      
      r1=r1_index*r_step+r_min;r2=r2_index*r_step+r_min;
      z1=z1_index*z_step+z_min;z2=z2_index*z_step+z_min;
      
   
      //interpolate Potential!  http://en.wikipedia.org/wiki/Bilinear_interpolation
      double dz=z2-z1;
      double dr=r2-r1;
      double dzr=dr*dz;
      double Vzr= V11*(r2-r)*(z2-z)/dzr
                 +V21*(r-r1)*(z2-z)/dzr
                 +V12*(r2-r)*(z-z1)/dzr
                 +V22*(r-r1)*(z-z1)/dzr;

      cout.precision(8);*/
      /*   r1,z2----------Vup----r2,z2
            |              |      |      
            |              |      |
          Vleft-----------Vrz----Vright                        
            |              |      |
            |              |      |
          r1,z1-----------Vdown--r2,z1        */
/*      double Vleft,Vright,Vup,Vdown;
      //interpolate potentials  
      Vleft =V11+(z-z1)*(V12-V11)/(z2-z1);
      Vright=V21+(z-z1)*(V22-V21)/(z2-z1);      
      Vup   =V12+(r-r1)*(V22-V12)/(r2-r1);
      Vdown =V11+(r-r1)*(V21-V11)/(r2-r1);
     //Calculate Electric Field left and right
      double Erleft=-(Vzr-Vleft)/(r-r1);
      double Erright=-(Vright-Vzr)/(r2-r);
      double Ezup=-(Vzr-Vup)/(z-z2);
      double Ezdown=-(Vdown-Vzr)/(z1-z);
      //Interpolate Electric Field
      double Er=Erleft+(r-r1)*(Erright-Erleft)/(r2-r1);
      double Ez=Ezup+(z-z2)*(Ezdown-Ezup)/(z1-z2);

      
 */     //methode a2) of via interpolatie Voltages op de zijassen berekenen, en dan enkel dichtste bij punten nemen!
   /*   Er=(Vright-Vzr)/(r2-r);
      Ez=(Vup-Vzr)/(z2-z);
   */   
      //methode b) of gewoon meer gemiddeldes nemen op de hoekpunten van omringen vierkand
   /*   
    */  
      //methode c) gewoon via hoekpunten, redelijk efficient!!
  /*    double Er=(-(V22-V12)/r_step-(V21-V11)/r_step)*0.5;
      double Ez=(-(V12-V11)/z_step-(V22-V21)/z_step)*0.5;
      double V=Vzr;   
   */ 
      
 /*     electric_field[0]=Er*(x/r);
      electric_field[1]=Er*(y/r);
      electric_field[2]=Ez;
      return electric_field;
}*/

/*
vector<double>& getElectrFast(const double& _x,const double& _y,const double& _z){
     //Interpolate potential+get potential on closest point.
     //Calculate Electric field due to this 2 potentials.
     //Not yet implemented
     
     double x,y,z,r;
     x = _x;y = _y;z= _z;
     r=sqrt(x*x+y*y);if(r == 0)r=1e-10; //otherwise division by zero!     
   
     
     //check of ze wel in of bound zitten he!
     if(x>r_max||y>r_max||z>z_max||z<z_min){
           cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;cin.get();}
     //scale both particles
      double z_index= z+(-z_min); 
      double r_index= r+(-r_min);

      double z1,z2,r1,r2;
      double z1_index,z2_index,r1_index,r2_index;

      //to be implemented
      
  //    electric_field[0]=Er*(x/r);
  //    electric_field[1]=Er*(y/r);
  //    electric_field[2]=Ez;
      return electric_field;
}

*/


/*

      Very old function

*/

/*
void ConvertBfile(char *_B_filename){
    ofstream testie;
    testie.open("9T01T_converted.txt"); 

    int x_scale_mags = 49, z_scale_mags = 1161;
    double grid_scale_mag = 0.0025;//0.0001 decay cooler 0.00005 doubled //m/gu
    double grid_shift_magn = -1.106; //in m, compared to ions coordinates

    
     double x,z,B;
     ifstream magneticfield;
     magneticfield.open(_B_filename); 
     if (!magneticfield)
     {cout<<"Error in opening magnet file. The program will be stopped."<<endl; exit(1);}
     
     for (int i = 0;i<z_scale_mags;i++){
         for (int j = 0;j<x_scale_mags;j++){
              magneticfield>>B;
              z=(1106-i)*0.25;//veelkans nog schalen en verschuiven
              x=20.0/50.0*j;//20.0/50.0=0.04 dus ongeveer 0.25^-1
              testie<<x<<"\t"<<z<<"\t"<<B<<endl;
         }         
     }
    testie.close(); 
    magneticfield.close();
    
    cout<<"Finished reading MagneticField"<<endl;   
     
}
   

     
void getNearestNeighbours(const double& _x,const double& _y,const double& _z){
     double x,y,z,r;
     x = _x;y = _y;z= _z;
        //check of ze wel in of bound zitten he!
     if(x>r_max||y>r_max||z>z_max||x<r_min||y<r_min||z<z_min){
           cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;}

      int rs[4],zs[4];
      double r_scale=1.0/r_step;
      double z_scale=1.0/z_step;
      

      r = sqrt(x*x + y*y);
      if(r == 0)r=1e-10;

//scale both particles he.
//uppercorner is(z_min,r_min)
      r = r+(-r_min);
      z= z+(-z_min);
      
           rs[0]=(int)(r/r_step);
           rs[1]=(int)(r/r_step-0.5); //the next one he: (r+0.5*r_step)/r_step
           rs[2]=(int)(r/r_step+0.5);
           zs[0]=(int)(z/z_step);
           zs[1]=(int)(z/z_step-0.5);
           zs[2]=(int)(z/z_step+0.5);
      cout<<"************"<<endl;
      //cout<<Vmatrix(rs[0],z);    
      cout<<"Voltage@ "<<sqrt(_x*_x+_y*_y)<<","<<_z<<" is: "<<Vmatrix(zs[0],rs[0])<<"\n";
    //  cout<< Vmatrix(zs[1],rs[1])<<" "<<Vmatrix(zs[1],rs[2])<<endl;
    //  cout<< Vmatrix(zs[2],rs[1])<<" "<<Vmatrix(zs[2],rs[2])<<endl;   
}     
  */   
