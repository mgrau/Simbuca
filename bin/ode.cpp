//ode.cpp  



#include "MPI_simbuca.h"
#include "ode.h"
#include "force.h"
#include <stdio.h>
#include <string.h>



//runga kutta waarden
/*
static const double DPc2=1.0/5.0,DPc3=3.0/10.0,DPc4=4.0/5.0,DPc5=8.0/9.0,DPc6=1.0,DPc7=1.0,
DPa21=1.0/5.0,
DPa31=3.0/40.0,DPa32=9.0/40.0,
DPa41=44.0/45.0,DPa42=-56.0/15.0,DPa43=32.0/9.0,
DPa51=19372.0/6561.0,DPa52=-25360.0/2187.0,DPa53=64448.0/6561.0,DPa54=-212.0/729.0,
DPa61=9017.0/3168.0,DPa62=-355.0/33.0,DPa63=46732.0/5247.0,DPa64=49.0/176.0,DPa65=-5103.0/18656.0,
DPa71=35.0/384.0,DPa73=500.0/1113.0,DPa74=125.0/192.0,DPa75=-2187.0/6784.0, DPa76=11.0/84.0,
DPb1=35.0/384.0 ,DPb2 =0.0 ,DPb3=500.0/1113.0,DPb4=125.0/192.0 ,DPb5=-2187.0/6784.0 ,DPb6= 11.0/84.0, DPb7=0.0,
//e=b-b*
DPe1=71.0/57600.0,DPe3=-71.0/16695.0,DPe4=71.0/1920.0,
DPe5=-17253.0/339200.0,DPe6=22.0/525.0,DPe7=-1.0/40.0;
*/


//DPe1=35.0/384.0-5179.0/57600.0,DPe3=500.0/1113.0-7571.0/16695.0,DPe4=125.0/192.0-393.0/640.0,
//DPe5=-2187.0/6784.0+92097.0/339200.0,DPe6=11.0/84.0-187.0/2100.0,DPe7=-1.0/40.0;



/*DPd1=-12715105075.0/11282082432.0, DPd2 = 0.0,
DPd3=87487479700.0/32700410799.0, DPd4=-10690763975.0/1880347072.0,
DPd5=701980252875.0/199316789632.0, DPd6= -1453857185.0/822651844.0,
DPd7=69997945.0/29380423.0;*/






LogFile ologger;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//////////////Dormand_Prince Numerical Integration//////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////



void InitAcceleraction(int nrParticles,IonCloud &_cloud,_ode_vars &odev){
/*    	Init forces for Gear Method
	a = (FB+FE)/m
	
	a_x = q/m * (vy*Bz - By*vz+Ex)
	a_y = q/m * (vz*Bx - Bz*vx+Ey)
	a_z = q/m * (vx*By - Bx*vy+Ez)
*/
	double q_div_mass =0;
	double vx,vy,vz=0;
	double x,y,z=0;
	vector <double> E,B;

     for(unsigned i = 0; i+1<nrParticles; i++){
       	q_div_mass = el_charge/(*_cloud.ions)[i].Getmass();
	x=_cloud.pos[i][0];vx=_cloud.vel[i][0];
	y=_cloud.pos[i][1];vy=_cloud.vel[i][1];
	z=_cloud.pos[i][2];vz=_cloud.vel[i][2];
	bool idealtrap=false;
//	cout<<" XXX change this. Take one Runga Kutta step for init!\n";
/*	
	if(!idealtrap){
		E=getElectr(x,y,z);
		B=getMagn(x,y,z);
       		derivs[i][0]=accel_new_suggestion[i][0]=q_div_mass*(vy*B[2]-vz*B[1]+E[0]);
       		derivs[i][1]=accel_new_suggestion[i][1]=q_div_mass*(vz*B[0]-vx*B[2]+E[1]);
       		derivs[i][2]=accel_new_suggestion[i][2]=q_div_mass*(vx*B[1]-vy*B[0]+E[2]);
	       	b[i][0]=q_div_mass*(derivs[i][1]*B[2]-derivs[i][2]*B[1]+E[0]);
       		b[i][1]=q_div_mass*(derivs[i][2]*B[0]-derivs[i][0]*B[2]+E[1]);
       		b[i][2]=q_div_mass*(derivs[i][0]*B[1]-derivs[i][1]*B[0]+E[2]);	
	       	c[i][0]=q_div_mass*(b[i][1]*B[2]-b[i][2]*B[1]+E[0]);
       		c[i][1]=q_div_mass*(b[i][2]*B[0]-b[i][0]*B[2]+E[1]);
       		c[i][2]=q_div_mass*(b[i][0]*B[1]-b[i][1]*B[0]+E[2]);
	       	d[i][0]=q_div_mass*(c[i][1]*B[2]-c[i][2]*B[1]+E[0]);
       		d[i][1]=q_div_mass*(c[i][2]*B[0]-c[i][0]*B[2]+E[1]);
       		d[i][2]=q_div_mass*(c[i][0]*B[1]-c[i][1]*B[0]+E[2]);
		
*/	/*	for(int j=0;j<3;j++){cout<<"derivs["<<j<<"]= "<<derivs[i][j]<<endl;}
		for(int j=0;j<3;j++){cout<<"E["<<j<<"]= "<<E[j]<<endl;}
		for(int j=0;j<3;j++){cout<<"B["<<j<<"]= "<<B[j]<<endl;}
	*/	

//	}else{
     	  	odev.forcev.derivs[i][0] = odev.accel_new_suggestion[i][0]=1.6e9;
       		odev.forcev.derivs[i][1] = odev.accel_new_suggestion[i][1]=1.1e9;
       		odev.forcev.derivs[i][2] = odev.accel_new_suggestion[i][2]=-5e8;
			
//	}
     } 

}


void InitVV(){
	VV_BField_off();
}



double Dormand_Prince_5(IonCloud &_cloud,_ode_vars &odev){
	/*	
	goes over all the particles he!!
	vtmp is the proposed output vector
	
	
	adapts the vector latest(6) and pre_latest(6) and odev.yerr(6)
	4e orde Runga Kutta methode voor y[1..n] en de afgeleiden dydx[1..n] op x. 
	De bedoeling is om de oplossing door te trekken over een interval h, 
	en de geincrementeerde waarde terug te geven als yout[1..n] . 
	De gebruiker gebruikt Force(x,y,dydx) om de afgeleide van dydx(=v)  dus a terug te geven
	
	new_suggestions are vectors off positions[0,1,2] and velocities[3,4,5]
	k`s are vectors with velocities[0,1,2] and accelerations[3,4,5] witch are later *h and counted 
	so they are vectors off positions[0,1,2] and velocities[3,4,5]! Then they are with their 
	respective weights the new vectors wich contain the information off the system.
	*/ 
    static const double DPc2=1.0/5.0,DPc3=3.0/10.0,DPc4=4.0/5.0,DPc5=8.0/9.0,DPc6=1.0,DPc7=1.0,
    DPa21=1.0/5.0,
    DPa31=3.0/40.0,DPa32=9.0/40.0,
    DPa41=44.0/45.0,DPa42=-56.0/15.0,DPa43=32.0/9.0,
    DPa51=19372.0/6561.0,DPa52=-25360.0/2187.0,DPa53=64448.0/6561.0,DPa54=-212.0/729.0,
    DPa61=9017.0/3168.0,DPa62=-355.0/33.0,DPa63=46732.0/5247.0,DPa64=49.0/176.0,DPa65=-5103.0/18656.0,
    DPa71=35.0/384.0,DPa73=500.0/1113.0,DPa74=125.0/192.0,DPa75=-2187.0/6784.0, DPa76=11.0/84.0,
    DPb1=35.0/384.0 ,DPb2 =0.0 ,DPb3=500.0/1113.0,DPb4=125.0/192.0 ,DPb5=-2187.0/6784.0 ,DPb6= 11.0/84.0, DPb7=0.0,
    //e=b-b*
    DPe1=71.0/57600.0,DPe3=-71.0/16695.0,DPe4=71.0/1920.0,
    DPe5=-17253.0/339200.0,DPe6=22.0/525.0,DPe7=-1.0/40.0;
    int nrparticles = _cloud.nrparticles;
    odev.err=odev.sk = 0.0;
    
    unsigned i;
    double h = odev.h;
//alles berekenen, pas later vermenigvuldigen met h

for(unsigned j=0; j< nrparticles; j++){
//stap 0
	//for (i=0; i<3;i++) odev.k1[i]=old_particles[j][i+3];				//Stap 1a: speed updating
	odev.k1[j][0]= _cloud.vel[j][0];
	odev.k1[j][1]= _cloud.vel[j][1];
	odev.k1[j][2]= _cloud.vel[j][2];

	
}
    // copy pos in pos2
    // copy vel in vel2
    //memcpy (dest, source, sizeof(source) );
    memcpy(_cloud.pos2,_cloud.pos,sizeof(double)*3*nrparticles);
    memcpy(_cloud.vel2,_cloud.vel,sizeof(double)*3*nrparticles);

	force(_cloud,odev); //1

for(unsigned j=0; j< nrparticles; j++){	
	odev.k1[j][3]= odev.forcev.derivs[j][0];	//Stap 1b: accel updaten
    odev.k1[j][4]= odev.forcev.derivs[j][1];
    odev.k1[j][5]= odev.forcev.derivs[j][2];
	//for (i=0; i<3;i++) new_suggestion[j][i]=old_particles[j][i]+c2*h;			//posities updaten
	_cloud.pos2[j][0]=_cloud.pos[j][0] +h*(DPa21*odev.k1[j][0]);
	_cloud.pos2[j][1]=_cloud.pos[j][1] +h*(DPa21*odev.k1[j][1]);
	_cloud.pos2[j][2]=_cloud.pos[j][2] +h*(DPa21*odev.k1[j][2]);
	//for (i=0; i<3;i++) odev.odev.k2[i]=new_suggestion[j][i+3]=old_particles[j][i+3]+h*(a21*odev.k1[i+3]);//Stap 2a: snelheden updaten
    odev.k2[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(DPa21*odev.k1[j][3]);
    odev.k2[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(DPa21*odev.k1[j][4]);
    odev.k2[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(DPa21*odev.k1[j][5]);
} 
     force(_cloud,odev);//2
for(unsigned j=0; j< nrparticles; j++){
    odev.k2[j][3]= odev.forcev.derivs[j][0];	//Stap 2b: accel updaten
    odev.k2[j][4]= odev.forcev.derivs[j][1];
    odev.k2[j][5]= odev.forcev.derivs[j][2];
	//for (i=0; i<3;i++) new_suggestion[j][i]=old_particles[j][i]+c3*h;	//posities updaten
	_cloud.pos2[j][0]=_cloud.pos[j][0]+h*(DPa31*odev.k1[j][0]+DPa32*odev.k2[j][0]);
    _cloud.pos2[j][1]=_cloud.pos[j][1]+h*(DPa31*odev.k1[j][1]+DPa32*odev.k2[j][1]);
   	_cloud.pos2[j][2]=_cloud.pos[j][2]+h*(DPa31*odev.k1[j][2]+DPa32*odev.k2[j][2]);
	//for (i=0; i<3;i++) odev.k3[i]=new_suggestion[j][i+3]=old_particles[j][i+3]+h*(a31*odev.k1[i+3]+a32*odev.k2[i+3]);//Stap 3a: snelheden updaten
	odev.k3[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(DPa31*odev.k1[j][3]+DPa32*odev.k2[j][3]);
	odev.k3[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(DPa31*odev.k1[j][4]+DPa32*odev.k2[j][4]);
   	odev.k3[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(DPa31*odev.k1[j][5]+DPa32*odev.k2[j][5]);
}   
    force(_cloud,odev);	//3
for(unsigned j=0; j< nrparticles; j++){
    odev.k3[j][3]= odev.forcev.derivs[j][0];	//Stap 1b: accel updaten
    odev.k3[j][4]= odev.forcev.derivs[j][1];
    odev.k3[j][5]= odev.forcev.derivs[j][2];

	//for (i=0; i<3;i++) new_suggestion[j][i]=old_particles[j][i]+c4*h;		//in temp steken voor later force te berekenen
    _cloud.pos2[j][0]=_cloud.pos[j][0]+h*(DPa41*odev.k1[j][0]+DPa42*odev.k2[j][0]+DPa43*odev.k3[j][0]);
    _cloud.pos2[j][1]=_cloud.pos[j][1]+h*(DPa41*odev.k1[j][1]+DPa42*odev.k2[j][1]+DPa43*odev.k3[j][1]);
    _cloud.pos2[j][2]=_cloud.pos[j][2]+h*(DPa41*odev.k1[j][2]+DPa42*odev.k2[j][2]+DPa43*odev.k3[j][2]);
	//for (i=0; i<3;i++) odev.k4[i]=new_suggestion[j][i+3]=old_particles[j][i+3]+h*(a41*odev.k1[i+3]+a42*odev.k2[i+3]+a43*odev.k3[i+3]);		//Stap 4a: posities updaten
    odev.k4[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(DPa41*odev.k1[j][3]+DPa42*odev.k2[j][3]+DPa43*odev.k3[j][3]);
    odev.k4[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(DPa41*odev.k1[j][4]+DPa42*odev.k2[j][4]+DPa43*odev.k3[j][4]);
    odev.k4[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(DPa41*odev.k1[j][5]+DPa42*odev.k2[j][5]+DPa43*odev.k3[j][5]);
}
    force(_cloud,odev);//4
for(unsigned j=0; j< nrparticles; j++){
	odev.k4[j][3]= odev.forcev.derivs[j][0];	//Stap 1b: accel updaten
    odev.k4[j][4]= odev.forcev.derivs[j][1];
    odev.k4[j][5]= odev.forcev.derivs[j][2];
    //for (i=0; i<3;i++) new_suggestion[j][i]=old_particles[j][i]+c5*h;		//in temp steken voor later force te berekenen
    _cloud.pos2[j][0]=_cloud.pos[j][0]+h*(DPa51*odev.k1[j][0]+DPa52*odev.k2[j][0]+DPa53*odev.k3[j][0]+DPa54*odev.k4[j][0]);
    _cloud.pos2[j][1]=_cloud.pos[j][1]+h*(DPa51*odev.k1[j][1]+DPa52*odev.k2[j][1]+DPa53*odev.k3[j][1]+DPa54*odev.k4[j][1]);
    _cloud.pos2[j][2]=_cloud.pos[j][2]+h*(DPa51*odev.k1[j][2]+DPa52*odev.k2[j][2]+DPa53*odev.k3[j][2]+DPa54*odev.k4[j][2]);
	//for (i=0; i<3;i++) odev.k5[i]=new_suggestion[j][i+3]=old_particles[j][i+3]+h*(a51*odev.k1[i+3]+a52*odev.k2[i+3]+a53*odev.k3[i+3]+a54*odev.k4[i+3]);		//Stap 5a: posities updaten
	odev.k5[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(DPa51*odev.k1[j][3]+DPa52*odev.k2[j][3]+DPa53*odev.k3[j][3]+DPa54*odev.k4[j][3]);
	odev.k5[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(DPa51*odev.k1[j][4]+DPa52*odev.k2[j][4]+DPa53*odev.k3[j][4]+DPa54*odev.k4[j][4]);
	odev.k5[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(DPa51*odev.k1[j][5]+DPa52*odev.k2[j][5]+DPa53*odev.k3[j][5]+DPa54*odev.k4[j][5]);
}
    force(_cloud,odev);    	//5
for(unsigned j=0; j< nrparticles; j++){
	odev.k5[j][3]= odev.forcev.derivs[j][0];	//Stap 1b: accel updaten
    odev.k5[j][4]= odev.forcev.derivs[j][1];
    odev.k5[j][5]= odev.forcev.derivs[j][2];						//Stap 5b: snelheden updaten
	//for (i=0; i<3;i++) new_suggestion[j][i]=old_particles[j][i]+h;//DPc6=1		in temp steken voor later force te berekenen
    _cloud.pos2[j][0]=_cloud.pos[j][0]+h*(DPa61*odev.k1[j][0]+DPa62*odev.k2[j][0]+DPa63*odev.k3[j][0]+DPa64*odev.k4[j][0]+DPa65*odev.k5[j][0]);
    _cloud.pos2[j][1]=_cloud.pos[j][1]+h*(DPa61*odev.k1[j][1]+DPa62*odev.k2[j][1]+DPa63*odev.k3[j][1]+DPa64*odev.k4[j][1]+DPa65*odev.k5[j][1]);
    _cloud.pos2[j][2]=_cloud.pos[j][2]+h*(DPa61*odev.k1[j][2]+DPa62*odev.k2[j][2]+DPa63*odev.k3[j][2]+DPa64*odev.k4[j][2]+DPa65*odev.k5[j][2]);
    //for (i=0; i<3;i++) odev.k6[i]=new_suggestion[j][i+3]=old_particles[j][i+3]+h*(a61*odev.k1[i+3]+a62*odev.k2[i+3]+a63*odev.k3[i+3]+a64*odev.k4[i+3]+a65*odev.k5[i+3]);		//Stap 6a: posities updaten
	odev.k6[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(DPa61*odev.k1[j][3]+DPa62*odev.k2[j][3]+DPa63*odev.k3[j][3]+DPa64*odev.k4[j][3]+DPa65*odev.k5[j][3]);
	odev.k6[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(DPa61*odev.k1[j][4]+DPa62*odev.k2[j][4]+DPa63*odev.k3[j][4]+DPa64*odev.k4[j][4]+DPa65*odev.k5[j][4]);
	odev.k6[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(DPa61*odev.k1[j][5]+DPa62*odev.k2[j][5]+DPa63*odev.k3[j][5]+DPa64*odev.k4[j][5]+DPa65*odev.k5[j][5]);
}
	force(_cloud,odev);//6
for(unsigned j=0; j< nrparticles; j++){
    odev.k6[j][3]= odev.forcev.derivs[j][0];	//Stap 1b: accel updaten
    odev.k6[j][4]= odev.forcev.derivs[j][1];
    odev.k6[j][5]= odev.forcev.derivs[j][2];
	//for (i=0; i<3;i++) new_suggestion[j][i]=old_particles[j][i]+h;//DPc7= 1 in temp steken voor later force te berekenen
	_cloud.pos2[j][0]=_cloud.pos[j][0]+h*(DPa71*odev.k1[j][0]+DPa73*odev.k3[j][0]+DPa74*odev.k4[j][0]+DPa75*odev.k5[j][0]+DPa76*odev.k6[j][0]);
	_cloud.pos2[j][1]=_cloud.pos[j][1]+h*(DPa71*odev.k1[j][1]+DPa73*odev.k3[j][1]+DPa74*odev.k4[j][1]+DPa75*odev.k5[j][1]+DPa76*odev.k6[j][1]);
	_cloud.pos2[j][2]=_cloud.pos[j][2]+h*(DPa71*odev.k1[j][2]+DPa73*odev.k3[j][2]+DPa74*odev.k4[j][2]+DPa75*odev.k5[j][2]+DPa76*odev.k6[j][2]);
    //for (i=0; i<3;i++) odev.k7[i]=new_suggestion[j][i+3]=old_particles[j][i+3]+h*(a71*odev.k1[i+3]+a73*odev.k3[i+3]+a74*odev.k4[i+3]+a75*odev.k5[i+3]+a76*odev.k6[i+3]);		//Stap 6a: posities updaten	
    odev.k7[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(DPa71*odev.k1[j][3]+DPa73*odev.k3[j][3]+DPa74*odev.k4[j][3]+DPa75*odev.k5[j][3]+DPa76*odev.k6[j][3]);
	odev.k7[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(DPa71*odev.k1[j][4]+DPa73*odev.k3[j][4]+DPa74*odev.k4[j][4]+DPa75*odev.k5[j][4]+DPa76*odev.k6[j][4]);
	odev.k7[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(DPa71*odev.k1[j][5]+DPa73*odev.k3[j][5]+DPa74*odev.k4[j][5]+DPa75*odev.k5[j][5]+DPa76*odev.k6[j][5]);
}
	force(_cloud,odev);//7
	
for(unsigned j=0; j< nrparticles; j++){
    odev.k7[j][3]= odev.forcev.derivs[j][0];	//Stap 1b: accel updaten
    odev.k7[j][4]= odev.forcev.derivs[j][1];
    odev.k7[j][5]= odev.forcev.derivs[j][2];  
    //for (i=0; i<6;i++){old_particles[j][i]=old_particles[j][i]+h*(b1*odev.k1[i]+b3*odev.k3[i]+b4*odev.k4[i]+b5*odev.k5[i]+b6*odev.k6[i]);} 
    _cloud.pos2[j][0]=_cloud.pos[j][0]+h*(DPb1*odev.k1[j][0]+DPb3*odev.k3[j][0]+DPb4*odev.k4[j][0]+DPb5*odev.k5[j][0]+DPb6*odev.k6[j][0]);
    _cloud.pos2[j][1]=_cloud.pos[j][1]+h*(DPb1*odev.k1[j][1]+DPb3*odev.k3[j][1]+DPb4*odev.k4[j][1]+DPb5*odev.k5[j][1]+DPb6*odev.k6[j][1]);
    _cloud.pos2[j][2]=_cloud.pos[j][2]+h*(DPb1*odev.k1[j][2]+DPb3*odev.k3[j][2]+DPb4*odev.k4[j][2]+DPb5*odev.k5[j][2]+DPb6*odev.k6[j][2]);
    _cloud.vel2[j][0]=_cloud.vel[j][0]+h*(DPb1*odev.k1[j][3]+DPb3*odev.k3[j][3]+DPb4*odev.k4[j][3]+DPb5*odev.k5[j][3]+DPb6*odev.k6[j][3]);    
    _cloud.vel2[j][1]=_cloud.vel[j][1]+h*(DPb1*odev.k1[j][4]+DPb3*odev.k3[j][4]+DPb4*odev.k4[j][4]+DPb5*odev.k5[j][4]+DPb6*odev.k6[j][4]);
    _cloud.vel2[j][2]=_cloud.vel[j][2]+h*(DPb1*odev.k1[j][5]+DPb3*odev.k3[j][5]+DPb4*odev.k4[j][5]+DPb5*odev.k5[j][5]+DPb6*odev.k6[j][5]);        
}


//nu nog ff den error berekenen, moet eigenlijk alleen maar met v_vectoren!!! is lastigste he!!
/*    for (i=0; i<6;i++){
         odev.yerr[i]=h*(e1*odev.k1[i]+e3*odev.k3[i]+e4*odev.k4[i]+e5*odev.k5[i]+e6*odev.k6[i]+e7*odev.k7[i]);
         odev.sk=rk_atol+rk_rtol*max(fabs(old_particles[j][i]),fabs(new_suggestion[j][i]));
         err = err+(odev.yerr[i]/odev.sk)*(odev.yerr[i]/odev.sk); //tot de 2e
     } 
*/ 
for(unsigned j=0; j< nrparticles; j++){
    for (i=0; i<6;i++){ 
        odev.yerr[i]=h*(DPe1*odev.k1[j][i]+DPe3*odev.k3[j][i]+DPe4*odev.k4[j][i]+DPe5*odev.k5[j][i]+DPe6*odev.k6[j][i]+DPe7*odev.k7[j][i]);
    }
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.pos[j][0]),fabs(_cloud.pos2[j][0]));
   odev.err = odev.err+(odev.yerr[0]/odev.sk)*(odev.yerr[0]/odev.sk); //0
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.pos[j][1]),fabs(_cloud.pos2[j][1]));
   odev.err = odev.err+(odev.yerr[1]/odev.sk)*(odev.yerr[1]/odev.sk); //1
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.pos[j][2]),fabs(_cloud.pos2[j][2]));
   odev.err = odev.err+(odev.yerr[2]/odev.sk)*(odev.yerr[2]/odev.sk); //2
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.vel[j][0]),fabs(_cloud.vel2[j][0]));
   odev.err = odev.err+(odev.yerr[3]/odev.sk)*(odev.yerr[3]/odev.sk); //3
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.vel[j][1]),fabs(_cloud.vel2[j][1]));
   odev.err = odev.err+(odev.yerr[4]/odev.sk)*(odev.yerr[4]/odev.sk); //4
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.vel[j][2]),fabs(_cloud.vel2[j][2]));
   odev.err = odev.err+(odev.yerr[5]/odev.sk)*(odev.yerr[5]/odev.sk); //5
}//end off loops over particles
    
#ifdef __MPI_ON__
    double dum_ode;
    int ode_npart_tot;
  MPI::COMM_WORLD.Allreduce(&odev.err,&dum_ode,1,MPI::DOUBLE,MPI_SUM);
  MPI::COMM_WORLD.Allreduce(&nrparticles,&ode_npart_tot,1,MPI::INT,MPI_SUM);
  odev.err = dum_ode;
  odev.err=  sqrt(odev.err/(6*ode_npart_tot)); 
#else 
  odev.err= sqrt(odev.err/(6*_cloud.nrparticles));  
#endif // __MPI_ON__
  
    return odev.err;       
}




double RungaKutta4(IonCloud &_cloud,_ode_vars &odev){
		
//	4e orde Runga Kutta methode voor y[1..n] en de afgeleiden dydx[1..n] op x. 
//	De bedoeling is om de oplossing door te trekken over een interval h, 
//	en de geincrementeerde waarde terug te geven als yout[1..n] . 
//	De gebruiker gebruikt Force(x,y,dydx) om de afgeleide van dydx terug te geven
	

//alles berekenen, pas later vermenigvuldigen met h
//elke stap is`t zelfde:	1)in new_suggestion steken om later force te berekenen
//				2)posities updaten
//				3)kracht berekenen
//				4)snelheden updaten

//stap 1
	//for (i=0; i<3;i++) odev.k1[i]=x_n[i+3];
	
    
    
    double h = odev.h;
    unsigned i,j;
    static const double	RKb21=0.2,
    RKb31=3.0/40.0,RKb32=9.0/40.0,RKb41=0.3,RKb42 = -0.9,RKb43=1.2,
    RKb51 = -11.0/54.0, RKb52=2.5,RKb53 = -70.0/27.0,RKb54=35.0/27.0,
    RKb61=1631.0/55296.0,RKb62=175.0/512.0,RKb63=575.0/13824.0,
    RKb64=44275.0/110592.0,RKb65=253.0/4096.0,RKc1=37.0/378.0,RKc3=250.0/621.0,RKc4=125.0/594.0,
    RKc6=512.0/1771.0;
    //error stuff
    static const double RKdc5 = -277.00/14336.0,RKdc1=RKc1-2825.0/27648.0,RKdc3=RKc3-18575.0/48384.0,
    RKdc4=RKc4-13525.0/55296.0,RKdc6=RKc6-0.25;

	int     nrparticles = _cloud.nrparticles;
for(j=0; j< nrparticles; j++){    				
    odev.k1[j][0]= _cloud.vel[j][0];
    odev.k1[j][1]= _cloud.vel[j][1];
    odev.k1[j][2]= _cloud.vel[j][2];

	
	

}
    memcpy(_cloud.pos2,_cloud.pos,sizeof(double)*3*nrparticles);
    memcpy(_cloud.vel2,_cloud.vel,sizeof(double)*3*nrparticles);
	force(_cloud,odev);

for(j=0; j< nrparticles; j++){
	for (i=3; i<6;i++) odev.k1[j][i]= odev.forcev.derivs[j][i-3];			
	_cloud.pos2[j][0]=_cloud.pos[j][0]+h*(RKb21*odev.k1[j][0]);
	_cloud.pos2[j][1]=_cloud.pos[j][1]+h*(RKb21*odev.k1[j][1]);
	_cloud.pos2[j][2]=_cloud.pos[j][2]+h*(RKb21*odev.k1[j][2]);

	//for (i=0; i<3;i++) odev.k2[i]=temp[i+3]=x_n[i+3]+delta_t*(b21*odev.k1[i+3]);
	odev.k2[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(RKb21*odev.k1[j][3]);
	odev.k2[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(RKb21*odev.k1[j][4]);
	odev.k2[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(RKb21*odev.k1[j][5]);
}
	force(_cloud,odev);									
for(j=0; j< nrparticles; j++){
	for (i=3; i<6;i++) odev.k2[j][i]=odev.forcev.derivs[j][i-3];

	//for (i=0; i<3;i++) temp[i]=x_n[i]+delta_t*(b31*odev.k1[i]+b32*odev.k2[i]);
    _cloud.pos2[j][0]=_cloud.pos[j][0]+h*(RKb31*odev.k1[j][0]+RKb32*odev.k2[j][0]);
    _cloud.pos2[j][1]=_cloud.pos[j][1]+h*(RKb31*odev.k1[j][1]+RKb32*odev.k2[j][1]);
    _cloud.pos2[j][2]=_cloud.pos[j][2]+h*(RKb31*odev.k1[j][2]+RKb32*odev.k2[j][2]);
	//for (i=0; i<3;i++) odev.k3[i]=temp[i+3]=x_n[i+3]+delta_t*(b31*odev.k1[i+3]+b32*odev.k2[i+3]);
    odev.k3[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(RKb31*odev.k1[j][3]+RKb32*odev.k2[j][3]);				
    odev.k3[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(RKb31*odev.k1[j][4]+RKb32*odev.k2[j][4]);
    odev.k3[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(RKb31*odev.k1[j][5]+RKb32*odev.k2[j][5]);
}
	force(_cloud,odev);									
for(j=0; j< nrparticles; j++){
	for (i=3; i<6;i++) odev.k3[j][i]=odev.forcev.derivs[j][i-3];

	//for (i=0; i<3;i++) temp[i]=x_n[i]+delta_t*(b41*odev.k1[i]+b42*odev.k2[i]+b43*odev.k3[i]);		
	_cloud.pos2[j][0]=_cloud.pos[j][0]+h*(RKb41*odev.k1[j][0]+RKb42*odev.k2[j][0]+RKb43*odev.k3[j][0]);		
	_cloud.pos2[j][1]=_cloud.pos[j][1]+h*(RKb41*odev.k1[j][1]+RKb42*odev.k2[j][1]+RKb43*odev.k3[j][1]);
	_cloud.pos2[j][2]=_cloud.pos[j][2]+h*(RKb41*odev.k1[j][2]+RKb42*odev.k2[j][2]+RKb43*odev.k3[j][2]);
	//for (i=0; i<3;i++) odev.k4[i]=temp[i+3]=x_n[i+3]+delta_t*(b41*odev.k1[i+3]+b42*odev.k2[i+3]+b43*odev.k3[i+3]);		
	odev.k4[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(RKb41*odev.k1[j][3]+RKb42*odev.k2[j][3]+RKb43*odev.k3[j][3]);		
	odev.k4[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(RKb41*odev.k1[j][4]+RKb42*odev.k2[j][4]+RKb43*odev.k3[j][4]);	
	odev.k4[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(RKb41*odev.k1[j][5]+RKb42*odev.k2[j][5]+RKb43*odev.k3[j][5]);	
}
	force(_cloud,odev);									
for(j=0; j< nrparticles; j++){
	for (i=3; i<6;i++) odev.k4[j][i]=odev.forcev.derivs[j][i-3];						

	//for (i=0; i<3;i++) temp[i]=x_n[i]+delta_t*(b51*odev.k1[i]+b52*odev.k2[i]+b53*odev.k3[i]+b54*odev.k4[i]);		
    _cloud.pos2[j][0]=_cloud.pos[j][0]+h*(RKb51*odev.k1[j][0]+RKb52*odev.k2[j][0]+RKb53*odev.k3[j][0]+RKb54*odev.k4[j][0]);	
    _cloud.pos2[j][1]=_cloud.pos[j][1]+h*(RKb51*odev.k1[j][1]+RKb52*odev.k2[j][1]+RKb53*odev.k3[j][1]+RKb54*odev.k4[j][1]);	
    _cloud.pos2[j][2]=_cloud.pos[j][2]+h*(RKb51*odev.k1[j][2]+RKb52*odev.k2[j][2]+RKb53*odev.k3[j][2]+RKb54*odev.k4[j][2]);	
//	for (i=0; i<3;i++) odev.k5[i]=temp[i+3]=x_n[i+3]+delta_t*(b51*odev.k1[i+3]+b52*odev.k2[i+3]+b53*odev.k3[i+3]+b54*odev.k4[i+3]);		
    odev.k5[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(RKb51*odev.k1[j][3]+RKb52*odev.k2[j][3]+RKb53*odev.k3[j][3]+RKb54*odev.k4[j][3]);	
    odev.k5[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(RKb51*odev.k1[j][4]+RKb52*odev.k2[j][4]+RKb53*odev.k3[j][4]+RKb54*odev.k4[j][4]);	
    odev.k5[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(RKb51*odev.k1[j][5]+RKb52*odev.k2[j][5]+RKb53*odev.k3[j][5]+RKb54*odev.k4[j][5]);	
}
	force(_cloud,odev);									
for(j=0; j< nrparticles; j++){
	for (i=3; i<6;i++) odev.k5[j][i]=odev.forcev.derivs[j][i-3];

	//for (i=0; i<3;i++) temp[i]=x_n[i]+delta_t*(b61*odev.k1[i]+b62*odev.k2[i]+b63*odev.k3[i]+b64*odev.k4[i]+b65*odev.k5[i]);		
	_cloud.pos2[j][0]=_cloud.pos[j][0]+h*(RKb61*odev.k1[j][0]+RKb62*odev.k2[j][0]+RKb63*odev.k3[j][0]+RKb64*odev.k4[j][0]+RKb65*odev.k5[j][0]);		
	_cloud.pos2[j][1]=_cloud.pos[j][1]+h*(RKb61*odev.k1[j][1]+RKb62*odev.k2[j][1]+RKb63*odev.k3[j][1]+RKb64*odev.k4[j][1]+RKb65*odev.k5[j][1]);		
	_cloud.pos2[j][2]=_cloud.pos[j][2]+h*(RKb61*odev.k1[j][2]+RKb62*odev.k2[j][2]+RKb63*odev.k3[j][2]+RKb64*odev.k4[j][2]+RKb65*odev.k5[j][2]);		
	//for (i=0; i<3;i++) odev.k6[i]=temp[i+3]=x_n[i+3]+delta_t*(b61*odev.k1[i+3]+b62*odev.k2[i+3]+b63*odev.k3[i+3]+b64*odev.k4[i+3]+b65*odev.k5[i+3]);		
    odev.k6[j][0]=_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(RKb61*odev.k1[j][3]+RKb62*odev.k2[j][3]+RKb63*odev.k3[j][3]+RKb64*odev.k4[j][3]+RKb65*odev.k5[j][3]);		
    odev.k6[j][1]=_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(RKb61*odev.k1[j][4]+RKb62*odev.k2[j][4]+RKb63*odev.k3[j][4]+RKb64*odev.k4[j][4]+RKb65*odev.k5[j][4]);		
    odev.k6[j][2]=_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(RKb61*odev.k1[j][5]+RKb62*odev.k2[j][5]+RKb63*odev.k3[j][5]+RKb64*odev.k4[j][5]+RKb65*odev.k5[j][5]);		
}
	force(_cloud,odev);
for(j=0; j< nrparticles; j++){
	for ( i=3; i<6;i++) odev.k6[j][i]=odev.forcev.derivs[j][i-3];
	_cloud.pos2[j][0]=_cloud.pos[j][0]+h*(RKc1*odev.k1[j][0]+RKc3*odev.k3[j][0]+RKc4*odev.k4[j][0]+RKc6*odev.k6[j][0]);
	_cloud.pos2[j][1]=_cloud.pos[j][1]+h*(RKc1*odev.k1[j][1]+RKc3*odev.k3[j][1]+RKc4*odev.k4[j][1]+RKc6*odev.k6[j][1]);
	_cloud.pos2[j][2]=_cloud.pos[j][2]+h*(RKc1*odev.k1[j][2]+RKc3*odev.k3[j][2]+RKc4*odev.k4[j][2]+RKc6*odev.k6[j][2]);
	_cloud.vel2[j][0]=_cloud.vel[j][0]+h*(RKc1*odev.k1[j][3]+RKc3*odev.k3[j][3]+RKc4*odev.k4[j][3]+RKc6*odev.k6[j][3]);
	_cloud.vel2[j][1]=_cloud.vel[j][1]+h*(RKc1*odev.k1[j][4]+RKc3*odev.k3[j][4]+RKc4*odev.k4[j][4]+RKc6*odev.k6[j][4]);
	_cloud.vel2[j][2]=_cloud.vel[j][2]+h*(RKc1*odev.k1[j][5]+RKc3*odev.k3[j][5]+RKc4*odev.k4[j][5]+RKc6*odev.k6[j][5]);
    odev.yerr[0] = h*(RKdc1*odev.k1[j][0]+RKdc3*odev.k3[j][0]+RKdc4*odev.k4[j][0]+RKdc5*odev.k5[j][0]+RKdc6*odev.k6[j][0]);
    odev.yerr[1] = h*(RKdc1*odev.k1[j][1]+RKdc3*odev.k3[j][1]+RKdc4*odev.k4[j][1]+RKdc5*odev.k5[j][1]+RKdc6*odev.k6[j][1]);
    odev.yerr[2] = h*(RKdc1*odev.k1[j][2]+RKdc3*odev.k3[j][2]+RKdc4*odev.k4[j][2]+RKdc5*odev.k5[j][2]+RKdc6*odev.k6[j][2]);
    odev.yerr[3] = h*(RKdc1*odev.k1[j][3]+RKdc3*odev.k3[j][3]+RKdc4*odev.k4[j][3]+RKdc5*odev.k5[j][3]+RKdc6*odev.k6[j][3]);
    odev.yerr[4] = h*(RKdc1*odev.k1[j][4]+RKdc3*odev.k3[j][4]+RKdc4*odev.k4[j][4]+RKdc5*odev.k5[j][4]+RKdc6*odev.k6[j][4]);
    odev.yerr[5] = h*(RKdc1*odev.k1[j][5]+RKdc3*odev.k3[j][5]+RKdc4*odev.k4[j][5]+RKdc5*odev.k5[j][5]+RKdc6*odev.k6[j][5]);   
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.pos[j][0]),fabs(_cloud.pos2[j][0]));
   odev.err = odev.err+(odev.yerr[0]/odev.sk)*(odev.yerr[0]/odev.sk); //0
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.pos[j][1]),fabs(_cloud.pos2[j][1]));
   odev.err = odev.err+(odev.yerr[1]/odev.sk)*(odev.yerr[1]/odev.sk); //1
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.pos[j][2]),fabs(_cloud.pos2[j][2]));
   odev.err = odev.err+(odev.yerr[2]/odev.sk)*(odev.yerr[2]/odev.sk); //2
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.vel[j][0]),fabs(_cloud.vel2[j][0]));
   odev.err = odev.err+(odev.yerr[3]/odev.sk)*(odev.yerr[3]/odev.sk); //3
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.vel[j][1]),fabs(_cloud.vel2[j][1]));
   odev.err = odev.err+(odev.yerr[4]/odev.sk)*(odev.yerr[4]/odev.sk); //4
   odev.sk=odev.rk_atol+odev.rk_rtol*max(fabs(_cloud.vel[j][2]),fabs(_cloud.vel2[j][2]));
   odev.err = odev.err+(odev.yerr[5]/odev.sk)*(odev.yerr[5]/odev.sk); //5
}//end off loops over particles

#ifdef __MPI_ON__
    double dum_ode;
    int ode_npart_tot;
  MPI::COMM_WORLD.Allreduce(&odev.err,&dum_ode,1,MPI::DOUBLE,MPI_SUM);
  MPI::COMM_WORLD.Allreduce(&nrparticles,&ode_npart_tot,1,MPI::INT,MPI_SUM);
  odev.err = dum_ode;
  odev.err=  sqrt(odev.err/(6*ode_npart_tot)); 
#else 
  odev.err= sqrt(odev.err/(6*_cloud.nrparticles));  
#endif // __MPI_ON__
   
    return odev.err; 

}


void GearMethod(IonCloud &_cloud,_ode_vars &odev){
   int nrparticles = _cloud.nrparticles;
    double h = odev.h;
   //note: To init everything, you can also take a Runga Kutta step... Smart he!
    
    double GEa0 = 3./20.;      //(f\FCr 7. Ordnung: 1925./14112.;)        //(f\FCr 6. Ordnung: 863./6048.;)       //(f\FCr 5. Ordnung: 3./20.;)
    double GEa1 = 251./360.;     //(f\FCr 7. Ordnung: 19087./30240.;)       //(f\FCr 6. Ordnung:665./1008.;)        //(f\FCr 5. Ordnung: 251./360.;)
    double GEa2 = 1./1.;
    double GEa3 = 11./18.;        //(f\FCr 7. Ordnung: 137./180.;)           //(f\FCr 6. Ordnung: 25./36.;)          //(f\FCr 5. Ordnung: 11./18.;)
    double GEa4 = 1./6.;        //(f\FCr 7. Ordnung: 5./16.;)              //(f\FCr 6. Ordnung: 35./144.;)         //(f\FCr 5. Ordnung: 1./6.;)
    double GEa5 = 1./60.;         //(f\FCr 7. Ordnung: 17./240.;)            //(f\FCr 6. Ordnung: 1./24.;)           //(f\FCr 5. Ordnung: 1./60.;)
    double GEa6 = 1./360.;         //(f\FCr 7. Ordnung: 1./120.;)             //(f\FCr 6. Ordnung: 1./360.;)
    //double a7 = 1./2520.;
    
    double GearError5=-10./137.;
	if(!odev.Gear_initialized){
		InitAcceleraction(nrparticles,_cloud,odev);odev.Gear_initialized=true;
		//cout<<"Gear arrays initialized. cin\n";
	}
		
	
 
    
//    cout<<"Check if vectors are initialized! XXX .\n";

  for(int i=0;i<nrparticles;i++) {
	//printf("t=%e:\n",t);
	_cloud.pos2[i][0] = _cloud.pos[i][0] + _cloud.vel[i][0]*h + odev.forcev.derivs[i][0]*odev.h2f + odev.b[i][0]*odev.h3f + odev.c[i][0]*odev.h4f + odev.d[i][0]*odev.h5f;// + e[i][0]*odev.h6f;// + f[i].x*h7f;
 	_cloud.pos2[i][1] = _cloud.pos[i][1] + _cloud.vel[i][1]*h + odev.forcev.derivs[i][1]*odev.h2f + odev.b[i][1]*odev.h3f + odev.c[i][1]*odev.h4f + odev.d[i][1]*odev.h5f;// + e[i][1]*odev.h6f;// + f[i].y*h7f;
 	_cloud.pos2[i][2] = _cloud.pos[i][2] + _cloud.vel[i][2]*h + odev.forcev.derivs[i][2]*odev.h2f + odev.b[i][2]*odev.h3f + odev.c[i][2]*odev.h4f + odev.d[i][2]*odev.h5f;// + e[i][2]*odev.h6f;// + f[i].z*h7f;
 	//printf("rn=(%e,%e,%e)\n",rn.x,rn.y,rn.z);
 	_cloud.vel2[i][0] = _cloud.vel[i][0] + odev.forcev.derivs[i][0]*h + odev.b[i][0]*odev.h2f + odev.c[i][0]*odev.h3f + odev.d[i][0]*odev.h4f;// + e[i][0]*h5f;// + f[i].x*odev.h6f;
 	_cloud.vel2[i][1] = _cloud.vel[i][1] + odev.forcev.derivs[i][1]*h + odev.b[i][1]*odev.h2f + odev.c[i][1]*odev.h3f + odev.d[i][1]*odev.h4f;// + e[i][1]*h5f;// + f[i].y*odev.h6f;
 	_cloud.vel2[i][2] = _cloud.vel[i][2] + odev.forcev.derivs[i][2]*h + odev.b[i][2]*odev.h2f + odev.c[i][2]*odev.h3f + odev.d[i][2]*odev.h4f;// + e[i][2]*h5f;// + f[i].z*odev.h6f;
 	//printf("vn=(%e,%e,%e)\n",vn.x,vn.y,vn.z);
 	odev.accel_new_suggestion[i][0] = odev.forcev.derivs[i][0] + odev.b[i][0]*h + odev.c[i][0]*odev.h2f + odev.d[i][0]*odev.h3f;// + e[i][0]*h4f;// + f[i].x*h5f;
 	odev.accel_new_suggestion[i][1] = odev.forcev.derivs[i][1] + odev.b[i][1]*h + odev.c[i][1]*odev.h2f + odev.d[i][1]*odev.h3f;// + e[i][1]*h4f;// + f[i].y*h5f;
 	odev.accel_new_suggestion[i][2] = odev.forcev.derivs[i][2] + odev.b[i][2]*h + odev.c[i][2]*odev.h2f + odev.d[i][2]*odev.h3f;// + e[i][2]*h4f;// + f[i].z*h5f;
 	//printf("an=(%e,%e,%e)\n",an.x,an.y,an.z);
 	odev.bn[i][0] = odev.b[i][0] + odev.c[i][0]*h + odev.d[i][0]*odev.h2f;// + e[i][0]*h3f;// + f[i].x*h4f;
 	odev.bn[i][1] = odev.b[i][1] + odev.c[i][1]*h + odev.d[i][1]*odev.h2f;// + e[i][1]*h3f;// + f[i].y*h4f;
 	odev.bn[i][2] = odev.b[i][2] + odev.c[i][2]*h + odev.d[i][2]*odev.h2f;// + e[i][2]*h3f;// + f[i].z*h4f;
 	//printf("bn=(%e,%e,%e)\n",bn.x,bn.y,bn.z);
 	odev.cn[i][0] = odev.c[i][0] + odev.d[i][0]*h;// + e[i][0]*odev.h2f;// + f[i].x*h3f;
 	odev.cn[i][1] = odev.c[i][1] + odev.d[i][1]*h;// + e[i][1]*odev.h2f;// + f[i].y*h3f;
 	odev.cn[i][2] = odev.c[i][2] + odev.d[i][2]*h;// + e[i][2]*odev.h2f;// + f[i].z*h3f;
 	//printf("cn=(%e,%e,%e)\n",cn.x,cn.y,cn.z);
 	odev.dn[i][0] = odev.d[i][0];// + e[i][0]*h;// + f[i].x*odev.h2f;
 	odev.dn[i][1] = odev.d[i][1];// + e[i][1]*h;// + f[i].y*odev.h2f;
 	odev.dn[i][2] = odev.d[i][2];// + e[i][2]*h;// + f[i].z*odev.h2f;

	//for error calc atm.
	odev.en[i][0]=0.0;
	odev.en[i][1]=0.0;
	odev.en[i][2]=0.0;
	
	//en[i].x = e[i].x;// + f[i].x*h;
	//en[i].y = e[i].y;// + f[i].y*h;
	//en[i].z = e[i].z;// + f[i].z*h;

	//fn[i].x = f[i].x;
	//fn[i].y = f[i].y;
	//fn[i].z = f[i].z;
 	//printf("dn=(%e,%e,%e)\n",dn.x,dn.y,dn.z);
/*	place[i].x=rn[i].x;
	place[i].y=rn[i].y;
	place[i].z=rn[i].z;
	velo[i].x=vn[i].x;
	velo[i].y=vn[i].y;
	velo[i].z=vn[i].z;
	printf("v_x=%lf\n", velo[i].x);
	printf("v_y=%lf\n", velo[i].y);
	printf("v_z=%lf\n", velo[i].z);*/
    }

// gesamtbeschleunigung(t);
force(_cloud,odev);

   for(int i=0;i<nrparticles;i++) {

 	//printf("a_kraft=(%e,%e,%e)\n",a_kraft.x,a_kraft.y,a_kraft.z);
 	odev.a_diff[i][0] = odev.h2f*(odev.forcev.derivs[i][0]-odev.accel_new_suggestion[i][0]);//h*h/2.*
 	odev.a_diff[i][1] = odev.h2f*(odev.forcev.derivs[i][1]-odev.accel_new_suggestion[i][1]);//h*h/2.*
 	odev.a_diff[i][2] = odev.h2f*(odev.forcev.derivs[i][2]-odev.accel_new_suggestion[i][2]);//h*h/2.*
 	//printf("diff=(%e,%e,%e)\n",diff.x,diff.y,diff.z);
 	_cloud.pos2[i][0] +=  GEa0*odev.a_diff[i][0];
 	_cloud.pos2[i][1] +=  GEa0*odev.a_diff[i][1];
 	_cloud.pos2[i][2] +=  GEa0*odev.a_diff[i][2];
  	_cloud.vel2[i][0] +=  GEa1*odev.a_diff[i][0]/h;
 	_cloud.vel2[i][1] +=  GEa1*odev.a_diff[i][1]/h;
 	_cloud.vel2[i][2] +=  GEa1*odev.a_diff[i][2]/h;

 	//acc[i].x = an[i].x + a2*diff[i].x/odev.h2f;
 	//acc[i].y = an[i].y + a2*diff[i].y/odev.h2f;
 	//acc[i].z = an[i].z + a2*diff[i].z/odev.h2f;
 	odev.b[i][0] = odev.bn[i][0] + GEa3*odev.a_diff[i][0]/odev.h3f;
 	odev.b[i][1] = odev.bn[i][1] + GEa3*odev.a_diff[i][1]/odev.h3f;
 	odev.b[i][2] = odev.bn[i][2] + GEa3*odev.a_diff[i][2]/odev.h3f;
 	odev.c[i][0] = odev.cn[i][0] + GEa4*odev.a_diff[i][0]/odev.h4f;
 	odev.c[i][1] = odev.cn[i][1] + GEa4*odev.a_diff[i][1]/odev.h4f;
 	odev.c[i][2] = odev.cn[i][2] + GEa4*odev.a_diff[i][2]/odev.h4f;
 	odev.d[i][0] = odev.dn[i][0] + GEa5*odev.a_diff[i][0]/odev.h5f;
 	odev.d[i][1] = odev.dn[i][1] + GEa5*odev.a_diff[i][1]/odev.h5f;
 	odev.d[i][2] = odev.dn[i][2] + GEa5*odev.a_diff[i][2]/odev.h5f;
	odev.e[i][0] = odev.en[i][0] + GEa6*odev.a_diff[i][0]/odev.h6f;
	odev.e[i][1] = odev.en[i][1] + GEa6*odev.a_diff[i][1]/odev.h6f;
	odev.e[i][2] = odev.en[i][2] + GEa6*odev.a_diff[i][2]/odev.h6f;
	//f[i].x = fn[i].x + a7*odev.a_diff[i][0]/h7f;
	//f[i].y = fn[i].y + a7*odev.a_diff[i][1]/h7f;
	//f[i].z = fn[i].z + a7*odev.a_diff[i][2]/h7f;
   }//end off loops over particles	
   // GearError=-10.0/137*pow(h,5)*e[i][x]
   
   
   
   /*
   for(int i=0;i<nrparticles;i++) {
   	odev.yerr[0]=GearError5*pow(h,5)*e[i][0];
    	odev.yerr[1]=GearError5*pow(h,5)*e[i][1];
    	odev.yerr[2]=GearError5*pow(h,5)*e[i][2];
    	
	odev.sk=rk_atol+rk_rtol*max(fabs(_cloud.pos[j][0]),fabs(_cloud.pos2[j][0]));
    	err = err+(odev.yerr[0]/odev.sk)*(odev.yerr[0]/odev.sk); //0
    	odev.sk=rk_atol+rk_rtol*max(fabs(_cloud.pos[j][1]),fabs(_cloud.pos2[j][1]));
    	err = err+(odev.yerr[1]/odev.sk)*(odev.yerr[1]/odev.sk); //1
    	odev.sk=rk_atol+rk_rtol*max(fabs(_cloud.pos[j][2]),fabs(_cloud.pos2[j][2]));
    	err = err+(odev.yerr[2]/odev.sk)*(odev.yerr[2]/odev.sk); //2
   }
        err= sqrt(err/(6*_cloud.nrparticles));       
       	cout<<"error: "<<err<<endl;      
   	cin.get();	
	
    */
  // return err; 	
}

void VelocityVerlet(IonCloud &_cloud,_ode_vars &odev){
    	// Velocity Verlet algo to approximate the B field since wc=qB/m is the fastest motion to resolve
	// Based on Boris approach: 
	// 1. Add half an electric impulse to the positions
	// 2. rotate the velocity vector
	// 3. Add the second half of the electric impulse to the positions
	// Implementation done as in the PhD thesis of G. Maero.
	// limitations: assumes a uniform B field.

    int nrparticles = _cloud.nrparticles;
    double h = odev.h;
    if(!odev.VV_initialized){
	InitVV();odev.VV_initialized=true;
	InitAcceleraction(nrparticles,_cloud,odev);
    }

    double DeltaT = odev.h; // 1 ns
    double ax,ay,az = 0.;
    double axo,ayo,azo = 0; //the old accelerations.
    long double omega, tanTheta05, sinTheta, cosTheta, S,C = 0.;
    double vxtmp;
    for (unsigned int j_=0; j_< nrparticles; j_++){ 

    	omega= _cloud.wc[j_];
    	tanTheta05=tan(0.5*omega*DeltaT);
    	sinTheta=2*tanTheta05/(1+tanTheta05*tanTheta05);
    	cosTheta=(1-tanTheta05*tanTheta05)/(1+tanTheta05*tanTheta05);
    	S=sinTheta-(omega*DeltaT);
    	C=cosTheta-1;

	//assign the old accelerations
        axo = odev.forcev.derivs[j_][0];
        ayo = odev.forcev.derivs[j_][1];
        azo = odev.forcev.derivs[j_][2];

        // First get electric field end "get fields", advance POSITION now. 
	_cloud.pos2[j_][0] = _cloud.pos[j_][0]+1/omega*( _cloud.vel[j_][0]*sinTheta-_cloud.vel[j_][1]*C)+1/omega/omega*(-ax*C-ay*S);         
        _cloud.pos2[j_][1] = _cloud.pos[j_][1]-1/omega*(-_cloud.vel[j_][1]*sinTheta-_cloud.vel[j_][0]*C)+1/omega/omega*(-ay*C+ax*S);          
        _cloud.pos2[j_][2] = _cloud.pos[j_][2]+DeltaT*_cloud.vel[j_][2]+0.5*DeltaT*DeltaT*az;

       // general step, calculate velocity, advance position
	force(_cloud,odev);
	ax = odev.forcev.derivs[j_][0];
	ay = odev.forcev.derivs[j_][1];
	az = odev.forcev.derivs[j_][2];  
        
         // end "get fields", calculate VELOCITY now at time T
        vxtmp = _cloud.vel[j_][0]*cosTheta+_cloud.vel[j_][1]*sinTheta+1/omega*(-ayo*C+axo*sinTheta)+
		+1/omega/omega*(-(ax-axo)/DeltaT*C-(ay-ayo)/DeltaT*S);
        _cloud.vel2[j_][1] = _cloud.vel[j_][1]*cosTheta-_cloud.vel[j_][0]*sinTheta-1/omega*(-axo*C-ayo*sinTheta)+
		+1/omega/omega*(-(ay-ayo)/DeltaT*C+(ax-axo)/DeltaT*S);

        _cloud.vel2[j_][0] = vxtmp;
        _cloud.vel2[j_][2] = _cloud.vel[j_][2]+0.5*DeltaT*(az+azo);

        // advance POSITION to T+DeltaT
	// Now the second part of leap frog changing the position   
        _cloud.pos2[j_][0] = _cloud.pos[j_][0]+1/omega*(_cloud.vel[j_][0]*sinTheta-_cloud.vel[j_][1]*C)+1/omega/omega*(-axo*C-ayo*S); 
        _cloud.pos2[j_][1] = _cloud.pos[j_][1]-1/omega*(-_cloud.vel[j_][1]*sinTheta-_cloud.vel[j_][0]*C)+1/omega/omega*(-ayo*C+axo*S);          
        _cloud.pos2[j_][2] = _cloud.pos[j_][2]+DeltaT*_cloud.vel[j_][2]+0.5*DeltaT*DeltaT*azo;

    }  // loop over ions, j_s
  
}

bool error_succes(const double err,_ode_vars &odev){
   //Finally, the controller tests whether err 1 and adjusts the stepsize. The
   //default setting is beta = 0 (no PI control). Set beta to 0.04 or 0.08 to turn on PI control
   //Returns true if err <= 1, false otherwise. If step was successful, sets hnext to the estimated
   //optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.
   static const double beta=0.08,alpha=0.2-beta*0.75,safe=0.9, //beta=0.4/k en k=5 hier want 5e orde
   minscale=0.2,maxscale=10.0;
   //Set beta to a nonzero value for PI control. beta D 0:04\960.08 is a good default. !!!
   double scale;
   if (err <= 1.0) { //Step succeeded. Compute hnext.
      if (err == 0.0)
         scale=maxscale;
      else { //PI control if beta \A4 0.
          scale=safe*pow(err,-alpha)*pow(odev.errold,beta);
          if (scale<minscale) scale=minscale; //Ensure minscale <= hnext/h <= maxscale.
          if (scale>maxscale) scale=maxscale;
      }
      if (odev.reject) //Dont let step increase if last one was odev.rejected
          odev.hnext=odev.h*min(scale,1.0); //odev.rejected.
      else
          odev.hnext=odev.h*scale;
      odev.errold=max(err,1.0e-4); //Bookkeeping for next call.
      odev.reject=false;
      return true;
   } 
   else { //Truncation error too large, reduce stepsize.
       scale=max(safe*pow(err,-alpha),minscale);
       odev.h *= scale;
       odev.reject=true;
       return false;
   }
}      
     


void step(IonCloud &_cloud,_ode_vars &odev){
    /* 
     Attempts a step with stepsize htry. On output, positions (x) and velocities (y) are replaced by their new values, hdid
     is the stepsize that was actually accomplished, and hnext is the estimated next stepsize.
     Doub h=htry; Set stepsize to the initial trial value.
	 */
    
    int nrparticles = _cloud.nrparticles;
	
	double error;
	//if(!poolvectorsInitialized){Initpoolvectors(_cloud.nrparticles);poolvectorsInitialized = true;}
	
    double particletime = _cloud.lifetime - odev.time_ini_ope;
    if(odev.sweep_flag)
    {
        odev.w_exc = odev.sweep_wi + (odev.sweep_wf-odev.sweep_wi)*(_cloud.lifetime - odev.initial_time)/(odev.final_time-odev.initial_time);
        if(odev.sweep_par!=0)
        {
            odev.w_exc = odev.sweep_wi + odev.sweep_par*(_cloud.lifetime - odev.initial_time) + odev.sweep_wf*pow(_cloud.lifetime - odev.initial_time,2);
            
        }
        
        //odev.U_exc = odev.sweep_wi + (odev.sweep_wf-odev.sweep_wi)*(_cloud.lifetime - odev.initial_time)/(odev.final_time-odev.initial_time);
    }
    
    if(odev.swift_flag)
    {
        if(odev.swift_RW)
        {
            odev.forcev.coswTtimeTU = odev.Interpolate_SWIFT(particletime,0);
            odev.forcev.sinwTtimeTU = odev.Interpolate_SWIFT(particletime,1);
        }
        else
        {
            odev.forcev.coswTtimeTU = odev.Interpolate_SWIFT(particletime);
        }
    }
    else
    {
        odev.forcev.coswTtimeTU=cos(odev.w_exc*particletime)*odev.U_exc;
        odev.forcev.sinwTtimeTU=sin(odev.w_exc*particletime)*odev.U_exc;
        odev.forcev.coswTtimeTU2=cos(odev.w_exc2*particletime)*odev.U_exc2;
        odev.forcev.coswTtimeTU3=cos(odev.w_exc3*particletime)*odev.U_exc3;
        odev.forcev.coswTtimeTU4=cos(odev.w_exc4*particletime)*odev.U_exc4;
        odev.forcev.cos2wTtimeTU=cos(2.*odev.w_exc*particletime)*odev.U_exc;
        odev.forcev.sin2wTtimeTU=sin(2.*odev.w_exc*particletime)*odev.U_exc;
    }
    
    
    for (;;) {   
		
        if(odev.RK4 == false && odev.DP5 == true && odev.VV == false) error = Dormand_Prince_5(_cloud,odev);
        if(odev.RK4 == true  && odev.DP5 == false && odev.VV == false) error = RungaKutta4(_cloud,odev);
        if(odev.RK4 == false && odev.DP5 == false && odev.VV == true){
		VelocityVerlet(_cloud,odev);
		odev.change_stepsize = false;
	}
        if(odev.RK4 == false && odev.DP5 == false && odev.VV == false) { // in the case of the Gear Method!
			GearMethod(_cloud,odev);
			//VerletMethod();
			odev.change_stepsize = false;
        }
        
        if(odev.change_stepsize == true){
/*	  
#ifdef __MPI_ON__
			MPI::COMM_WORLD.Allreduce(&error,&dum_ode,1,MPI::DOUBLE,MPI_MAX);
			error = dum_ode;
#endif // __MPI_ON__
*/
			odev.countstep++;
			if (error_succes(error,odev))
				break; //Step odev.rejected. Try again with reduced h set		 
			if (fabs(odev.h) <= fabs(_cloud.pos[0][0])*odev.EPS) //It was...fabs(h) <= fabs(particles[0].px)*EPS XXX this [0][0] should be different for all
				throw("stepsize underflow created");
        }else{
			break;//don`t change stepsize, go out off loop.
			//hnext stays the same, as initialized in InitOde.
        }
		
    }
	// copy pos2 in pos
    // copy vel2 in vel
    memcpy(_cloud.pos,_cloud.pos2,sizeof(double)*3*nrparticles);
    memcpy(_cloud.vel,_cloud.vel2,sizeof(double)*3*nrparticles);
    
    
    
    //update particles time he.
    _cloud.lifetime = _cloud.lifetime+odev.h;        
    odev.h=odev.hnext; //hnext is calculated in error stuff
   
    //printf("step =%e\n",h);
}


double GetTimeStep(_ode_vars &odev){
       return odev.h;
}

//int GetCountStep(){
//     return countstep;
//}

_ode_vars::_ode_vars()
{
     //init the vectors for 1 particle:
      k1 = new double[1][6];
      k2 = new double[1][6];
      k3 = new double[1][6];
      k4 = new double[1][6];
      k5 = new double[1][6];
      k6 = new double[1][6];

     httry=1e-9; //maybe -9 is beter? Sam also used 1e-9
   
     RK4 = true;
     DP5 = false;
     change_stepsize = true;
     Gear_initialized = false;
     VV_initialized=false;
     rk_atol=1e-7;
     rk_rtol=1e-7;

     Gdt= 1e-9;
       
    EPS=numeric_limits<double>::epsilon();
    poolvectorsInitialized = false;      //standard false he
    withCharge = false;
    PrintAfterOperation = false;
    PrintatZpos_bool =false;
    PrintZpos = 0;
    yerr.resize(6);
    // force
    
    
    // sweep
    sweep_flag = false;
    swift_flag = false;
    swift_RW = false;


}

_ode_vars::~_ode_vars()
{
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
    delete [] k6;
    delete [] k7;

}

void _ode_vars::Init_ode(int _ode_order, double _timestep, bool _adaptive_stepsize)
{
    errold=1.0e-4;
    reject = false;
    httry = _timestep;
    //XXX hopefully odev.httry = _timestep still OK for runga kutta!
    h=httry;
    countstep = 0;
    hnext=_timestep;reset_timestep= _timestep;
    
    h2f=h*h/2.;
    h3f=pow(h,3)/6.;
    h4f=pow(h,4)/24.;
    h5f=pow(h,5)/120.;
    h6f=pow(h,6)/720.;
    
    change_stepsize = _adaptive_stepsize;
    
    switch (_ode_order) {
        case 0: {
            ologger<<"Choose Velocity Verlet integration method.";
            RK4 = false;
            DP5 = false;
	    VV = true;
            h=_timestep;
            hnext=_timestep;
            break;
        }
        case 1: {
            ologger<<"Choose Gear Method 5th order. Leapfrog";
            RK4 = false;
            DP5 = false;
	    VV = false;
            h=_timestep;
            hnext=_timestep;
            break;
        }
        case 4: {
            ologger<<"Choose Runga Kutta 4th order.";
            RK4 = true;
            DP5 = false;
	    VV = false;
            break;
        }
        case 5: {
            ologger<<"Choose Dormand Prince 5th order.";
            DP5 = true;
            RK4 = false;
            VV = false;
            break;
        }
        default: {
            ologger<<"Choose a possible order to perform ode's 0,1,4 or 5.\n";
            ologger<<"By default Runge Kutta 4th order is choosen.";
            RK4 = true;
            DP5 = false;
        }
    }
    if(_adaptive_stepsize){ologger<<" With adaptive stepsize ";}else{ologger<<" With fixed stepsize ";}
    ologger<<"and abs error: ";ologger<<rk_atol;ologger<<", rel error: ";ologger<<rk_rtol;ologger<<"\n";
    poolvectorsInitialized = false;
}

void _ode_vars::Reset_ode(){
    h=httry;
    hnext=reset_timestep;
    h2f=h*h/2.;
    h3f=pow(h,3)/6.;
    h4f=pow(h,4)/24.;
    h5f=pow(h,5)/120.;
    h6f=pow(h,6)/720.;
    
    
}

void _ode_vars::Initpoolvectors(int nrParticles){
    //Init kvectors and the returning force vectors. This has to be done because
    //nrParticles is a dynamic number!!
//    k1.resize(nrParticles);
//    k2.resize(nrParticles);
//    k3.resize(nrParticles);
//    k4.resize(nrParticles);
//    k5.resize(nrParticles);
//    k6.resize(nrParticles);
//    k7.resize(nrParticles);
    
    k1 = new double[nrParticles][6];
    k2 = new double[nrParticles][6];
    k3 = new double[nrParticles][6];
    k4 = new double[nrParticles][6];
    k5 = new double[nrParticles][6];
    k6 = new double[nrParticles][6];
    k7 = new double[nrParticles][6];

    
    bn.resize(nrParticles);
    cn.resize(nrParticles);
    dn.resize(nrParticles);
    en.resize(nrParticles);
    fn.resize(nrParticles);
    a_diff.resize(nrParticles);
    
    b.resize(nrParticles);
    c.resize(nrParticles);
    d.resize(nrParticles);
    e.resize(nrParticles);
    
    //forcev.derivs.resize(nrParticles);
    forcev.derivs = new double[nrParticles][3];
    
    accel_new_suggestion.resize(nrParticles);
    Ion tmpion;
    
    for (unsigned int i = 0; i < nrParticles; ++i){
//        k1[i].resize(6);
//        k2[i].resize(6);
//        k3[i].resize(6);
//        k4[i].resize(6);
//        k5[i].resize(6);
//        k6[i].resize(6);
//        k7[i].resize(6);
//        forcev.derivs[i].resize(3);
        accel_new_suggestion[i].resize(3);
        a_diff[i].resize(3);
        bn[i].resize(3);
        cn[i].resize(3);
        dn[i].resize(3);
        b[i].resize(3);
        c[i].resize(3);
        d[i].resize(3);
        e[i].resize(3);
        en[i].resize(3);
        fn[i].resize(3);
    }
}

// OLD FORCE FUNCTION
void _ode_vars::InitExcitationVars(int _order, bool _rotatingwall, bool _antirotatingwall, double _U_exc, double _w_exc){
	for(unsigned int i=0;i<15;i++)
    {
        forcev.excitation_type[i] = false ;
    }
    //order should be 1,2 or 3;
    switch ( _order ){
        case 0:
            //the case of no-excitation
            break;
        case 1:
            forcev.excitation_type[0] =true; //dip
            break;
        case 2:
            forcev.excitation_type[1] =true; //quad
            break;
        case 3:
            forcev.excitation_type[2] =true;  // oct
            break;
        case 4:
            forcev.excitation_type[3] =true; // axial
            break;
        case 7:
            forcev.excitation_type[6] = true; //SIMCO
            break;
        case 8:
            forcev.excitation_type[7] =true;
            break;
	case 9:
	    forcev.excitation_type[8] = true; //AR
	    break;
	case 10:
	    forcev.excitation_type[9] = true; //AR
	    break;
        default:
            SLogger slogger("trapparameters");
            slogger << ERROR << "error in ode.cpp. Internally a wrong excitation order!" << SLogger::endmsg;exit(1);
    }
    forcev.excitation_type[4] =   _rotatingwall;
    forcev.excitation_type[5] =_antirotatingwall;
    U_exc = _U_exc;
    w_exc = _w_exc;
#ifdef __MPI_ON__
    MPI_Barrier(MPI::COMM_WORLD);
#endif // __MPI_ON__
}


double _ode_vars::SetSWIFT_function(string nf){
    ifstream file;
    double dum,dum2;
    bool first_line= true;
    double t0;
    file.open(nf.c_str(),ios::in);
    if(!file)
    {
        cout << "SWIFT FILE " << nf << " not found" << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        while(!file.eof())
        {
            file >> dum >> dum2;
            if(first_line)
            {
                t0 = dum;
                first_line =false;
            }
            swift_time.push_back(dum-t0);
            swift_amp.push_back(dum2);
            //cout << swift_time.back() << " " << swift_amp.back() << endl;
        }
        swift_time.pop_back();
        swift_amp.pop_back();
        file.close();
    }
    swift_nsample = swift_time.size();
    swift_dt = swift_time[1] - swift_time[0];
    //cout << swift_dt << endl;
    return swift_time.back();
    
}

double _ode_vars::SetSWIFT_function_RW(string nf){
    ifstream file;
    double dum,dum2,dum3;
    bool first_line= true;
    double t0;
    file.open(nf.c_str(),ios::in);
    if(!file)
    {
        cout << "SWIFT FILE " << nf << " not found" << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        while(!file.eof())
        {
            file >> dum >> dum2 >> dum3;
            if(first_line)
            {
                t0 = dum;
                first_line =false;
            }
            swift_time.push_back(dum-t0);
            swift_amp_sin.push_back(dum2);
            swift_amp_cos.push_back(dum3);
            //cout << swift_time.back() << " " << swift_amp.back() << endl;
        }
        swift_time.pop_back();
        swift_amp.pop_back();
        file.close();
    }
    swift_nsample = swift_time.size();
    swift_dt = swift_time[1] - swift_time[0];
    //cout << swift_dt << endl;
    return swift_time.back();
    
}


double _ode_vars::Interpolate_SWIFT(double t_){
    //linear interpolation
    int index = floor(t_/swift_dt);
    double res;
    if(index<swift_nsample)
    {
        res = swift_amp[index] + (t_ - swift_time[index])*(swift_amp[index+1]-swift_amp[index])/swift_dt;
    }
    else
    {
        res = swift_amp[swift_nsample-1];
    }
    return res;
}


void _ode_vars::UnsetSWIFT()
{
    swift_time.clear();
    swift_amp_sin.clear();
    swift_amp_cos.clear();
    swift_amp.clear();
}

double _ode_vars::Interpolate_SWIFT(double t_,int sc){
    //linear interpolation
    int index = floor(t_/swift_dt);
    double res;
    if(sc==0)
    {
        if(index<swift_nsample)
        {
            res = swift_amp_cos[index] + (t_ - swift_time[index])*(swift_amp_cos[index+1]-swift_amp_cos[index])/swift_dt;
        }
        else
        {
            res = swift_amp_cos[swift_nsample-1];
        }
    }
    else
    {
        if(index<swift_nsample)
        {
            res = swift_amp_sin[index] + (t_ - swift_time[index])*(swift_amp_sin[index+1]-swift_amp_sin[index])/swift_dt;
        }
        else
        {
            res = swift_amp_sin[swift_nsample-1];
        }
    }
    return res;
}

//  ************** _force_vars
_force_vars::_force_vars()
{
    
    for(unsigned i=0;i<15;i++){
        excitation_type[i] = false;} // dip, quad, oct , axial , RW, AW,SIMCO,SIMCOAxial
    scaledCoulombFactor =1.;
    coulombinteraction = false;

}

_force_vars::~_force_vars()
{
    
}
void _force_vars::Reset_excitation_type(){
    for(int i=0;i<15;i++)
        excitation_type[i] = false;
    return;
}
void _force_vars::Load_EXC_EMAP(string _file_name,double _factor){
    exc_emap = new fieldmap(400,400);
    exc_emap->ReadField(&_file_name[0]);
    exc_emap->SetFactor(_factor);
    exc_emap->SetRmin(-20);
    exc_emap->SetZmin(-20);
    return;
}

void _force_vars::Reset_EXC_EMAP(){
    delete exc_emap;
    return;
}
