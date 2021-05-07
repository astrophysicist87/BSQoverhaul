#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

using namespace std;

#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "rungekutta4.h"
#include "hydrosim.h" 
#include "eos.h"
#include "output.h"


void Simulation(double dt,LinkList<2> &linklist)
{
	cout << "Ready to start hydrodynamics\n";
	
	
	 linklist.frzc=0;
	 
	 linklist.cf=0;
	 	
	
	Output<2> out(linklist); // sets up output
	
	linklist.t=linklist.t0;
	
	linklist.Ez=0;
	
	if ((linklist.qmf==1)||(linklist.qmf==3)) {
	out.eprofile(linklist);
	cout << "printed first timestep" << endl;
	}
	else if (linklist.qmf==2){
	out.eprofile(linklist);
	cout << "printed first timestep" << endl;
	exit(1);
	}
	else if (linklist.qmf==4){
	  out.eccout(linklist);
          cout << "eccentricity printed" << endl;
        }
	

	while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {
		
		 linklist.cfon=1;
 	    rungeKutta2<2>(dt,&idealhydro3<2>,linklist);  // slower method, more accurate	
 	    
	    if (linklist.cf>0) out.FOprint(linklist); // prints out frozen out SPH particles
	    
 	    //if you add more points to print off you must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
 
	    if (linklist.qmf==3){
	    double tsub=linklist.t-floor(linklist.t);
 	    if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check 	    
 	    {
		cout << "t=" << linklist.t <<endl;  // outputs time step
 	    	out.eprofile(linklist);   // energy density profile
 	    	
 	    } 	
 	    else if ((tsub<(0.5+dt*0.5))&&(tsub>=(0.5-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
 	    {
		cout << "t=" <<  linklist.t <<endl;  // outputs time step
 	    	out.eprofile(linklist);   // energy density profile
 	    	
 	    }   

	    }
  	  
  	  
	}

	cout << "done" << endl;
	linklist.endEV();
}

void svSimulation(double dt,LinkList<2> &linklist)
{
	cout << "Ready to start hydrodynamics\n";

	
	 int cc=0;
	 linklist.frzc=0;
	 linklist.cf=0;

	
	Output<2> out(linklist);
	
	linklist.t=linklist.t0;
	
	if (linklist.qmf==1||linklist.qmf==3){
	  out.sveprofile(linklist);
	  cout << "printed first timestep" << endl;
	  linklist.conservation_entropy();
	  cout << "t=" << linklist.t << " S=" << linklist.S << endl; 
	
	} 
	else if(linklist.qmf==4){
          out.eccout(linklist);
	  cout << "eccentricity printed" << endl;
	  exit(0);
        }
	
	linklist.Ez=0;
	while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {
		linklist.cfon=1;
		
 		svrungeKutta2<2>(dt,&shear<2>,linklist);
	   
 	    if (linklist.cf>0) out.svFOprint(linklist);
 	    
 	    if (linklist.qmf==3){
 	    double tsub=linklist.t-floor(linklist.t);
	    cout << "t=" << linklist.t <<endl;                                                                         
	          out.sveprofile(linklist);
	    }
  	  
	}

	
	linklist.endEV();
}


template <int D> 
void idealhydro3(LinkList<D>  &linklist) // ideal Equations of motion, only set up completely for 2+1 at the moment
{
     
     
     linklist.initiate();  // initiates linklist
     
     	
     	for(int i=0; i<linklist.n();i++)
	{
		 linklist.optimization(i); // calculates entropy density
	}
	 
	int curfrz=0;	
	for(int i=0; i<linklist.n(); i++) 
	{
                 //  Computes gamma and velocity
                 
        	linklist._p[i].calc(linklist.t); // calculates gamma, flow vectors, and updates EOS
        	
        	linklist._p[i].returnA(); // returns A form EOM
//        	if(linklist._p[i].EOS.s() < 0)
        	if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n()); // checks if SPH particles have frozen out
        	
		
	}
	
	if (linklist.cfon==1) // creates list of frozen out SPH particles
	{
	linklist.number_part+=curfrz;	
	linklist.list.resize(curfrz);
	}
	if (linklist.rk2==1) linklist.conservation();	// calculates conservation of energy, only on first Runge kutta step
			
	int m=0;
	for(int i=0; i<linklist.n();i++)
	{	
		//      Computes gradients to obtain dsigma/dt
		linklist.optimization2(i);
		linklist._p[i].sigset(linklist.t); // returns another value for EOM
		if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1)) // if SPH particles have remained frozen out after two time steps, allows to freezeout
		{
			linklist.list[m]=i;
			linklist._p[i].Freeze=4;
			++m;
		}
		
	}
	
	linklist.conservation_Ez(); // calculates dEz for conservation of Energy
	if (linklist.first==1) //checks if t=t0, if so sets the initial energy E0
	{
		linklist.conservation_entropy();
		linklist.conservation_E();				
	}
	
	
	
		//calculate matrix elements		
	for(int i=0; i<linklist.n();i++)
	{	
	
		// set the Mass and the Force matrix
		double *M,*F;
     		M=new double[D*D];
     		F=new double[D];
		M[0]=linklist._p[i].Agam*linklist._p[i].u.x[0]*linklist._p[i].u.x[0]+linklist._p[i].EOS.w()*linklist._p[i].gamma;
		M[3]=linklist._p[i].Agam*linklist._p[i].u.x[1]*linklist._p[i].u.x[1]+linklist._p[i].EOS.w()*linklist._p[i].gamma;
		M[1]=linklist._p[i].Agam*linklist._p[i].u.x[0]*linklist._p[i].u.x[1];
		
		F[0]=linklist._p[i].Agam2*linklist._p[i].u.x[0]-linklist._p[i].gradP.x[0];
		F[1]=linklist._p[i].Agam2*linklist._p[i].u.x[1]-linklist._p[i].gradP.x[1];
		
	
		
		
		double det=M[0]*M[3]-M[1]*M[1]; // inverts matrix
		double MI[4];
		MI[0]=M[3]/det;
		MI[1]=-M[1]/det;
		MI[3]=M[0]/det;
		linklist._p[i].du_dt.x[0]=F[0]*MI[0]+F[1]*MI[1];
		linklist._p[i].du_dt.x[1]=F[0]*MI[1]+F[1]*MI[3];
		
		
		
		linklist._p[i].div_u = (1./ linklist._p[i].gamma)*inner( linklist._p[i].u, linklist._p[i].du_dt) - ( linklist._p[i].gamma/ linklist._p[i].sigma)* linklist._p[i].dsigma_dt ;
		
		delete [] M;
		delete [] F;
			
	}
	if (linklist.cfon==1) linklist.freezeout(curfrz); // calculates normals from freezeout hypersurface
	
	linklist.destroy();
 }
 

 template <int D> 
void shear(LinkList<D>  &linklist)  // shear+bulk Equations of motion, only set up completely for 2+1 at the moment
{
     linklist.setshear();
     linklist.initiate();
     
     	for(int i=0; i<linklist.n();i++)
	{
		linklist.optimization(i);
		

		//		if ((linklist._p[i].eta<0)||isnan(linklist._p[i].eta)) 
		if (linklist._p[i].eta<0)
		{
			cout << i <<  " neg entropy " <<  linklist._p[i].EOS.T()*197.3   << " " << linklist._p[i].eta << endl;
			
			linklist._p[i].eta=0;

		}
		
	}
	
	int curfrz=0;
	for(int i=0; i<linklist.n(); i++) 
	{
                 //  Computes gamma and velocity
                 
        	linklist._p[i].calc(linklist.t);
		linklist._p[i].setvisc(linklist.etaconst,linklist.bvf,linklist.svf,linklist.zTc,linklist.sTc,linklist.zwidth,linklist.visc);
		if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n());
        		
	}
	if (linklist.cfon==1)
	{
	linklist.number_part+=curfrz;		
	linklist.list.resize(curfrz);
	}
	
	int m=0;	
	for(int i=0; i<linklist.n();i++)
	{	
		//      Computes gradients to obtain dsigma/dt
		linklist.svoptimization2(i,linklist.t,curfrz);
		linklist._p[i].dsigma_dt = -linklist._p[i].sigma*(linklist._p[i].gradV.x[0][0]+linklist._p[i].gradV.x[1][1]) ;
		
		linklist._p[i].svsigset(linklist.t,i);
		if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1)) 
		{
			linklist.list[m]=i;
			linklist._p[i].Freeze=4;
			++m;
		}
			
	}
	
	if (linklist.rk2==1) linklist.svconservation();
	linklist.svconservation_Ez();
	
	
		//calculate matrix elements		
	for(int i=0; i<linklist.n();i++)
	{	
		double gamt=1./linklist._p[i].gamma/linklist._p[i].stauRelax;
		double pre=linklist._p[i].eta_o_tau/2./linklist._p[i].gamma;
		//p4=gamt-linklist._p[i].sigl*4./3.;
		double p1=gamt-4./3./linklist._p[i].sigma*linklist._p[i].dsigma_dt+1./linklist.t/3.;
		Vector<double,D>  minshv=rowp1(0,linklist._p[i].shv);		
		//p2=linklist._p[i].setas*gamt;
		Matrix <double,D,D> partU=linklist._p[i].gradU+transpose(linklist._p[i].gradU);

		// set the Mass and the Force
		Matrix <double,D,D> M=linklist._p[i].Msub(i);
		Vector<double,D> F=linklist._p[i].Btot*linklist._p[i].u+ linklist._p[i].gradshear -(linklist._p[i].gradP+linklist._p[i].gradBulk+ linklist._p[i].divshear);
		// shear contribution
		F+=pre*linklist._p[i].v*partU+p1*minshv;
		
		
		double det=deter(M);
		Matrix <double,D,D> MI;
		MI.x[0][0]=M.x[1][1]/det;
		MI.x[0][1]=-M.x[0][1]/det;
		MI.x[1][0]=-M.x[1][0]/det;
		MI.x[1][1]=M.x[0][0]/det;
		linklist._p[i].du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
		linklist._p[i].du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];
		
		Matrix <double,D,D> ulpi=linklist._p[i].u*colp1(0,linklist._p[i].shv);

		double vduk=inner(linklist._p[i].v,linklist._p[i].du_dt);

		Matrix <double,D,D> Ipi=-linklist._p[i].eta_o_tau/3.*(linklist._p[i].Imat+linklist._p[i].uu)+4./3.*linklist._p[i].pimin;
		
        	linklist._p[i].div_u = (1./ linklist._p[i].gamma)*inner( linklist._p[i].u, linklist._p[i].du_dt) - ( linklist._p[i].gamma/ linklist._p[i].sigma)* linklist._p[i].dsigma_dt ;
		linklist._p[i].bigtheta=linklist._p[i].div_u*linklist.t+linklist._p[i].gamma;
		
		
		Matrix <double,D,D> sub=linklist._p[i].pimin+linklist._p[i].shv.x[0][0]/linklist._p[i].g2*linklist._p[i].uu-1./linklist._p[i].gamma*linklist._p[i].piutot;
		
		
		linklist._p[i].inside=linklist.t*(inner((-minshv+linklist._p[i].shv.x[0][0]*linklist._p[i].v),linklist._p[i].du_dt)- con2(sub,linklist._p[i].gradU)    -linklist._p[i].gamma*linklist.t*linklist._p[i].shv33);
		linklist._p[i].detasigma_dt =1./linklist._p[i].sigma/linklist._p[i].EOS.T()*( -linklist._p[i].bigPI*linklist._p[i].bigtheta+linklist._p[i].inside);
		
		
		
		
        	linklist._p[i].dBulk_dt = (-linklist._p[i].zeta/linklist._p[i].sigma*linklist._p[i].bigtheta - linklist._p[i].Bulk/linklist._p[i].gamma )/linklist._p[i].tauRelax;
        	
        	Matrix <double,D,D> ududt=linklist._p[i].u*linklist._p[i].du_dt;
         	 

        	linklist._p[i].dshv_dt= -gamt*(linklist._p[i].pimin+linklist._p[i].setas*0.5*partU)-0.5*linklist._p[i].eta_o_tau*(ududt+transpose(ududt))+linklist._p[i].dpidtsub()-vduk*(ulpi+transpose(ulpi)+(1/linklist._p[i].gamma)*Ipi)+linklist._p[i].sigl*Ipi;
		
	}
	
	
	if (linklist.cfon==1) linklist.svfreezeout(curfrz);
	
	
	linklist.destroy();
	
 }


