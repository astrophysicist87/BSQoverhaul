#ifndef ENTERIC_H_ 
#define ENTERIC_H_

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
 
#include "mathdef.h"
#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "LinkList.h"
#include "eostables.h"
#include "eos.h"
#include <stdlib.h>
#include <vector>

using namespace std;

void readEOStable();
void readEOS_lowS(string &firstry);


void readEOS_T(string &firstry);

void readEOS_p(string &firstry);

void open_dir(string eventfolder, string *&filelist,int &count);

string convertInt(int number);

std::vector<std::string> split(const std::string &s, char delim);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

// 2 dimension functions
void readICs_tnt(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& efcheck, int & numpart, eos EOS); //trento (entropy density)

// 3 dimension functions
void readICs_tnt(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor,double const& efcheck, int & numpart, eos EOS); //trento (entropy density)



template <int D>
void manualenter(_inputIC &ics, LinkList<D> &linklist)
{
   
   string manf = ifolder+ics.man;
   istringstream manis (manf.c_str());
   //manenter:
   FILE * openmanf = fopen (manf.c_str(),"r");	
    double h,factor;
   double it0;
   int start,end; 
   string eos_s,eos_p,ic,outf,ictype,ic2,eostype,fvisc,low_switch;
   string table ("table");
   string trento ("trento");
   int df;
  	if ((!manf.empty())&&(openmanf!= NULL)){
  	   fscanf(openmanf,"%*s %lf	%*s   	%lf	    \n",&h,&ics.dt);
  	   char charin[100],charin2[100];
  	   fscanf(openmanf,"%*s %s %*s  %s \n",charin2,charin);  
  	   fvisc=charin2;   
  	   eostype=charin;
  	   if (eostype==table)
  	   {
  	   	fscanf(openmanf,"%s \n",charin);  	   
  	   	eos_s=charin;
  	   	fscanf(openmanf,"%s \n",charin);  	   
  	   	eos_p=charin;
  	   	char charin3[100];
  	   	fscanf(openmanf,"%*s %s\n",charin3);  	   
  	   	low_switch=charin3;
  	   }
  	   fscanf(openmanf,"%*s %i \n",&linklist.etaconst); //if 1, eta/s const, else temp dependent
  	   if (linklist.etaconst!=1) fscanf(openmanf,"%*s %lf  %lf  %lf \n",&linklist.zwidth,&linklist.sTc,&linklist.zTc); // guassian width zeta/s, eta/s T_c, zeta/s T_c
	   fscanf(openmanf,"%*s %lf %lf \n",&linklist.bvf,&linklist.svf);
	   fscanf(openmanf,"%*s %lf \n",&it0);  
	   fscanf(openmanf,"%*s %lf \n",&freezeoutT); 		    
	   freezeoutT/=197.3; 
	   
  	   fscanf(openmanf,"%s \n",charin);  	   
  	   ictype=charin;
  	   
  	   
  	   hcor=0;
  	   
  	
  	  
  	   
  	   if (ictype==trento)
  	   {
  	   	fscanf(openmanf,"%lf \n",&factor);  
  	   	linklist.factor=factor;
  	   	
        }
  	   
  	   fscanf(openmanf,"%*s %s \n",charin);  	   
  	   ic=charin;
  	    int qmin;
  	   fscanf(openmanf,"%*s %i \n",&qmin);  
  	   linklist.qmf=qmin;  
  	   
  	   fscanf(openmanf,"%*s %s \n",charin);  	   
  	   outf=charin;
  	   
  	   if (ictype==trento)
  	   {
  	   	fscanf(openmanf,"%d \n",&start);  
  	   	fscanf(openmanf,"%d \n",&end); 
  	   	
  	   	
  	   	if (ics.on==0){
  	   	linklist.start=start; 
  	   	linklist.end=end;}
  	   	else
  	   	{
  	   	linklist.start=ics.start; 
  	   	linklist.end=ics.end;
  	   	}
  	   	cout << "Event # " << linklist.start << " with fac="<< factor << endl;
           }
           
	   //           cout << "Event # " << linklist.start << " with fac="<< factor << endl;           
           fscanf(openmanf,"%*s %i \n",&df);  
           if (df!=0) // for df, analytical solution - no decays
  	   {
  	   	fscanf(openmanf,"%*s %s \n",charin);
  	   	string dfpre,dffile=charin;

        // Added by Christopher Plumberg - 04/14/2021
        // If df < 0, assume df is one level higher in directory structure
        string directoryPrefix = (df < 0) ? "../" : "";
  	   	if (abs(df)==1) dfpre="df/input/";  
  	   	else if (abs(df)==2) dfpre="sampling/input/"; 
  	   	else cout << "Error: undefined df calcualtion type " << df <<  endl;
  	   	string dfout=directoryPrefix+dfpre+dffile;
  	   	
  	   	
  	   	
  	   	ofstream DFO;
  	   	DFO.open(dfout.c_str());
  	   	if (!DFO.is_open())
		{cout << "Error: cannot open "<< dfout << " file!" << endl;
		getchar();}
		
		DFO << "typeofequations: " << fvisc <<  endl;
		DFO << "folder: " << outf << endl;
		DFO << "range(pt,phi): input/gl15.dat input/gq20.dat" << endl;
		DFO << "decays: 1" << endl;
		DFO << "rangeofevents: " << linklist.start << " " << linklist.end << endl;
  	   	DFO << "freezeouttemp: " << freezeoutT*0.1973 << endl;
  	   	
  	   	
  	   	while (fgets (charin, 100 ,openmanf)!= NULL )  DFO << charin ;
  	   	
  	   	DFO.close();
  	   	
 
           }
          
         
        
  	   fclose(openmanf);
	   cout << manf.c_str() << ": Input sucessful!\n";
  	}
  	else {
  		
  	
  		cout << "Error: " << manf.c_str() << " does not exist.  Please enter new a file name that includes h,timestep,dimensions and file names for EOS/IC's \n";
		
		exit(1); 	
  	}

  	linklist.setv(fvisc);
  	linklist.eost=eostype;
  	linklist.cevent=0;
//    cout << "h=" << h << " timestep=" << timestep << " dimensions=" << D << " dt=" << dt << " Number of Loops=" << Nloop << " Output loops=" << Outnum << "\n";
    cout << fvisc << " hydro, h=" << h <<  " dimensions=" << D << " dt=" << ics.dt << " QM fluc:  " <<  linklist.qmf << "\n";
    
    if (eostype==table)
    {
    	cout << "Using Equation of State table from: " << eos_s << " and " << eos_p << "\n";
    	
    	//       Start reading EoS table           //
    	string on ("on");   	
		
    	readEOS_T(eos_s);  
    	readEOS_p(eos_p);
    	
    	if (low_switch==on)
    	{
    		
    		linklist.lowT=1;
    		string lowsT ("lowsT.dat");
    		readEOS_lowS(lowsT);
    	}
    	else
    		linklist.lowT=0;
    }
    else
    {
    	cout << "Using " << eostype << " Equation of State." << endl;
    }
   
   
    eos EOS;
    EOS.eosin(eostype);   
    double efcheck=EOS.efreeze();
    double sfcheck=EOS.sfreeze();
  
    linklist.efcheck=efcheck;    
    linklist.sfcheck=sfcheck; 
    linklist.fcount=0;    
    linklist.average=0;
    //       Start reading ICs          //    
    Particle<D> *_p;
    int numpart,_Ntable3;
    
    linklist.gtyp=0;
	if (ictype==trento)
	{
		
		int count=linklist.end-linklist.start+1;
		linklist.ebe_folder=outf;
		string *filelist;
		filelist=new string[count];
		
		int j=0;
		for (int i=linklist.start;i<=linklist.end;i++)
		{
			filelist[j]= ic+"/ic"+convertInt(i)+".dat";
			j++;
		}
		linklist.filenames=filelist;
		linklist.fcount=count;
		linklist.fnum=linklist.start;
		
		
		readICs_tnt(linklist.filenames[0],  _Ntable3, _p,factor,sfcheck,  numpart, EOS);
		
		_p[0].start(eostype);
		linklist.setup(it0,_Ntable3,h,_p,ics.dt,numpart);
		
		cout << "number of sph particles=" << _Ntable3 << endl;
		linklist.gtyp=5;
		
	}

    string name3=ofolder+linklist.ebe_folder+"/temp_ic"+convertInt(linklist.start)+".dat";
    ofstream PRINT(name3.c_str());
    if (!PRINT.is_open())
      {
	cout << "Can't open " << name3 << endl;
	exit(1);
      }
    
    PRINT << _Ntable3-numpart << endl;
    
    PRINT.close();

   
 // cout << "Setting up EOS" << endl;
    
   // sets up EOS
	for (int i=0;i<linklist.n();i++)
	{	
		linklist._p[i].start(eostype);
	}
	
   
	if (ictype==trento)
   {
   	
   	linklist.updateIC();
   }

   linklist.freezeset();
   
   if ((fvisc=="bulk+shear")||(fvisc=="shear+bulk"))  
   {
	   linklist.sv_set();
   }
    
    linklist.first=1;
    
}




template <int D>
void nextevent(int i, LinkList<D> &linklist)
{

	int numpart, _Ntable3;
	Particle<D> *_p;
	if (linklist.gtyp==5) readICs_tnt(linklist.filenames[i],  _Ntable3, _p,linklist.factor,linklist.sfcheck,numpart,linklist._p[0].EOS);
    linklist.setupnext(_Ntable3,_p,numpart);
    	
    		
	for (int i=0;i<linklist.n();i++)
	{
		linklist._p[i].start(linklist.eost);
	}
	
    
	linklist.updateIC();
	

  	 if (linklist.visc==1) linklist.etas_set();
  	 
  	 
  	 linklist.first=1;
  	 linklist.fnum+=1;
	
}


#endif
