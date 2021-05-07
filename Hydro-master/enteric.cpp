#ifndef _ENTERIC_CPP_
#define _ENTERIC_CPP_

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <vector>
#include "mathdef.h"
#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "enteric.h"
#include<dirent.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>

using namespace std;

void readEOStable()
{

    /////////////////////////////////////////////////////
    ////       Start reading EoS table           /////
    /////////////////////////////////////////////////////
  	

  	cout << "The EOS should be contained in two files (either .txt or .dat files).  The first contains the Temperature (in MeV's) and the first line should include the number of steps.\n";
  	cout << "Enter file name of the first file:\n";
	cin.clear();
	cin.ignore (cin.rdbuf( )->in_avail( ));
	string filen;
  	getline (cin,filen);  
  	
  	readEOS_T(filen); 
  	
  	cout << "The second EOS (either .txt or .dat files) contains the  Energy Density, Pressure, Entropy,   and c_s^2 (in MeV's) and the first line should include the number of steps.\n";
  	cout << "Enter file name of the first file:\n";
	cin.clear();
        cin.ignore (cin.rdbuf( )->in_avail( ));
        string filen2;
  	getline (cin,filen2);  
  	
  	readEOS_p(filen2);
  	
    /////////////////////////////////////////////////////
    ////       Finished reading EoS table           /////
    /////////////////////////////////////////////////////
}


void readEOS_T(string &firstry)
{
	string filename;
	filename = ifolder+firstry;
	//nameenter:
	
        FILE * myfile = fopen (filename.c_str(),"r");
        
        
        if ((!firstry.empty())&&(myfile!= NULL)){
  	   fscanf(myfile,"%i \n",&nETH);
           
           ETH=new eostableshigh[nETH];
           int i=0;
          while (i<nETH ){
                   fscanf(myfile,"%lf     %*f  %*f \n",&ETH[i].T); //assumes T is in the formate MeV
                   
                   ETH[i].T/=197.3;                 
                   ++i;
               }
  	  fclose(myfile);
	 // cout << filename.c_str() << ": Input sucessful!\n";
  	}
  	else {
  		cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";
		
		exit(1);	
  	}	
  		

}


void readEOS_p(string &firstry)
{
	 
	string filename;
	filename = ifolder+firstry;
	//nameenter:
        FILE * myfile = fopen (filename.c_str(),"r");	
  	if ((!firstry.empty())&&(myfile!= NULL)){
  	   fscanf(myfile,"%*i \n");
           int i=0;
          while (i<nETH ){
                   fscanf(myfile,"%lf     %lf    %lf    %lf \n",&ETH[i].e,&ETH[i].p,&ETH[i].s,&ETH[i].dtds); //assumes that e,p is in the format GeV/fm^3 and s in 1/fm^3
                   ETH[i].e/=0.1973;
                   ETH[i].p/=0.1973;
                   ++i;
               }
  	  fclose(myfile);
	 // cout << filename.c_str() << ": Input sucessful!\n";
  	}
  	else {
  		cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";
		
		exit(1);	
  	}	
}

void readEOS_lowS(string &firstry)
{

	
	string filename;
	filename = ifolder+firstry;
	//nameenter:
        FILE * myfile = fopen (filename.c_str(),"r");
        
        if ((!firstry.empty())&&(myfile!= NULL)){
  	   fscanf(myfile,"%i \n",&nETL);
           ETL = new eostableslow[nETL];
           int i=0;
          while (i<nETL ){
                   fscanf(myfile,"%lf     %lf  %lf    %lf   %lf\n",&ETL[i].T,&ETL[i].e,&ETL[i].p,&ETL[i].s,&ETL[i].dtds); //assumes T is in the formate MeV
                   ETL[i].T/=197.3;
                   i++;
               }
  	  fclose(myfile);
	//  cout << filename.c_str() << ": Input sucessful!\n";
  	}
  	else {
  		cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";
		
		exit(1);
  	}	
  		

}


//2 dimensions

//trento
void readICs_tnt(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& sfcheck, int & numpart, eos EOS)
{
  
  string filename;
  filename = ifolder+firstry;
  //nameenter:
  ifstream input(filename.c_str());
	if (!input.is_open())
 	{
 	cout << "Can't open " << filename << endl;
 	exit(1);
 	}

	string line;
	vector<double> xsub,ysub,esub;
	
	getline(input,line);
	std::vector<std::string> gx = split(line, ' ');
	
	double stepx,stepy;
	stringstream s;
	s << gx[1];
	s >> stepx;
	
	
	stringstream s1;
	s1 << gx[2];
	s1 >> stepy;
	
	cout << "dx=dy=" << stepx << " " << stepy << endl;
	
	
	while (getline(input,line)) {
	std::vector<double> y (3,0) ;
	
	std::vector<std::string> x = split(line, ' ');
	
	
	for(int j=0;j<3;j++)
	{
	stringstream ss;
	ss << x[j];
	ss >> y[j];
	}
	
	if ((factor*y[2])>0.05){
		xsub.push_back(y[0]);
		ysub.push_back(y[1]);
		esub.push_back(y[2]);
	}
	
	}
	input.close();
         
    
    _Ntable3=xsub.size();
    _p= new Particle<2>[_Ntable3];
          
    cout << "After e-cutoff=" << _Ntable3 << endl;
          
          
    int kk=_Ntable3;
    numpart=0;
          
    
        
    for(int j=0;j<_Ntable3;j++){
      _p[j].r.x[0]=xsub[j];
      _p[j].r.x[1]=ysub[j];
      // _p[j].e_sub=EOS.e_out(factor*esub[j]);
      _p[j].s_an=factor*esub[j];
      _p[j].u.x[0]=0;
      _p[j].u.x[1]=0;        
      _p[j].eta_sigma  = 1;
      _p[j].sigmaweight=stepx*stepy;
      _p[j].Bulk = 0;
                 
       

      if (_p[j].s_an>sfcheck)
	{
	  _p[j].Freeze=0;
	}
      else
	{
	  _p[j].Freeze=4;
	  --kk;
	  ++numpart;
	}    
    }
          
    cout << "After freezeout=" << _Ntable3-numpart << endl;
        
  
}


void readICs_tnt(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor,double const& efcheck, int & numpart, eos EOS)
{
}

void open_dir(string eventfolder, string *&filelist, int &count)
{

    string dotslash ("./");
    string dash ("/");
    string dir_in=dotslash+ifolder+("/")+eventfolder+("/");
    string dat("dat");
    string d1 ("dens1.dat");
    string d2 ("dens2.dat");
    DIR *pdir;
    count=0;
    struct dirent *pent;
    pdir=opendir(dir_in.c_str());
    if (!pdir)
    {
    	cout<<"Directory doesnot exist";
    }
    int errno4=0; //errno.h
    while ((pent=readdir(pdir)))
    {
    	string file_open=pent->d_name;
    	int j=file_open.length();
    	if (j>4)
    	{
    	string datcheck=file_open.substr((j-3),3);
    	if (dat==datcheck&&file_open!=d1&&file_open!=d2) count++;
    	}
    
    	
    }
    closedir(pdir);
    
    filelist=new string[count];
    pdir=opendir(dir_in.c_str());
    int i=0;
    while ((pent=readdir(pdir)))
    {
    	string file_open=pent->d_name;
    	int j=file_open.length();
    	
    	if (j>4)
    	{
    	string datcheck=file_open.substr((j-3),3);
    	if (dat==datcheck&&file_open!=d1&&file_open!=d2)
    	{
    	filelist[i]= eventfolder+dash+file_open;
    	i++;
    	}
    	}
    }
    if (errno4)
    {
    	cout<<"Error while accessing directory";
    }
    closedir(pdir);
    
}


string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {    	
        if (!item.empty()) elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

#endif
