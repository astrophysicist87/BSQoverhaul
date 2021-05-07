#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include <fstream>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <dirent.h>
#include <sstream>
#include <vector>
#include "particle.h"
#include "LinkList.h"
#include <errno.h>
#include <sys/stat.h>
#include <string>

using namespace std;

template <int D>
class Output
{
	private:
		ofstream ESC,FO,SC;
	   	string ESCname,FOname;
	   	string Scons;
	   	int countEP,countGC;
	   	string convertInt(int number);
	   	string cfolder;
	   	typedef struct stat Stat;
	   	
	public:
		Output<D>(LinkList<D> &linklist);
		void cleardir(LinkList<D> &linklist);
		void conservation(LinkList<D> &linklist);
		void eprofile(LinkList<D> &linklist);
		void eccout(LinkList<D> &linklist);
	    double ecc(LinkList<D> &linklist, double  & psi,double & rout , double & Etot, int m,int n);
		void sveprofile(LinkList<D> &linklist);
		void do_mkdir(const char *path, int mode);
		void averages(LinkList<D> &linklist);
		void FOstart(LinkList<D> &linklist);
		void FOprint(LinkList<D> &linklist);
		void svFOprint(LinkList<D> &linklist);
		void SCprint(LinkList<D> &linklist);
};

template <int D>
Output<D>::Output(LinkList<D> &linklist)
{
	countEP = 0;
	FOstart( linklist );
}


template <int D>
void Output<D>::FOstart(LinkList<D> &linklist)
{
	
   	
   	string FOin;
   	
   	if (linklist.average==1)
   	 {
   	 	FOname=cfolder+"average_freezeout.dat";
   	 }   	    	
   	 else if (linklist.fcount==0&&linklist.average!=1)
   	{
		if  (linklist.visc==3)
   	 	{
   	 		FOin="sbvfreezeout";
   	 	}
   	
   	 	FOname=ofolder+FOin+"_weight.dat";
   	 }
   	 else
   	 {
   	 	string under ("_ev");
   	 	string even;
   	 	even=convertInt(linklist.fnum);
   	 	if (linklist.visc==3)
   	 	{
   	 		FOin="sbvfreezeout";
   	 	}
   	 	
   	 	FOname=cfolder+FOin+under+even+".dat";
   	 }

  	FO.open(FOname.c_str());
  	if (!FO.is_open())
	{
		cout << "Error: cannot open E_S_conserv.dat file!" << endl;
		exit(1);
	}
	
  	FO.close();
  	
}


template <int D>
void Output<D>::eprofile(LinkList<D> &linklist)
{
	ofstream EPN;
	
   	countEP+=1;
   	string epin ("eprofile");
   	string dat (".dat");
   	string sCEP;
   	sCEP=convertInt(countEP);
   	string epname;
   	
   	if (linklist.average==1)
   	 {
   	 	epname=cfolder+"average_eprof"+ sCEP + ".dat";
   	 }   
   	else if (linklist.fcount==0&&linklist.average!=1)
   	{
   	 	epname=ofolder+epin+sCEP+dat;
   	 }
   	 else
   	 {
   	 	string under ("_ev");
   	 	string even;
   	 	even=convertInt(linklist.fnum);
   	 	epname=cfolder+epin+sCEP+under+even+dat;
   	 }
   	
   	//cout << epname << endl;
  	EPN.open(epname.c_str());
  	if (!EPN.is_open())
	{
		
		cout << "Error: cannot open eprofile" << sCEP <<".dat file!" << endl;
		exit(1);
	}
	else 
	{
		EPN << linklist.t << endl;
		
		
		for (int i=0;i<linklist.n();i++)
		{
    		EPN << linklist._p[i].r   << " " << linklist._p[i].EOS.e() <<  endl;	   	
		}
	}
	

	
  	EPN.close();
   	
}


template <int D>
void Output<D>::sveprofile(LinkList<D> &linklist)
{
	ofstream EPN;
	
   	countEP+=1;
   	string epin ("sveprofile");
   	string dat (".dat");
   	string sCEP;
   	sCEP=convertInt(countEP);
   	string epname;
   	
   	if (linklist.average==1)
   	 {
   	 	epname=cfolder+"average_eprof"+ sCEP + ".dat";
   	 }   
   	else if (linklist.fcount==0&&linklist.average!=1)
   	{
   	 	epname=ofolder+epin+sCEP+dat;
   	 }
   	 else
   	 {
   	 	string under ("_ev");
   	 	string even;
   	 	even=convertInt(linklist.fnum);
   	 	epname=cfolder+epin+sCEP+under+even+dat;
   	 }
   	
   	//cout << epname << endl;
  	EPN.open(epname.c_str());
  	if (!EPN.is_open())
	{
		
		cout << "Error: cannot open eprofile" << sCEP <<".dat file!" << endl;
		exit(1);
	}
	else 
	{
		EPN << linklist.t << endl;
		
		
		
		for (int i=0;i<linklist.n();i++)
		{
		EPN << linklist._p[i].r   << " " << linklist._p[i].EOS.e() << " " << linklist._p[i].EOS.T()*0.1973 << " " << linklist._p[i].v    <<endl;	
		}
	}
	

	
  	EPN.close();
   	
}

template <int D>
void Output<D>::eccout(LinkList<D> &linklist)
{
  ofstream EPN;
  
  countEP+=1;
  string epin ("eccCM");
  string dat (".dat");
  string sCEP;
  sCEP=convertInt(linklist.fnum);
  string epname;
  
  if (linklist.fcount==0&&linklist.average!=1)
    {
      epname=ofolder+epin+sCEP+dat;
    }
  else
    {
      string even=convertInt(linklist.fnum);
      epname=cfolder+epin+sCEP+dat;
    }
  
  cout << epname << endl;  
  EPN.open(epname.c_str());
  if (!EPN.is_open())
    {
      
      cout << "Error: cannot open eprofile" << sCEP <<".dat file!" << endl;
      exit(1);
    }
  
  double psi,rout,etot;
  double e22=ecc(linklist,psi,rout,etot,2,2);
  EPN << convertInt(linklist.fnum) << " " << etot << " " ;
  EPN << e22 << " " << psi << " "  ;
  double e33=ecc(linklist,psi,rout,etot,3,3);
  EPN << e33 << " " << psi  << " "  ;
  double e44=ecc(linklist,psi,rout,etot,4,4);
  EPN << e44 << " " << psi  << " "  ;
  double e55=ecc(linklist,psi,rout,etot,5,5);
  EPN << e55 << " " << psi  << " "  ;
  double e66=ecc(linklist,psi,rout,etot,6,6);
  EPN << e66 << " " << psi  << " "  ;
  
  EPN << rout << endl;
  EPN.close();
  
}

template <int D>
double Output<D>::ecc(LinkList<D> &linklist, double  & psi,double & rout , double & Etot, int m,int n){

  

  int max=linklist.n();
  
  double xcm=0,ycm=0,etot=0;
  for (int i=0;i<max;i++){
   
    xcm+=linklist._p[i].r.x[0]*linklist._p[i].EOS.e();
    ycm+=linklist._p[i].r.x[1]*linklist._p[i].EOS.e();
    etot+=linklist._p[i].EOS.e();
    
  }
  xcm/=etot;
  ycm/=etot;

  Etot=etot/max;
  
  vector <double> r2,phi;
  r2.resize(max);
  phi.resize(max);
  double psit=0,psib=0,rb=0;
  for (int s=0;s<max;s++){
    
    double xsub=(linklist._p[s].r.x[0]-xcm);
    double ysub=(linklist._p[s].r.x[1]-ycm);
    r2[s]=xsub*xsub+ysub*ysub;
    phi[s]=atan2(ysub,xsub);
    double rv=linklist._p[s].EOS.e()*pow(r2[s],(m/2.));
    psit+=rv*sin(1.0*n*phi[s]);
    psib+=rv*cos(1.0*n*phi[s]);
    
    rb+=rv;
    
  }
  psit/=max;
  psib/=max;
  
  psi=1./(1.0*n)*atan2(psit,psib);
  double ec=0;
  for (int s=0;s<max;s++) ec+=linklist._p[s].EOS.e()*pow(r2[s],m/2.)*cos(n*(phi[s]-psi)); 

  ec/=rb;
  rout=rb/etot;
  
  return ec;

}







template <int D>
string Output<D>::convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

template <int D>
void Output<D>::conservation(LinkList<D> &linklist)
{ 	


  	ESC.open(ESCname.c_str(), ios::out | ios::app );
	ESC << linklist.t << " " <<  linklist.Eloss << endl;
  	ESC.close();
  	
  	if (ESC.is_open())
  	{
  		cout << "error: still open" << endl;
 		exit(1);
  	}
  	
  	
}


template <int D>
void Output<D>::FOprint(LinkList<D> &linklist)
{ 	

  	FO.open(FOname.c_str(), ios::out | ios::app );
  	for (int i=0;i< linklist.cf;i++)
  	{
	FO <<  linklist.divTtemp[i] << " " << linklist.divT[i] << " " << linklist.gsub[i] << " " << linklist.uout[i] << " " << linklist.swsub[i] << " " <<   linklist.tlist[i] <<  " " << linklist.rsub[i] << " " << linklist.sFO[i] <<  " " << linklist.Tfluc[i] <<  endl;
	}
  	FO.close();
  	
  	delete [] linklist.divTtemp;
  	delete [] linklist.divT;
  	delete [] linklist.gsub;
  	delete [] linklist.uout;
  	delete [] linklist.swsub;
  	delete [] linklist.rsub;
  	delete [] linklist.tlist;
  	
  	if (FO.is_open())
  	{
  		cout << "error: still open" << endl;
 		exit(1);
  	}
  	
  	
}


template <int D>
void Output<D>::svFOprint(LinkList<D> &linklist)
{ 	

  	FO.open(FOname.c_str(), ios::out | ios::app );
  	for (int i=0;i< linklist.cf;i++)
  	{
  	FO <<  linklist.divTtemp[i] << " " << linklist.divT[i] << " " << linklist.gsub[i] << " " << linklist.uout[i] << " " << linklist.swsub[i] << " " <<  linklist.bulksub[i] << " "  <<  linklist.shearsub[i].x[0][0] << " "  <<  linklist.shearsub[i].x[1][1] << " "  <<  linklist.shearsub[i].x[2][2] <<" " <<linklist.shear33sub[i] << " "  <<  linklist.shearsub[i].x[1][2] <<  " "  <<  linklist.tlist[i] <<  " " << linklist.rsub[i] <<  " " << linklist.sFO[i] <<   " " << linklist.Tfluc[i] <<endl; // removed entropy, added in tau
	}
  	FO.close();
  	
  	delete [] linklist.divTtemp;
  	delete [] linklist.divT;
  	delete [] linklist.gsub;
  	delete [] linklist.uout;
  	delete [] linklist.swsub;
  	delete [] linklist.bulksub;
  	delete [] linklist.shearsub;
  	delete [] linklist.rsub;
  	delete [] linklist.tlist;
  	delete [] linklist.shear33sub;
  	
  	if (FO.is_open())
  	{
  		cout << "error: still open" << endl;
 		exit(1);
  	}
  	
  	
}


template <int D>
void Output<D>::averages(LinkList<D> &linklist)
{ 	

	ofstream AVG;
	string avgname=cfolder+"anisotropy.dat";
	
	AVG.open(avgname.c_str());
  	if (!AVG.is_open())
	{
		
		cout << "Error: cannot open anisotropy.dat file!" << endl;
		exit(1);
	}
	else 
	{
  		for (int i=0;i<linklist.steps;i++)
  		{
		AVG << linklist.tsave[i] << " " <<  linklist.pan[i] <<  " " <<  linklist.san[i] << endl;
		}
  	}
  	
  	AVG.close();
  	
  	delete [] linklist.pan;
  	delete [] linklist.san;
  	
}

template <int D>
void Output<D>::cleardir(LinkList<D> &linklist)
{

    
    DIR *pdir;
    struct dirent *pent;
    if (linklist.fcount==0)
    {
    pdir=opendir(ofolder.c_str());
    }
    else
    {
    pdir=opendir(cfolder.c_str());
    }
    if (!pdir)
    {
    cout<<"Error: " << cfolder << " does not exist";
    }
    int errno3=0; //errno.h
    while ((pent=readdir(pdir)))
    {
    //cout<<pent->d_name;
    string file_delete=pent->d_name;
    file_delete= ofolder+file_delete;
//    DeleteFile(file_delete.c_str());
    std::remove(file_delete.c_str());
    }
    if (errno3)
    {
    cout<<"Error while accessing directory";
    }
    closedir(pdir);
    
}

template <int D>
void Output<D>::do_mkdir(const char *path, int mode)
{
    Stat            st;

    if (stat(path, &st) != 0)
    {
	cout << "Directory path doesn't exist" << endl;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        cout << "Couldn't make the directory for output files" << endl;
        exit(1);
    }

}

#endif
