
#ifndef _LINKLIST_H_
#define _LINKLIST_H_

#include "vector.h"
#include "matrix.h"
#include "particle.h"
#include "mathdef.h"
#include "random.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

template <int D>
class LinkList {
private:
	int _n;
	
	Vector<double,D> min;   	
   	Vector<double,D> max;	
   	int range; //range is number of boxes from left ot right extra	   	
	int Size;
	static constexpr double e0=1.;
	static constexpr double q=1.;
	double step;
	static constexpr double dTemp=0.00001;
	
	
	double _h;
	Vector<double,D> uni;
	
	
	void qmflow();
	
	
	
	double gradPressure_weight(int a, int b) {return _p[b].sigmaweight*_p[a].sigma*( _p[b].EOS.p()/(_p[b].sigma*_p[b].sigma) + _p[a].EOS.p()/(_p[a].sigma*_p[a].sigma) );}
//	double gradBulk_weight(int a, int b) {	return _p[b].sigmaweight*_p[a].sigma*( (_p[b].Bulk_N+_p[b].Bulk)/_p[b].sigma/_p[b].gamma/t + (_p[a].Bulk_N+_p[a].Bulk)/_p[a].sigma/_p[a].gamma/t );}
	
	Vector<double,D> gradKernel (Vector<double,D> a);
	
	double linint(double x, double x1, double x2, double y1, double y2);
	Vector<double,D> linint(double x, double x1, double x2, Vector<double,D> y1, Vector<double,D> y2);
	double linfac(double fac, double y1, double y2);
	Vector<double,D> linfac(double fac, Vector<double,D> y1, Vector<double,D> y2);
	Matrix<double,D+1,D+1> linfac(double fac, Matrix<double,D+1,D+1> y1, Matrix<double,D+1,D+1> y2);
public:
	int start,end,fnum;
	vector<double> sFO,Tfluc; //entropy at freezeout
	int qmf; //if==1 quantum mechanicanical corrections to the flow or added, if==0 no corrections are included

	double kernel (Vector<double,D> kin);

	double knorm,knorm2,kgrad,kgrad2;
	int number_part;
	double t0;
	static constexpr double tend=50.02;
	double t,dt;
	double factor;
	int frzc;
	double tau,taup, taupp;
	int rk2;
	double gd2;
	
	int gtyp;
	//double van,avgT,avgvt;

	int cfon;
	vector<int> list;
	int *lead;	
	int *link;
	int cf;
	int visc; //visc=0 for ideal, visc=1 for bulk, visc=2 for shear, visc=3 for bulk+shear
	double efcheck,sfcheck;
	int steps;
	//double *ep1,*ep2,*ec1,*ec2,*tsave,*pan,*san;
	
	int first;
	int average;
	int lowT;
	
	double *divTtemp,*gsub,*bulksub,*swsub,*shear33sub,*tlist;
	double avgetasig; // possibly not needed?
	Matrix<double,D+1,D+1> *shearsub;
	Vector<double,D> *divT,*rsub;
	Vector<double,D> *uout;
	double wfz,cs2;

	Vector<int,D> size;
	
	Vector<int,D> *dael;
	
	Particle<D> *_p;
	
	double dEz,E,Ez,E0,Etot,S,S0,Eloss,Esubb;
	
	double E1,E2;
	
	LinkList<D>();

	int etaconst;
	double bvf, svf,zwidth,sTc,zTc;
	
	void setup(double it0, int ntot,double h,Particle<D> *_p,double dtsave, int & numpart);
	void setupnext(int ntot, Particle<D> *_pin, int & numpart);
	
	string eos_s,eos_p;
	string *filenames;
	string ebe_folder;
	
	int fcount,cevent;
	string eost;
	
	~LinkList<D>();
	
	void initiate();
	
	void guess();
	
	int triToSum(Vector<int,D> dael, Vector<int,D> size);
	
	void destroy();
	
	void optimization(int a);
	void optimization2(int a);
	void svoptimization2(int a,double tin,int & count);
	
	void conservation_entropy();
	void conservation();
	void conservation_E();
	void conservation_Ez();
	void svconservation();
	void svconservation_E();
	void svconservation_Ez();
	int n() {return _n;}
//	void freezeout(int curfrz);
//	void interpolate(int curfrz);
	void svfreezeout(int curfrz);
	void svinterpolate(int curfrz);
	void freezeset();
	void optint(int a, double & T00,  double & Tx0);
	
	
	double eanal2(double t, double x, double y);
	double g(double t, double x, double y);
	double ux(double t, double x, double y);
	double uy(double t, double x, double y);
	void updateIC();
	void setv(string vtype);
	void etas_set();
	void sv_set();
	void sv_setb();
	void endEV();
	void setshear();
	void prints(); //possibly not needed?
	
};

template <int D>
LinkList<D>::LinkList()
{
	range=2; //number of boxes on the sides
	for(int i=0; i<D; i++) uni.x[i]=1.0;
	
	
}

template <int D>
void LinkList<D>::destroy()
{
delete [] lead;

}

template <int D>
void LinkList<D>::setup(double it0, int ntot,double h, Particle<D> *_pin,double dtsave,int & numpart)
{
	t0=it0;
	_h=h;
	_n=ntot;
	_p=_pin;
	knorm=10/7./PI/(_h*_h);
	knorm2=knorm*0.25;
	kgrad=-10/7./PI/pow(_h,3)*3/4.;
	kgrad2=10/7./PI/pow(_h,3)/_h;
	link=new int[_n];
	dael=new Vector<int,D>[_n];
	steps=100*(floor(tend-t0)+1);
	
	dt=dtsave;
	
	number_part=numpart;
	avgetasig=0.;	
}



template <int D>
void LinkList<D>::setupnext(int ntot, Particle<D > *_pin,int & numpart)
{
	_n=ntot;
	_p=_pin;
	number_part=numpart;
	link=new int[_n];
	dael=new Vector<int,D>[_n];
	avgetasig=0.;
}

template <int D>
void LinkList<D>::prints( )
{
	for (int i=0; i<_n; i++) {
    	_p[i].saves=_p[i].EOS.s()*_p[i].sigmaweight/_p[i].sigma/_p[i].gamma/t0;
    	}
}

template <int D>
void LinkList<D>::endEV()
{	
	delete [] link;
	delete [] dael;	
}

template <int D>
LinkList<D>::~LinkList()
{
}

template <int D>
int LinkList<D>::triToSum(Vector<int,D>  dael, Vector<int,D>  size)
{
	if (D==1)
	{
		return dael.x[0];
	}
	else if (D==2)
	{
		return dael.x[0] + dael.x[1]*size.x[0];
	}
	else if (D==3)
	{
		return dael.x[0] + dael.x[1]*size.x[0]+dael.x[2]*size.x[0]*size.x[1];
	}
	else 
	{
		cout << "Error: in LinkList triToSum, dimensions are not 1,2,3" << endl;
	}
	
	
}



//template <int D>
//void LinkList<D>::freezeout(int curfrz)
//{
//    if (frzc==0)
//    {
//    	taupp=t;	
//    	frzc=1;
//    	for (int i=0; i<_n; i++) {
//    	
//    	
//    	_p[i].frz2.r=_p[i].r;
//    	_p[i].frz2.u=_p[i].u;
//    	_p[i].frz2.sigma=_p[i].sigma;
//    	_p[i].frz2.T=_p[i].EOS.T();
//    	_p[i].frz2.theta=_p[i].div_u+_p[i].gamma/t;
//    	_p[i].frz2.gradP=_p[i].gradP;
//    	}
//    	
//    	
//    }
//    else if (frzc==1)
//    {    		
//   	taup=t;
//    	frzc=2;
//    	for (int i=0; i<_n; i++) {
//    	    	
//    	_p[i].frz1.r=_p[i].r;
//    	_p[i].frz1.u=_p[i].u;
//    	_p[i].frz1.sigma=_p[i].sigma;
//    	_p[i].frz1.T=_p[i].EOS.T();
//    	_p[i].frz1.theta=_p[i].div_u+_p[i].gamma/t;
//    	_p[i].frz1.gradP=_p[i].gradP;
//    	}
//
//	divTtemp=new double [curfrz];
//    	divT=new Vector<double,D> [curfrz];
//    	gsub=new double [curfrz];
//    	uout=new Vector<double,D> [curfrz];
//    	swsub=new double [curfrz];
//    	tlist=new double [curfrz];
//    	rsub=new Vector<double,D> [curfrz];
//
//	if (curfrz>0) interpolate(curfrz);
//    	else cf=0;
//    	 
//	
//    }
//    else
//    {
//    	
//    	
//    	tau=t;
//    	
//    	divTtemp=new double [curfrz];
//    	divT=new Vector<double,D> [curfrz];
//    	gsub=new double [curfrz];
//    	uout=new Vector<double,D> [curfrz];
//    	swsub=new double [curfrz];
//    	tlist=new double [curfrz];
//    	rsub=new Vector<double,D> [curfrz];
//    	
//    	
//    	
//    	
//    	if (curfrz>0) interpolate(curfrz);
//    	else cf=0;
//    	
//    	
//    	//sets up the variables for the next time step
//    	for (int i=0; i<_n; i++) {
//    	_p[i].frz2=_p[i].frz1;
//    	
//    	
//    	_p[i].frz1.r=_p[i].r;
//    	_p[i].frz1.u=_p[i].u;
//    	_p[i].frz1.sigma=_p[i].sigma;
//    	_p[i].frz1.T=_p[i].EOS.T();
//    	_p[i].frz1.theta=_p[i].div_u+_p[i].gamma/t;
//    	_p[i].frz1.gradP=_p[i].gradP;
//    	}
//    	taupp=taup;
//    	taup=tau;
//    }
//    cfon=0;
//}
//
//template <int D>
//void LinkList<D>::interpolate(int curfrz)
//{
//	
//	sFO.resize(curfrz,0);
//	Tfluc.resize(curfrz,0);
//	for (int j=0;j<curfrz;j++)
//	{
//		int i=list[j];
//		
//		int swit=0;
//		if (abs(_p[i].frz1.T-freezeoutT)<abs(_p[i].frz2.T-freezeoutT)) swit=1;
//		else swit=2;
//		
//		Vector<double,D> gradPsub;
//		double sigsub,thetasub;
//		
//		if (swit==1){
//			tlist[j]=taup;
//			rsub[j]=_p[i].frz1.r;
//			uout[j]=_p[i].frz1.u;		
//		
//			gradPsub=_p[i].frz1.gradP;			
//			sigsub=_p[i].frz1.sigma;
//			thetasub=_p[i].frz1.theta;
//			Tfluc[j]=_p[i].frz1.T;
//		}
//		else if (swit==2){
//			tlist[j]=taupp;
//			rsub[j]=_p[i].frz2.r;
//			uout[j]=_p[i].frz2.u;
//			bulksub[j]=_p[i].frz2.bulk;
//					
//			gradPsub=_p[i].frz2.gradP;			
//			sigsub=_p[i].frz2.sigma;
//			thetasub=_p[i].frz2.theta;
//			Tfluc[j]=_p[i].frz2.T;
//		}
//
//		
//		
//		
//		gsub[j]=sqrt( Norm2(uout[j]) + 1. );
//
//		sigsub/=gsub[j]*tlist[j];
//		
//		swsub[j]=_p[i].sigmaweight/sigsub;
//
//		
//		sFO[j]=_p[0].EOS.s_terms_T(Tfluc[j]);
//    		divT[j]=(1/sFO[j])*gradPsub;
//    		divTtemp[j]=-(1/(gsub[j]*sFO[j]))*(cs2*wfz*thetasub+inner(uout[j],gradPsub));
//                 
//                 
//              
//                
//		double insub=divTtemp[j]*divTtemp[j]-Norm2(divT[j]);
//		double norm=-sqrt(abs(insub));
//		divTtemp[j]/=norm;
//		divT[j]=(1/norm)*divT[j];
//
//		avgetasig+=sFO[j]/sigsub;
//		
//		sFO[j]*=pow(Tfluc[j]*0.1973,3);
//		Tfluc[j]*=0.1973;
//	
//	}
//	cf=curfrz;
//}


template <int D>
void LinkList<D>::svfreezeout(int curfrz)
{
    
    if (frzc==0)
    {
    	taupp=t;	
    	frzc=1;
    	for (int i=0; i<_n; i++) {
    	
    	
    	_p[i].frz2.r=_p[i].r;
    	_p[i].frz2.u=_p[i].u;
    	_p[i].frz2.sigma=_p[i].sigma;
    	_p[i].frz2.T=_p[i].EOS.T();
    	_p[i].frz2.bulk=_p[i].bigPI ;
    	_p[i].frz2.theta=_p[i].div_u+_p[i].gamma/t;
    	_p[i].frz2.gradP=_p[i].gradP;
    	_p[i].frz2.shear=_p[i].shv;
    	_p[i].frz2.shear33=_p[i].shv33;
    	_p[i].frz2.inside=_p[i].inside;
    	}
    	
    }
    else if (frzc==1)
    {    		
   	taup=t;
    	frzc=2;
    	for (int i=0; i<_n; i++) {
    	    	
    	_p[i].frz1.r=_p[i].r;
    	_p[i].frz1.u=_p[i].u;
    	_p[i].frz1.sigma=_p[i].sigma;
    	_p[i].frz1.T=_p[i].EOS.T();
    	_p[i].frz1.bulk=_p[i].bigPI ;
    	_p[i].frz1.theta=_p[i].div_u+_p[i].gamma/t;
    	_p[i].frz1.gradP=_p[i].gradP;
    	_p[i].frz1.shear=_p[i].shv;
    	_p[i].frz1.shear33=_p[i].shv33;
    	_p[i].frz1.inside=_p[i].inside;
    	}

	divTtemp=new double [curfrz];
    	divT=new Vector<double,D> [curfrz];
    	gsub=new double [curfrz];
    	uout=new Vector<double,D> [curfrz];
    	swsub=new double [curfrz];
    	bulksub=new double [curfrz];
    	shearsub=new Matrix<double,D+1,D+1> [curfrz];
    	shear33sub=new double [curfrz];
    	tlist=new double [curfrz];
    	rsub=new Vector<double,D> [curfrz];

	if (curfrz>0) 
    		svinterpolate(curfrz);
    	else 
    		cf=0;
    	 
    }
    else
    {
    	
    	for (int i=0; i<_n; i++) {
    	       if (_p[i].Freeze<4){
	    	if ((_p[i].btrack<=3)&&(_p[i].btrack>0)){
	    		_p[i].fback4=_p[i].fback2;
	    		_p[i].fback3=_p[i].fback; 
	    		_p[i].fback2=_p[i].frz2;
	    		_p[i].fback=_p[i].frz1;   	
	    	}
	    	else if (_p[i].btrack==0){
	    		if (_p[i].fback.gradP.x[0]!=0){
		    		_p[i].frz2=_p[i].fback2;
		    		_p[i].frz1=_p[i].fback;  
	    		}
	    		else{
	    			_p[i].frz2=_p[i].fback4;
		    		_p[i].frz1=_p[i].fback3; 
		    		cout << "back second"  << endl;
	    		} 
	    		
	    		
	    		curfrz++;
	    		list.push_back(i);
	    		_p[i].Freeze=4;	
	    		_p[i].btrack=-1;
	    	}
	    	}
    	}
    	
    	tau=t;
    	
    	divTtemp=new double [curfrz];
    	divT=new Vector<double,D> [curfrz];
    	gsub=new double [curfrz];
    	uout=new Vector<double,D> [curfrz];
    	swsub=new double [curfrz];
    	bulksub=new double [curfrz];
    	shearsub=new Matrix<double,D+1,D+1> [curfrz];
    	shear33sub=new double [curfrz];
    	tlist=new double [curfrz];
    	rsub=new Vector<double,D> [curfrz];
    	
    	
    	
    	if (curfrz>0) 
    		svinterpolate(curfrz);
    	else 
    		cf=0;
    	
    	
    	//sets up the variables for the next time step
    	for (int i=0; i<_n; i++) {
    	_p[i].frz2=_p[i].frz1;
    	
    	
    	_p[i].frz1.r=_p[i].r;
    	_p[i].frz1.u=_p[i].u;
    	_p[i].frz1.sigma=_p[i].sigma;
    	_p[i].frz1.T=_p[i].EOS.T();
    	_p[i].frz1.bulk=_p[i].bigPI ;
    	_p[i].frz1.theta=_p[i].div_u+_p[i].gamma/t;
    	_p[i].frz1.gradP=_p[i].gradP;
    	_p[i].frz1.shear=_p[i].shv;
    	_p[i].frz1.shear33=_p[i].shv33;
    	_p[i].frz1.inside=_p[i].inside;
    	}
    	taupp=taup;
    	taup=tau;
    }
    cfon=0;
}




template <int D>
void LinkList<D>::svinterpolate(int curfrz)
{
	
	sFO.resize(curfrz,0);
	Tfluc.resize(curfrz,0);
	for (int j=0;j<curfrz;j++)
	{
		
		
		int i=list[j];
		
		
		int swit=0;
		if (abs(_p[i].frz1.T-freezeoutT)<abs(_p[i].frz2.T-freezeoutT)) swit=1;
		else swit=2;
		
		double sigsub,thetasub,inside;
		Vector<double,D> gradPsub;
		if (swit==1){
			if (_p[i].btrack!=-1) tlist[j]=taup;
			else tlist[j]=taup-dt;
			rsub[j]=_p[i].frz1.r;
			uout[j]=_p[i].frz1.u;
			bulksub[j]=_p[i].frz1.bulk;
			shearsub[j]=_p[i].frz1.shear;
			shear33sub[j]=_p[i].frz1.shear33;
		
			gradPsub=_p[i].frz1.gradP;
			inside=_p[i].frz1.inside;
			sigsub=_p[i].frz1.sigma;
			thetasub=_p[i].frz1.theta;
			Tfluc[j]=_p[i].frz1.T;
		}
		else if (swit==2){
			if (_p[i].btrack!=-1) tlist[j]=taupp;
			else tlist[j]=taupp-dt;
			rsub[j]=_p[i].frz2.r;
			uout[j]=_p[i].frz2.u;
			bulksub[j]=_p[i].frz2.bulk;
			shearsub[j]=_p[i].frz2.shear;
			shear33sub[j]=_p[i].frz2.shear33;
		
			gradPsub=_p[i].frz2.gradP;
			inside=_p[i].frz2.inside;
			sigsub=_p[i].frz2.sigma;
			thetasub=_p[i].frz2.theta;
			Tfluc[j]=_p[i].frz2.T;
		}
		else {
		cout << "LinkList.h: Not at freeze-out temperature" << endl;
		
		}
		
		
		
		sFO[j]=_p[0].EOS.s_terms_T(Tfluc[j]);
		
		gsub[j]=sqrt( Norm2(uout[j]) + 1 );

		
		sigsub/=gsub[j]*tlist[j];
		swsub[j]=_p[i].sigmaweight/sigsub;
		
    		divT[j]=(1/sFO[j])*gradPsub;
    		divTtemp[j]=-(1/(gsub[j]*sFO[j]))*(cs2*(wfz+bulksub[j])*thetasub-cs2*inside+inner(uout[j],gradPsub));
                 
                
		double insub=divTtemp[j]*divTtemp[j]-Norm2(divT[j]);
		double norm=-sqrt(abs(insub));
		divTtemp[j]/=norm;
		divT[j]=(1/norm)*divT[j];
		
		avgetasig+=sFO[j]/sigsub;
			
		sFO[j]*=pow(Tfluc[j]*0.1973,3);
		Tfluc[j]*=0.1973;
	
	}
	cf=curfrz;
}


template <int D>
double LinkList<D>::linint(double x, double x1, double x2, double y1, double y2)
{
	return  y1+(x-x1)*(y2-y1)/(x2-x1);
}

template <int D>
Vector<double,D> LinkList<D>::linint(double x, double x1, double x2, Vector<double,D> y1, Vector<double,D> y2)
{
	return  y1+((x-x1)/(x2-x1))*(y2-y1);
}

template <int D>
double LinkList<D>::linfac(double fac, double y1, double y2)
{
	return  y1+fac*(y2-y1);
}

template <int D>
Vector<double,D> LinkList<D>::linfac(double fac, Vector<double,D> y1, Vector<double,D> y2)
{
	return  y1+fac*(y2-y1);
}

template <int D>
Matrix<double,D+1,D+1> LinkList<D>::linfac(double fac, Matrix<double,D+1,D+1> y1, Matrix<double,D+1,D+1> y2)
{
	return  y1+fac*(y2-y1);
}


template <int D>
void LinkList<D>::conservation()
{
    
    conservation_entropy();
    conservation_E();
    Etot=E+Ez;
   Eloss= (E0-Etot)/E0*100;
    rk2=0;
  
}


template <int D>
void LinkList<D>::svconservation()
{
    
//    conservation_entropy();
    svconservation_E();
    Etot=E+Ez;
   Eloss= (E0-Etot)/E0*100;
   rk2=0;
 //   cout << Eloss << "% of Energy loss at time t=" << t << endl;
//    cout << (S0-S)/S0 << "% of Entropy loss at time t=" << t << endl;
  
}

template <int D>
void LinkList<D>::conservation_entropy()
{
    
    S=0.;

    for (int i=0; i<_n; i++) {
    S+= _p[i].eta_sigma*_p[i].sigmaweight ;
    }
    
    if (first==1)
    {
    	
    	S0=S;
    }
}

template <int D>
void LinkList<D>::conservation_E()
{
    
    E=0.;
	
    for (int i=0; i<_n; i++) {
    E+= (_p[i].EOS.w()* _p[i].gamma* _p[i].gamma-_p[i].EOS.p())/_p[i].sigma*_p[i].sigmaweight*t;
    }
    
    
    
    if (first==1)
    {
    	first=0;    	
    	E0=E;
    }
    

}

template <int D>
void LinkList<D>::freezeset()
{
	cs2=_p[0].EOS.cs2out(freezeoutT);
	wfz=_p[0].EOS.wfz(freezeoutT);
}

template <int D>
void LinkList<D>::sv_set()
{
        
	for (int i=0; i<_n; i++) {
	double gg=_p[i].gamcalc();
	_p[i].g2=gg*gg;
    	_p[i].shv33=0;
    	}

}

template <int D>
void LinkList<D>::sv_setb()
{
        
	for (int i=0; i<_n; i++) {
	double gg=_p[i].gamcalc();
	_p[i].g2=gg*gg;
    	}

}

template <int D>
void LinkList<D>::svconservation_E()
{
    
    E=0.;
    for (int i=0; i<_n; i++) {
    E+= (_p[i].C* _p[i].g2-_p[i].EOS.p()-_p[i].bigPI+_p[i].shv.x[0][0])/_p[i].sigma*_p[i].sigmaweight*t;
    }
    
    if (first==1)
    {
    	
    	first=0;
    	E0=E;
    }
    

}


template <int D>
void LinkList<D>::conservation_Ez()
{
    
    dEz=0.;

    
    for (int i=0; i<_n; i++) {
    dEz+= _p[i].EOS.p()/_p[i].sigma*_p[i].sigmaweight;
    }

}


template <int D>
void LinkList<D>::svconservation_Ez()
{
    
    dEz=0.;

    double t2=t*t;
    for (int i=0; i<_n; i++) {
    dEz+=( _p[i].EOS.p()+_p[i].bigPI+_p[i].shv33*t2)/_p[i].sigma*_p[i].sigmaweight;
    }

}

template <int D>
double LinkList<D>::kernel (Vector<double,D> a)
{
	double r=Norm(a);
	double q=r/_h;
	double qq=q*q;

	
	if(q>=2) return 0;
	
	if(q>=1) return knorm2*pow((2-q),3);
	
	return knorm*(1 - 1.5*qq + 0.75*q*qq);
}

template <int D>
Vector<double,D> LinkList<D>::gradKernel (Vector<double,D> a)
{
	Vector<double,D> tsubb;
	
	
	double r=Norm(a);
	double q=r/_h;
	
	

	
	if(q>2) 
		return tsubb;
	if(q>1)	
		return kgrad/r*pow((2-q),2)*a;
	
	return kgrad2*(-3+9*q/4.)*a;

}

template <int D>
void LinkList<D>::setshear()
{
    
	double t2=t*t;
    for (int i=0; i<_n; i++) {
    _p[i].sets(t2);         
    }
    
    
    

}

template <int D>
void LinkList<D>::initiate() {
	
	
	
// check what happens with particle separates by itself?  Where in fortran code?	


	//find system boundaries

		max=min=_p[0].r;
		
		
		

		for(int i=1; i<_n; i++)	
		{
		for(int j=0;j<D;j++) {

			if(_p[i].r.x[j]>max.x[j]) max.x[j]=_p[i].r.x[j];

			if(_p[i].r.x[j]<min.x[j]) min.x[j]=_p[i].r.x[j];

		}
		}
		
		
		
	
	//evaluate system size

		//2*range puts extra boxes on sides of grid
		double sub=1./_h;	
		 	
		size=sub*(max-min)+(2.*range+1.)*uni;


		//Size is the volume
		Size=1;
		
		

	
	Vector<double,D> dsub;
	
	// finds total volume of the system
        for(int i=0; i<D; i++) Size*=size.x[i];

	    
	//dael: relates every particle with its linklist cube 
	// also convert particle position to an integer

	for(int j=0; j<_n; j++) {

			
			
			dael[j]=sub*(_p[j].r-min)+(1.*range)*uni;
		
			
	}
	
	//lead: relates every linklist cube with one of the particles (leader) in it
	//link: links the leader particle of one cube with the others of the same cube
	// if only one particle in cube then it is the lead
	
		

	lead = new int[Size];

	for(int j=0; j<Size; j++) lead[j]=-1;
				

	for(int k=_n-1; k>=0; k--) {
		
			
			int tt=triToSum(dael[k],size); // need to understand still... seems like it needs another coordinate
		
			link[k]=lead[tt];

			lead[tt]=k;
			
			
			

	}
	
	return;
}

template <int D>
void LinkList<D>::optimization(int a)
{

    _p[a].sigma = 0;
    Vector<int,D> i;
	for(i.x[0]=-2; i.x[0]<=2; i.x[0]++)
        {
        for(i.x[1]=-2; i.x[1]<=2; i.x[1]++)
        {

        	

                 int b=lead[triToSum(dael[a]+i, size)];
                 while(b!=-1 )
                 {
                 	double kern=kernel(_p[a].r-_p[b].r);
                    	_p[a].sigma += _p[b].sigmaweight*kern;
			
                    	
                   	b=link[b];

			

                   	
                   	
                 }
        }
	}

	
	_p[a].eta =  _p[a].sigma*_p[a].eta_sigma; 
	
	
}



template <int D>
void LinkList<D>::optint(int a, double & ux0,  double & uy0)
{

    _p[a].sigma = 0;
    Vector<int,D> i;
  uy0=0;
  ux0=0;
	for(i.x[0]=-2; i.x[0]<=2; i.x[0]++)
        {
        for(i.x[1]=-2; i.x[1]<=2; i.x[1]++)
        {

        	

                 int b=lead[triToSum(dael[a]+i, size)];
                 while(b!=-1 )
                 {
                 	double kern=kernel(_p[a].r-_p[b].r);
                 	
                    	_p[a].sigma += _p[b].sigmaweight*kern;
                    	
                    	
			ux0+=_p[b].u.x[0]*gd2*kern;
			uy0+=_p[b].u.x[1]*gd2*kern;
                    	
                   	b=link[b];

			

                   	
                   	
                 }
        }
	}

	

	_p[a].eta =  _p[a].sigma*_p[a].eta_sigma; 
	
	
}



template <int D>
void LinkList<D>::optimization2(int a)
{

    _p[a].gradP=0;
    _p[a].dsigma_dt = 0;
    _p[a].gradsig=0;

    Vector<int,D> i;
    
		
	
	
	
	for(i.x[0]=-2; i.x[0]<=2; i.x[0]++)
        {

        for(i.x[1]=-2; i.x[1]<=2; i.x[1]++)
        {

                 int b=lead[triToSum(dael[a]+i, size)];


                    while(b!=-1 )
                    {

                    Vector<double,D> siggradK=_p[b].sigmaweight*gradKernel(_p[a].r-_p[b].r);

			
                    _p[a].dsigma_dt+= inner(siggradK,_p[a].v-_p[b].v);
                    _p[a].gradP += _p[a].sigma*( _p[b].EOS.p()/(_p[b].sigma*_p[b].sigma) + _p[a].EOS.p()/(_p[a].sigma*_p[a].sigma) )*siggradK;
                    _p[a].gradsig +=siggradK;
                    

                    b=link[b];
                    }
        }
	}
	

}



template <int D>
void LinkList<D>::svoptimization2(int a,double tin,int & count)
{

    _p[a].gradP=0.;
    _p[a].gradBulk = 0.;
    _p[a].gradV = 0.;
    _p[a].gradshear=0.;
    _p[a].divshear=0.;

    Vector<int,D> i;
		
	if (_p[a].btrack!=-1) _p[a].btrack=0;
	
	double rdis=0;

	for(i.x[0]=-2; i.x[0]<=2; i.x[0]++)
        {
        for(i.x[1]=-2; i.x[1]<=2; i.x[1]++)
        {

                 int b=lead[triToSum(dael[a]+i, size)];
                 
                
                 


                    while(b!=-1 )
                    {

                    
                    Vector<double,D> gradK=gradKernel(_p[a].r-_p[b].r);
                    Vector<double,D> va=rowp1(0,_p[a].shv);
                    Vector<double,D> vb=rowp1(0,_p[b].shv);
                    Matrix<double,D,D> vminia,vminib;
                    mini(vminia,_p[a].shv);
                    mini(vminib,_p[b].shv);
                    double sigsqra=1/(_p[a].sigma*_p[a].sigma);
                    double sigsqrb=1/(_p[b].sigma*_p[b].sigma);
                    Vector<double,D> sigsigK=_p[b].sigmaweight*_p[a].sigma*gradK;

                    _p[a].gradP +=( sigsqrb*_p[b].EOS.p()+ sigsqra*_p[a].EOS.p() )*sigsigK;
                    
                    if (((Norm(_p[a].r-_p[b].r)/_h)<=2)&&(a!=b)) {
                     	if (_p[a].btrack!=-1) _p[a].btrack++;
                     	if (_p[a].btrack==1) rdis=Norm(_p[a].r-_p[b].r)/_h;
                     }
                    
                    _p[a].gradBulk += ( _p[b].Bulk/_p[b].sigma/_p[b].gamma + _p[a].Bulk/_p[a].sigma/_p[a].gamma)/tin*sigsigK;
                    _p[a].gradV+=(_p[b].sigmaweight/_p[a].sigma*( _p[b].v -  _p[a].v ))*gradK;
                    
                   
                    
                    _p[a].gradshear+=(inner(sigsigK,_p[a].v))*(sigsqrb*vb+sigsqra*va);
                    _p[a].divshear+=(sigsqrb*(sigsigK*transpose(vminib))+sigsqra*(sigsigK*transpose(vminia)));

                    b=link[b];
                    }
        }
	}
	
	if ((_p[a].btrack==1)&&((_p[a].EOS.T()*197.3)>=150)) {	
    		_p[a].frz2.t=tin;
    	}
    	else if ((_p[a].btrack==0)&&((_p[a].EOS.T()*197.3)>=150)&&(_p[a].Freeze<4)){
    	cout <<"Missed " << a << " "<< tin << "  " << _p[a].EOS.T()*197.3 << " " << rdis << " " << cfon <<  endl;	
    	}
	

}


template <int D>
double LinkList<D>::eanal2(double t, double x, double y)
{
	double xperp=x*x+y*y;
	double t2=t*t;
	
	return e0/pow(t,4/3.)*pow(2*q,8/3.)/pow((1+2*q*q*(t2+xperp)+pow(q,4)*pow((t2-xperp),2)),4/3.);
}

template <int D>
double LinkList<D>::g(double t, double x, double y)
{
	double xperp=x*x+y*y;
	double q2=q*q;
	return atanh(2*q2*t*sqrt(xperp)/(1+q2*t*t+q2*xperp));
}

template <int D>
double LinkList<D>::ux(double t, double x, double y)
{
	return sinh(g(t,x,y))*x/sqrt(x*x+y*y);
}

template <int D>
double LinkList<D>::uy(double t, double x, double y)
{
	return sinh(g(t,x,y))*y/sqrt(x*x+y*y);
}


template <int D>
void LinkList<D>::updateIC()
{
	
	
	for (int i=0;i<_n;i++)
	{	
	        if (gtyp!=5) _p[i].s_an=_p[i].EOS.s_out(_p[i].e_sub);
		
                
		_p[i].EOS.update_s(_p[i].s_an);
		if (gtyp==5) _p[i].e_sub=_p[i].EOS.e();
                _p[i].gamma=_p[i].gamcalc();
		_p[i].sigmaweight  *= _p[i].s_an*_p[i].gamma*t0;
		
	}
	
	if (gtyp!=3) guess();
	
}





template <int D>
void LinkList<D>::guess()
{

	initiate();

        for (int i=0;i<_n;i++)
        {

		optimization(i);
	}
	
        double tmax=0;
        Vector<double,D> rmax;
	for (int j=0;j<D;j++)
        {
                rmax.x[j]=0;
	}
//      int imax;                                                                                                                                     


	int count1=0;

	for (int i=0;i<_n;i++)
        {

                _p[i].s_sub=_p[i].sigma/_p[i].gamma/t0;

                 _p[i].EOS.update_s(_p[i].s_sub);
		
		_p[i].sigsub=0;

		if (_p[i].EOS.T()>tmax)
		{

                        tmax=_p[i].EOS.T();

		        rmax=_p[i].r;

                }

	        _p[i].frzcheck(t0,count1,_n);

        }

	//cout << "largest temperature" << endl;                                                                                                      
	cout << tmax*197.3 << " " << rmax <<endl;

	//Ss=0.;                                                                                                                                      

        if (qmf==1) qmflow();
}




template <int D>
void LinkList<D>::qmflow()
{
	double pi2=2*M_PI;
	
	int run=0;
	for (int i=0;i<_n;i++)
	{
		double E=_p[i].EOS.e()/_p[i].sigma*_p[i].gamma*t0;
		double sub=-1./(2.*pow((E*_h),2));
		
		int on=0;
		double v;
		
		while (on==0){
		v=rand01(run);
		if (rand01(run)<exp(v*v*sub)) on=1;				
		}
		double phi=rand0(pi2,run);
		//cout << phi <<endl;
		_p[i].v.x[0]=_p[i].u.x[0]/_p[i].gamma+v*cos(phi);
		_p[i].v.x[1]=_p[i].u.x[1]/_p[i].gamma+v*sin(phi);
		
		double gamnew=1/sqrt(1-Norm2(_p[i].v));
		_p[i].sigmaweight  *= gamnew/_p[i].gamma;
		_p[i].gamma=gamnew;
		_p[i].u.x[0]=_p[i].gamma*_p[i].v.x[0];
        	_p[i].u.x[1]=_p[i].gamma*_p[i].v.x[1]; 
	
	}
	
	initiate();

	
	for (int i=0;i<_n;i++)
	{
		
		optimization(i);
	}

	for (int i=0;i<_n;i++)
	{
		
		_p[i].s_sub=_p[i].sigma/_p[i].gamma/t0;
		
		_p[i].EOS.update_s(_p[i].s_sub);
	}


}


template <int D>
void LinkList<D>::setv(string vtype)
{
   	string bulkshear ("bulk+shear");
   	string shearbulk ("shear+bulk");
   	string ideal ("ideal");
   	
   	if (vtype==ideal)
   	{
   		visc=0;
   	}
   	else if ((vtype==bulkshear) || (vtype==shearbulk))
   	{
   		visc=3;
   	}

}

#endif
