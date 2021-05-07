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
#include "bbmg.h"

void svSimulation(double dt,LinkList<2> &linklist)
{
    cout << "Ready to start hydrodynamics\n";

    linklist.frzc=0;
    linklist.cf=0;


    Output<2> out(linklist);

    BBMG<2> bbmg(linklist);
    bbmg.initial(linklist);
    cout << "started bbmg" << endl;

    linklist.t=linklist.t0;

    if (linklist.qmf==1||linklist.qmf==3) {
        out.sveprofile(linklist);
        cout << "printed first timestep" << endl;
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " S=" << linklist.S << endl;
        if (linklist.qmf==1) exit(0);
    }
    else if(linklist.qmf==4) {
        out.eccout(linklist);
        cout << "eccentricity printed" << endl;
        exit(0);
    }

    linklist.Ez=0;
    while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {
        linklist.cfon=1;


        svrungeKutta2<2>(dt,&shear<2>,linklist);
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " " <<  linklist.Eloss << " " << linklist.S <<  endl;
        out.sveprofile(linklist);

        if (linklist.cf>0) out.svFOprint(linklist);

        if (linklist.qmf==3) {
            double tsub=linklist.t-floor(linklist.t);
            // if you add more points to print, must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
            if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                linklist.conservation_entropy();
                cout << "t=" << linklist.t << " S=" << linklist.S << endl;  // outputs time step
                out.sveprofile(linklist);   // energy density profile
                cout << "eloss= " << linklist.t << " " <<  linklist.Eloss << endl;
                //out.conservation(linklist); // conservation of energy

            }

            else if ((tsub<(0.5+dt*0.5))&&(tsub>=(0.5-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                cout << "t=" <<  linklist.t <<endl;  // outputs time step
                out.sveprofile(linklist);   // energy density profile
            }

        }

    }

    linklist.endEV();
}

void BSQSimulation(double dt,LinkList<2> &linklist)
{
    cout << "Ready to start hydrodynamics\n";

    linklist.frzc=0;
    linklist.cf=0;

    Output<2> out(linklist);

    BBMG<2> bbmg(linklist);
    bbmg.initial(linklist);
    cout << "started bbmg" << endl;

    linklist.t=linklist.t0;

    if (linklist.qmf==1||linklist.qmf==3) {
        out.bsqsveprofile(linklist);
        cout << "printed first timestep" << endl;
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " S=" << linklist.S << endl;
        if (linklist.qmf==1) exit(0);
    }
    else if(linklist.qmf==4) {
        out.eccout(linklist);
        cout << "eccentricity printed" << endl;
        exit(0);
    }

    //out.sveprofile(linklist);
    linklist.Ez=0;
    while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {
        linklist.cfon=1;


        bsqrungeKutta2<2>(dt,&BSQshear<2>,linklist);
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " " <<  linklist.Eloss << " " << linklist.S <<  endl;
        out.bsqsveprofile(linklist);

        if (linklist.cf>0) out.bsqsvFOprint(linklist);

        if (linklist.qmf==3) {
            double tsub=linklist.t-floor(linklist.t);
            // if you add more points to print, must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
            if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                linklist.conservation_entropy();
                cout << "t=" << linklist.t << " S=" << linklist.S << endl;  // outputs time step
                out.bsqsveprofile(linklist);   // energy density profile
                cout << "eloss= " << linklist.t << " " <<  linklist.Eloss << endl;
                //out.conservation(linklist); // conservation of energy

            }

            else if ((tsub<(0.5+dt*0.5))&&(tsub>=(0.5-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                cout << "t=" <<  linklist.t <<endl;  // outputs time step
                out.bsqsveprofile(linklist);   // energy density profile
            }
        }
    }

    linklist.endEV();
}



template <int D>
void shear(LinkList<D>  &linklist)  // shear+bulk Equations of motion, only set up completely for 2+1 at the moment
{

    linklist.setshear();
    linklist.initiate();

    for(int i=0; i<linklist.n(); i++)
    {
        linklist.optimization(i);


        if ((linklist._p[i].eta<0)||isnan(linklist._p[i].eta))
        {
            cout << i << " neg entropy "
					<< linklist._p[i].EOS.T()*197.3 << " "
					<< linklist._p[i].eta << endl;
            linklist._p[i].eta=0;
        }
    }


    int curfrz=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //  Computes gamma and velocity
        linklist._p[i].calc(    linklist.t );
        linklist._p[i].setvisc( linklist.etaconst, linklist.bvf, linklist.svf,
                                linklist.zTc, linklist.sTc, linklist.zwidth,
                                linklist.visc );
        if (linklist.cfon==1) linklist._p[i].frzcheck( linklist.t, curfrz, linklist.n() );

    }
    if (linklist.cfon==1)
    {
        linklist.number_part+=curfrz;
        linklist.list.resize(curfrz);
    }

    int m=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
        linklist.svoptimization2(i, linklist.t, curfrz);
        linklist._p[i].dsigma_dt
			= -linklist._p[i].sigma * ( 
					linklist._p[i].gradV.x[0][0]
					+ linklist._p[i].gradV.x[1][1] );

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
    for(int i=0; i<linklist.n(); i++)
    {
        double gamt=1./linklist._p[i].gamma/linklist._p[i].stauRelax;
        double pre=linklist._p[i].eta_o_tau/2./linklist._p[i].gamma;
        double p1=gamt-4./3./linklist._p[i].sigma*linklist._p[i].dsigma_dt+1./linklist.t/3.;
        Vector<double,D> minshv=rowp1(0,linklist._p[i].shv);
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
        linklist._p[i].du_dt.x[0]=(F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1]);
        linklist._p[i].du_dt.x[1]=(F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1]);

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

template <int D>
void BSQshear(LinkList<D>  &linklist)  // shear+bulk Equations of motion, only set up completely for 2+1 at the moment
{

    linklist.setshear();
    linklist.initiate();

    for(int i=0; i<linklist.n(); i++)
    {
        int curfrz=0;//added by Christopher Plumberg to get compilation
        linklist.bsqsvoptimization(i);    // NOT bsqsvoptimization2!!! fix arguments accordingly!!!

        if ((linklist._p[i].eta<0)||isnan(linklist._p[i].eta))
        {
            cout << i <<  " neg entropy " <<  linklist._p[i].EOS.T()*197.3   << " " << linklist._p[i].eta << endl;

            linklist._p[i].eta=0;
        }

    }


    int curfrz=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //  Computes gamma and velocity

        linklist._p[i].calcbsq(linklist.t); //resets EOS!!
        linklist._p[i].setvisc(linklist.etaconst,linklist.bvf,linklist.svf,linklist.zTc,linklist.sTc,linklist.zwidth,linklist.visc);
        if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n());

    }
    if (linklist.cfon==1)
    {
        linklist.number_part+=curfrz;
        linklist.list.resize(curfrz);
    }

    int m=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
        linklist.bsqsvoptimization2(i,linklist.t,curfrz);

        linklist._p[i].dsigma_dt = -linklist._p[i].sigma*(linklist._p[i].gradV.x[0][0]+linklist._p[i].gradV.x[1][1]) ;

        linklist._p[i].bsqsvsigset(linklist.t,i);
        if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1))
        {
            linklist.list[m]=i;
            linklist._p[i].Freeze=4;
            ++m;
        }

    }

    if (linklist.rk2==1) linklist.bsqsvconservation();
    linklist.bsqsvconservation_Ez();


    //calculate matrix elements
    for(int i=0; i<linklist.n(); i++)
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


        linklist._p[i].inside=linklist.t*(inner((-minshv+linklist._p[i].shv.x[0][0]*linklist._p[i].v),linklist._p[i].du_dt)- con2(sub,linklist._p[i].gradU)    -      linklist._p[i].gamma*linklist.t*linklist._p[i].shv33);
        linklist._p[i].detasigma_dt =1./linklist._p[i].sigma/linklist._p[i].EOS.T()*( -linklist._p[i].bigPI*linklist._p[i].bigtheta+linklist._p[i].inside);


        linklist._p[i].dBulk_dt = (-linklist._p[i].zeta/linklist._p[i].sigma*linklist._p[i].bigtheta - linklist._p[i].Bulk/linklist._p[i].gamma )/linklist._p[i].tauRelax;

        Matrix <double,D,D> ududt=linklist._p[i].u*linklist._p[i].du_dt;


        linklist._p[i].dshv_dt= -gamt*(linklist._p[i].pimin+linklist._p[i].setas*0.5*partU)-0.5*linklist._p[i].eta_o_tau*(ududt+transpose(ududt))+linklist._p[i].dpidtsub()-vduk*(ulpi+transpose(ulpi)+(1/linklist._p[i].gamma)*Ipi)+linklist._p[i].sigl*Ipi;

        linklist._p[i].drhoB_dt=-linklist._p[i].rhoB*linklist._p[i].sigma*linklist._p[i].bigtheta;
        linklist._p[i].drhoS_dt=-linklist._p[i].rhoS*linklist._p[i].sigma*linklist._p[i].bigtheta;
        linklist._p[i].drhoQ_dt=-linklist._p[i].rhoQ*linklist._p[i].sigma*linklist._p[i].bigtheta;

    }


    if (linklist.cfon==1) linklist.bsqsvfreezeout(curfrz);


    linklist.destroy();
}
