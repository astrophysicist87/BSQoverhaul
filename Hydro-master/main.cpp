#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <string.h>

using namespace std;


#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "enteric.h"
#include "hydrosim.h"
#include "LinkList.h"
#include "eostables.h"
#include "matrix.h"


char ifolder []="inputfiles/";
string ofolder ("./outputfiles/");
double freezeoutT;
double zconst=1/(4*PI);
eostableshigh *ETH;
eostableslow *ETL;
int nETH,nETL;
Matrix <double,2,2> Imat;
int hcor;


int main (int argc, char *argv[])
{
	_inputIC ics;

	ics.man   = "settings.inp";
	ics.rnum  = "0";
	ics.start = 0;
	ics.end   = 0;
	ics.on    = 1;
	
	Imat.identity();

	LinkList<2> linklist;
	manualenter<2>(ics,linklist);

	switch ( linklist.visc )
	{
		case 0:	// ideal
			Simulation(   ics.dt, linklist );
			break;
		//case 1:	// bulk
		//	vSimulation(  ics.dt, linklist );
		//	break;
		case 3: // shear+bulk
			svSimulation( ics.dt, linklist );
			break;
		default:
			std::cerr << "Error: visc = " << linklist.visc
				  << " is not a valid option!" << std::endl;
			break;
	}

	string table ("table");
	if ( linklist.eost == table )
	{
		delete [] ETH;
		delete [] ETL;
	}

	return 0;
}
