#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#include "enteric.h"
#include "eostables.h"
#include "hydrosim.h"
#include "matrix.h"
#include "particle.h"
#include "tables.h"
#include "vector.h"
#include "LinkList.h"

//char ifolder []="inputfiles/";
//string ofolder ("./outputfiles/");
/*double freezeoutT;
double freezeoutB;
double freezeoutS;
double freezeoutQ;
double zconst=1/(4*PI);
eostableshigh *ETH;
eostableslow *ETL;
int nETH,nETL;
int hcor;
Matrix <double,2,2> Imat;*/


int main (int argc, char *argv[])
{
	// read in command-line arguments
    if (argc > 1)
    {
        ics.man = argv[1];
    }
    else
    {
		std::cerr << "Please specify input file." << std::endl;
    }

    //Imat.identity();
	
	LinkList<2> linklist;
	manualenter<2>(ics,linklist);

	switch ( linklist.visc )
	{
		case 3: // shear+bulk
			svSimulation(ics.dt,linklist);
			break;
		case 4: // shear+bulk+BSQ
			BSQSimulation(ics.dt,linklist);
			break;
		default;
			std::cerr << "Error: visc = " << linklist.visc
				  << " is not a valid option!" << std::endl;
			break;
	}

    return 0;
}
