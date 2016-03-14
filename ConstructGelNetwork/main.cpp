//
//  main.cpp
//  ConstructGelNetwork
//
//  Created by Ryohei Seto on 2012/08/17.
//
//
#include <iostream>
#include "threedim.h"
#include "twodim.h"
using namespace std;
double r;
int n;
double Lx; // box size
double Lz; // box size
double Lz_pd; // for periodic boundary: Lz_pd = Lz -2;

int main (int argc, const char * argv[])
{
    if (argc == 1){
		cerr << "1:dimension 2:filling 3:Lx 4:Lz 5:rank 6:rsd 7.phi " << endl;
        cerr << " (filling = 1 or 2 means sediment or random walk)" << endl;
		exit(1);
	}
    int dimension = atoi(argv[1]);
    int filling = atoi(argv[2]);
    /* box size L*Lz
     */
    Lx = atof(argv[3]);
    Lz = atof(argv[4]);
    /* periodic boundary condition:
     * 0 < x < L
     * 0 < z < Lz - 2
     */
    //Lz_pd = Lz - 2;
	Lz_pd = Lz;
    int rank = atoi(argv[5]);
	int rsd = atoi(argv[6]);
	double vol_frac = atof( argv[7]);
    
	//	int total_number = 10000;
    //int m = (int)(log(total_number)/log(2))+1;
    switch (dimension) {
        case 2:
            //make2dNetwork(double a, double lo);
            buildSpaceFillingNetwork2D(filling, rank, rsd, vol_frac);
            break;
        case 3:
			//make3dNetwork(rank, m, rsd);
            break;
        default:
            cerr << "dimension = 2 or 3" << endl;
            return 1;
    }
}
