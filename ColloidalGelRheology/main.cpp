//
//  main.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//
/*
 * -DTWODIMENSION=1
 * TWODIMENSION=1
 * Apple LLVM compiler 4.0
 * Preprocessor Macros Not Used in Precompiled Headers
 */
#include <iostream>
#include <fstream>
#include <algorithm>
#include "output.h"
using namespace std;
void cerr_command_usage();
void rheologyTest(int argc, char *const argv[], System &sy);
void bendingLinearChain(int argc, char *const argv[], System &sy);
//void testParameters(int argc, char *const argv[], System &sy);
void rod_bending(int argc, char *const argv[], System &sy);
int main (int argc, char *const argv[]) {
	System sy;
	sy.simulation = argv[1];
	rheologyTest(argc, argv, sy);
	return 0;
}

void cerr_command_usage()
{
	cerr << "usage of ColloidalGelRheology" << endl;
	cerr << "arg: c1/c2/s/bt" << endl;
	cerr << "  c1: uniaxial compression" << endl;
	cerr << "  c2: biaxial compression" << endl;
	cerr << "  s:  shear" << endl;
	cerr << "  bt: bending test" << endl;
	cerr << "parameter_file ";
	cerr << "initial_file ";
	cerr << "version ";
	cerr << endl;
}

void rheologyTest(int argc, char *const argv[], System &sy)
{
	/* arguments */
	if (argc < 5) {
		cerr_command_usage();
		exit(1);
	}
	sy.setParameterFile(argv[2]);
	sy.setVersion(argv[4]);
	/* read files */
	sy.readParameterFile();
	sy.readBondParameter();
	sy.setInitClusterFile(argv[3]);
	sy.importPositions();
#ifndef TWODIMENSION
	cerr <<  "3D simulation" << endl;
#else
	cerr <<  "2D simulation" << endl;
#endif
	sy.strainControlSimulation();
	return;
}

void bendingLinearChain(int argc, char *const argv[], System &sy)
{
	sy.setParameterFile(argv[2]);
	/* read files */
	sy.readParameterFile();
	sy.readBondParameter();
	sy.lx = 50;
	sy.ly = 0;
	sy.lz_init = 50;
	sy.lz = sy.lz_init;
	sy.setLinearChain(11);
	sy.setExternalforce(atof(argv[3]));
	//	sy.importPositions();
#ifndef TWODIMENSION
	cerr <<  "3D simulation" << endl;
#else
	cerr <<  "2D simulation" << endl;
#endif
	//	sy.strainControlSimulation();
	sy.bendingSimulation();
	return;
}
