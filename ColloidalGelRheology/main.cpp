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
//void testParameters(int argc, char *const argv[], System &sy);
void rod_bending(int argc, char *const argv[], System &sy);

int main (int argc, char *const argv[])
{
	System sy;
	sy.simulation = argv[1][0];
	if (sy.simulation == 't'){
		//testParameters(argc, argv, sy);
	} else {
		rheologyTest(argc, argv, sy);
	}
	return 0;
}

void cerr_command_usage()
{
	cerr << "usage of ColloidalGelRheology" << endl;
	cerr << "arg: c/s ";
	cerr << "parameter_file ";
	cerr << "initial_file ";
	cerr << "version ";
	cerr << endl;
}

void rheologyTest(int argc, char *const argv[], System &sy){
	/* arguments */
	if (argc < 5){
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

//void testParameters(int argc, char *const argv[], System &sy){
//	sy.setParameterFile(argv[2]);
//	/* read files */
//	sy.readParameterFile();
//	sy.readBondParameter();
//	sy.lx = 50.;
//	sy.ly = 50.;
//	sy.lz = 50.;
//	sy.lx0 = 0.5*sy.lx;
//	sy.ly0 = 0.5*sy.ly;
//	sy.lz0 = 0.5*sy.lz;
//	
//	sy.initialprocess = true;
//	
//	sy.particle.push_back(new Particle(0, vec3d(0.,0.,0.), 0, sy));
//	sy.particle.push_back(new Particle(1, vec3d(2.,0.,0.), 0, sy));
//	//particle.push_back(new Particle(2, vec3d(4,0,0), sy));
//	sy.n_particle = 2;
//	for (int i=0; i< sy.n_particle ; i++){
//		sy.particle[i]->p.x += sy.lx0;
//		sy.particle[i]->p.y += sy.ly0;
//		sy.particle[i]->p.z += sy.lz0;
//		sy.particle[i]->cerr();
//	}
//	sy.dt = sy.dt_max;
//	// sy.setTimeStep(0.00001);
//	/* prepare calcuration */
//	sy.initGrid();
//	sy.setBondGenerationDistance(2.0);
//	sy.setSimulationViscosity();
//	sy.preparationOutput();
//	
//	sy.initDEM();
//	sy.makeInitialBond(2.0001);
//	sy.outputYaplot();
//	/******************************/
//	
//	vec3d rotate_axis;
//	double angle ;
//	///* tortion */
//	//angle = M_PI/6;
//	//rotate_axis.set(1,0,0);
//	//sy.particle[1]->setRotate(rotate_axis,angle);
//	///////////////////////
//	
//	/* bending */
//	angle = M_PI/6;
//	rotate_axis.set(0,-1.,0);
//	sy.particle[1]->setRotate(rotate_axis,angle);
//	///////////////////////
//	/* sliding */
//	//particle[1]->p.z += 0.5;
//	/* elongation */
//	// sy.particle[1]->p.x += -0.1;
//	sy.particle[1]->p.z += 0.5;
//	sy.outputYaplot();
//	sy.initialprocess = false;
//	
//	foreach( vector<Bond *>, sy.bond, it_bond){
//		(*it_bond)->addContactForce();
//	}
//	
//	sy.particle[0]->resetForce();
//	sy.particle[1]->resetForce();
//
//	//particle[2]->p.set(4.0,0,0.5);
//	cerr << "simulation start" << endl;
//	
//	for (int t = 0; t < 10000; t++){
//		//sy.outputConfiguration('n');
//		sy.outputYaplot();
//		sy.simuAdjustment();
//		foreach( vector<Bond *>, sy.bond, it_bond){
//			(*it_bond)->addContactForce();
//		}
//		//   sy.particle[0]->resetForce();
//		//////////////////////////////
//		// calculate
//		// x*_{j+1} and v*_{j+1}
//		foreach( vector<Particle *>, sy.particle, it_particle){
//			(*it_particle)->move_Euler();
//		}
//		cout << t << ' ' << sy.particle[1]->p.z << endl;
//	}
//	
//	/* particle_active and bond_active are effective.
//	 * The others are not effective;
//	 */
//	//sy.outputParameters();
//	/*****************************/
//	//compression_process(sy);
//	//relaxation_process(sy);
//	//shear_process(sy);
//	/*****************************/
//	
//}
