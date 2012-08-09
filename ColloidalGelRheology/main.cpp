//
//  main.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//
/*
 * -DTWODIMENSION=1
 */
#include <iostream>
#include <fstream>
#include <algorithm>
#include "output.h"
using namespace std;
void cerr_command_usage();
void rheologyTest(int argc, char *const argv[], System &sy);
void testParameters(int argc, char *const argv[], System &sy);
void rod_bending(int argc, char *const argv[], System &sy);

int main (int argc, char *const argv[])
{
	System sy;
	sy.simulation = argv[1][0];
#ifdef TWODIMENSION
	cerr << "2D simulation\n";
#else
	cerr << "3D simulation\n";
#endif
	if (sy.simulation == 't'){
		testParameters(argc, argv, sy);
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

void preProcesses(System &sy,
				  vector <Particle *> &particle_active,
				  vector <Bond *> &bond_active){
	sy.check_active(particle_active, bond_active);
	sy.calcVolumeFraction();
	sy.outputConfiguration('n');
	while(sy.checkPercolation() == false){
		sy.shiftForPercolation();
		sy.makeNeighbor();
		sy.generateBond( bond_active );
	}
	sy.dt = sy.dt_max;
	if (sy.simulation == 's'){
		sy.strain_target = sy.step_strain_x;
		cerr << "strain target = " << sy.strain_target << endl;
		sy.strain_x = 0;
	} else {
		sy.vf_target = sy.volume_fraction * sy.volumefraction_increment;
		sy.stress_z = 0.;
	}
	sy.checkState(particle_active, bond_active);
	sy.outputData();
	sy.outputConfiguration('e');
}

void timeEvolution(System &sy,
				   vector <Particle *> &particle_active,
				   vector <Bond *> &bond_active){
	double t_next = sy.time + sy.interval_convergence_check*sy.dt_max;
	sy.calc_count = 0;
	sy.counter_relax_for_restructuring = 0;
	while ( sy.time < t_next ){
		/* Proceed applying strain when
		 * (1) the strain does not reach to the target value.
		 * (2) it is not relaxation process after bond breaking.
		 */
		if (sy.reachStrainTarget() ||
			sy.counter_relax_for_restructuring < sy.relax_for_restructuring){
			sy.prog_strain = false;
		} else {
			sy.prog_strain = true;
		}
		/* Main time evolution */
		if (sy.simulation == 'c') {
			sy.TimeDevStrainControlCompactionEuler(particle_active, bond_active);
		} else {
			sy.TimeDevStrainControlShearEuler(particle_active, bond_active);
		}
		/* Breakup of bonds */
		if (sy.counter_relax_for_restructuring++ >= sy.relax_for_restructuring){
			sy.checkBondFailure(bond_active);
			if (!sy.regeneration_bond.empty()){
				sy.outputRestructuring();
				sy.regeneration_onebyone();
				sy.counter_relax_for_restructuring = 0;
			}
			if (!sy.rupture_bond.empty()){
				sy.outputRestructuring();
				sy.rupture(bond_active);
				sy.counter_relax_for_restructuring = 0;
			}
		}
		/* Bond generation */
		if ( sy.calc_count % sy.interval_makeNeighbor == 0 ){
			sy.makeNeighbor();
		}
		sy.generateBond( bond_active );
		sy.wl[0]->addNewContact(particle_active);
		sy.wl[1]->addNewContact(particle_active);
		
		/* Counter */
		sy.time += sy.dt;
		sy.calc_count ++;
	}
	if (sy.simulation == 'c'){
		sy.calcVolumeFraction();
	}
}

void middleProcedures(System &sy,
					  vector <Particle *> &particle_active,
					  vector <Bond *> &bond_active){
	sy.simuAdjustment();
	if (sy.simulation == 's'){
		sy.calcShearStress();
	} else {
		sy.calcStress();
	}
	sy.checkState(particle_active, bond_active);
	if (sy.simulation == 's'){
		if (sy.diff_stress_x == 0){
			sy.stress_x_change = 0;
			sy.stress_z_change = 0;
		} else {
			sy.stress_x_change = abs(sy.stress_x - sy.stress_x_before)/sy.stress_x;
			sy.stress_z_change = abs(sy.stress_z - sy.stress_z_before)/sy.stress_z;
		}
		sy.strain_x = abs( sy.wl[0]->x - sy.wl[1]->x )/sy.lz;
		sy.strain_z = (sy.lz_init - sy.lz)/sy.lz_init;
		
	} else {
		sy.stress_z_change = abs(sy.stress_z - sy.stress_z_before)/sy.stress_z;
		sy.stress_x_change = abs(sy.stress_x - sy.stress_x_before)/sy.stress_x;
		sy.strain_z = (sy.lz_init - sy.lz)/sy.lz_init;
	}
	sy.stress_z_before = sy.stress_z;
	sy.stress_x_before = sy.stress_x;
	sy.makeNeighbor();
	sy.optimalTimeStep();
	if (sy.simulation == 's'){
		if ( sy.strain_x > sy.strain_x_last_output + 1.0 ){
			sy.outputConfiguration('n');
			sy.lz_last_output = sy.lz;
		}
	} else {
		if ( sy.lz < sy.lz_last_output - 1.0 ){
			sy.outputConfiguration('n');
			sy.lz_last_output = sy.lz;
		}
	}
	sy.output_log();
}

void strainControlShear(System &sy,
						vector <Particle *> &particle_active,
						vector <Bond *> &bond_active){
	preProcesses(sy, particle_active, bond_active);
	while( sy.strain_x < sy.max_strain_x){
		timeEvolution(sy, particle_active, bond_active);
		middleProcedures(sy, particle_active, bond_active);
		cerr << sy.strain_x << " --> " << sy.strain_target  << endl;
		if ( sy.reachStrainTarget() ){
			if (sy.mechanicalEquilibrium()){
				cerr <<"Equilibrium!"<< endl;
				/*  Output */
				sy.outputData();
				sy.outputConfiguration('e');
				sy.monitorDeformation('e');
				/*  Prepare next step */
				sy.strain_target += sy.step_strain_x;
			}
		}
	}
	return;
}

void strainControlCompression(System &sy,
							  vector <Particle *> &particle_active,
							  vector <Bond *> &bond_active)
{
	preProcesses(sy, particle_active, bond_active);
	while( sy.volume_fraction < sy.max_volume_fraction ){
		timeEvolution(sy, particle_active, bond_active);
		middleProcedures(sy, particle_active, bond_active);
		cerr << sy.volume_fraction << " --> " << sy.vf_target  << endl;
		if ( sy.reachStrainTarget() ){
			if ( sy.mechanicalEquilibrium() ){
				cerr <<"Equilibrium!"<< endl;
				/*  Output */
				sy.outputData();
				sy.outputConfiguration('e');
				sy.monitorDeformation('e');
				/*  Prepare next step */
				sy.vf_target *= sy.volumefraction_increment;
			}
		}
	}
	return;
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
	/* prepare calcuration */
	sy.initialprocess = true;
	sy.setSimulationViscosity();
	sy.preparationOutput();
	sy.setBondGenerationDistance(2.0);
	sy.setWall();
	sy.initGrid();
	//outputParameter(sy, fout_yap);
	sy.initDEM();
	cerr << "N " << sy.n_particle << endl;
	sy.makeInitialBond(2.0001);
	sy.initialprocess = false;
	
	/******************************/
#ifndef TWODIMENSION
	cerr <<  "3D simulation" << endl;
#else
	cerr <<  "2D simulation" << endl;
#endif
	/******************************/
	
	sy.checkState(sy.particle, sy.bond);
	
	/* particle_active and bond_active are effective.
	 * The others are not effective;
	 */
	/*****************************/
	vector <Particle *> particle_active;
	vector <Bond *> bond_active;
	if (sy.simulation == 'c'){
		strainControlCompression(sy,particle_active,bond_active);
	} else if (sy.simulation == 's'){
		strainControlShear(sy,particle_active,bond_active);
	}
	/*****************************/
	return;
}


void testParameters(int argc, char *const argv[], System &sy){
	sy.setParameterFile(argv[2]);
	/* read files */
	sy.readParameterFile();
	sy.readBondParameter();
	sy.lx = 50.;
	sy.ly = 50.;
	sy.lz = 50.;
	sy.lx0 = 0.5*sy.lx;
	sy.ly0 = 0.5*sy.ly;
	sy.lz0 = 0.5*sy.lz;
	
	sy.initialprocess = true;
	
	sy.particle.push_back(new Particle(0, vec3d(0.,0.,0.), 0, sy));
	sy.particle.push_back(new Particle(1, vec3d(2.,0.,0.), 0, sy));
	//particle.push_back(new Particle(2, vec3d(4,0,0), sy));
	sy.n_particle = 2;
	for (int i=0; i< sy.n_particle ; i++){
		sy.particle[i]->p.x += sy.lx0;
		sy.particle[i]->p.y += sy.ly0;
		sy.particle[i]->p.z += sy.lz0;
		sy.particle[i]->cerr();
	}
	sy.dt = sy.dt_max;
	// sy.setTimeStep(0.00001);
	/* prepare calcuration */
	sy.initGrid();
	sy.setBondGenerationDistance(2.0);
	sy.setSimulationViscosity();
	sy.preparationOutput();
	
	sy.initDEM();
	sy.makeInitialBond(2.0001);
	sy.outputYaplot();
	/******************************/
	
	vec3d rotate_axis;
	double angle ;
	///* tortion */
	//angle = M_PI/6;
	//rotate_axis.set(1,0,0);
	//sy.particle[1]->setRotate(rotate_axis,angle);
	///////////////////////
	
	/* bending */
	angle = M_PI/6;
	rotate_axis.set(0,-1.,0);
	sy.particle[1]->setRotate(rotate_axis,angle);
	///////////////////////
	/* sliding */
	//particle[1]->p.z += 0.5;
	/* elongation */
	// sy.particle[1]->p.x += -0.1;
	sy.particle[1]->p.z += 0.5;
	sy.outputYaplot();
	sy.initialprocess = false;
	
	foreach( vector<Bond *>, sy.bond, it_bond){
		(*it_bond)->addContactForce();
	}
	
	sy.particle[0]->resetForce();
	sy.particle[1]->resetForce();

	//particle[2]->p.set(4.0,0,0.5);
	cerr << "simulation start" << endl;
	
	for (int t = 0; t < 10000; t++){
		//sy.outputConfiguration('n');
		sy.outputYaplot();
		sy.simuAdjustment();
		foreach( vector<Bond *>, sy.bond, it_bond){
			(*it_bond)->addContactForce();
		}
		//   sy.particle[0]->resetForce();
		//////////////////////////////
		// calculate
		// x*_{j+1} and v*_{j+1}
		foreach( vector<Particle *>, sy.particle, it_particle){
			(*it_particle)->move_Euler();
		}
		cout << t << ' ' << sy.particle[1]->p.z << endl;
	}
	
	/* particle_active and bond_active are effective.
	 * The others are not effective;
	 */
	//sy.outputParameters();
	/*****************************/
	//compression_process(sy);
	//relaxation_process(sy);
	//shear_process(sy);
	/*****************************/
	
}
