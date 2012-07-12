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
    //    particle.push_back(new Particle(2, vec3d(4,0,0), sy));
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

void cerr_command_usage()
{
	cerr << "usage of ColloidalGelRheology" << endl;
	cerr << "arg: c/s ";
	cerr << "parameter_file ";
	cerr << "initial_file ";
	cerr << "version ";
	cerr << endl;
}

void strainControlShear(System &sy){
    vector <Particle *> particle_active;
	vector <Bond *> bond_active;
	sy.check_active(particle_active, bond_active);
    sy.calcVolumeFraction();
    sy.outputConfiguration('n');
    while(sy.checkPercolation() == false){
        sy.shiftForPercolation();
        sy.makeNeighbor();
        sy.generateBond( bond_active );
    }
    sy.outputConfiguration('e');
	sy.makeNeighbor(); // To be checked.
    double strain_x_equilibrium = sy.step_strain_x;
    bool prog_strain = true;
    sy.dt = sy.dt_max;
    sy.strain_x = 0;
    int counter_relax_for_restructuring = 0;
    double stress_x_before = sy.stress_x;
    double stress_z_before = sy.stress_z;
    double strain_x_before = sy.strain_x;
    sy.checkState(particle_active, bond_active);
    sy.outputData();
	while( sy.strain_x < sy.max_strain_x){
        ////////////////////////////////////////////////////////////// 
        //// The information of state is kept to be compare 
        //// after time iteration.
		int counterRegenerate_before = sy.counterRegenerate;
        stress_x_before = sy.stress_x;
        stress_z_before = sy.stress_z;
        strain_x_before = sy.strain_x;
        ///////////////////////////////////////////////////////////////
        //// Set up the simulation of this time interval.
        double t_next = sy.time + sy.interval_convergence_check*sy.dt_max;
        sy.calc_count = 0;
        while ( sy.time < t_next ){
            if (counter_relax_for_restructuring > sy.relax_for_restructuring){
                prog_strain = true;
            } else{
                counter_relax_for_restructuring ++;
                prog_strain = false;
            }
            if ( sy.strain_x >= strain_x_equilibrium ){
                prog_strain = false;
            }
            sy.TimeDevStrainControlShearEuler(particle_active, 
                                              bond_active, prog_strain);
            if (counter_relax_for_restructuring > sy.relax_for_restructuring){
                sy.checkBondFailure(bond_active);
                if (!sy.regeneration_bond.empty()){
                    //sy.outputRestructuring();
                    sy.regeneration_onebyone();
                    counter_relax_for_restructuring = 0;
                    
                }
                if (!sy.rupture_bond.empty()){
                    //sy.outputRestructuring();
                    sy.rupture(bond_active);
                    counter_relax_for_restructuring = 0;
                }
            }
            sy.generateBond( bond_active );
            sy.wl[0]->addNewContact(particle_active);            
            sy.wl[1]->addNewContact(particle_active);
            if ( sy.calc_count % 1000 == 0 ){
                sy.makeNeighbor();
            }
            sy.time += sy.dt;
            sy.calc_count ++;
        }
        ///////////////////////////////////////////////////////////////
        sy.simuAdjustment();
        sy.checkState(particle_active, bond_active);
        sy.calcShearStress();
        if (sy.diff_stress_x == 0){
            sy.stress_x_change = 0;
            sy.stress_z_change = 0;
        } else {
            sy.stress_x_change = abs(sy.stress_x - stress_x_before)/sy.stress_x;
            sy.stress_z_change = abs(sy.stress_z - stress_z_before)/sy.stress_z;
        }
        sy.strain_x = abs( sy.wl[0]->x - sy.wl[1]->x )/sy.lz;
		sy.strain_z = (sy.lz_init - sy.lz)/sy.lz_init;
		sy.output_log();
		sy.makeNeighbor();
        sy.optimalTimeStep();
		if ( sy.strain_x >= strain_x_equilibrium ){
            /* Check the equilibrium condition
             * (1.1) The difference between the stresses at the two walls is enough small.
             * (1.2) The stress is less than the minimum stress.
             * (2) The change of the measured stress is enough small.
             * (3) The maximum velocity of the particles are enough small.
             * (4) The maximum angular velocity of the particles are enough small.
             *  ---------> Bond breakup requires (1)-(4).
             * (5) No restructuring.
             *  ---------> The equilibrium requires (1)-(5).
             */
            if ( (sy.diff_stress_x < sy.diff_stress_convergence
                  ||  sy.stress_x < sy.stress_minimum)
                && sy.stress_x_change < sy.stress_change_convergence
                && sy.max_velocity < sy.max_velocity_convergence
                && sy.max_ang_velocity < sy.max_ang_velocity_convergence
                && counterRegenerate_before == sy.counterRegenerate){
                cerr <<"Equilibrium!"<< endl;
                sy.outputData();
                sy.outputConfiguration('e');
                sy.monitorDeformation('e');
                /////////////////////////
                //counter_relax_for_restructuring = 0;
                strain_x_equilibrium += sy.step_strain_x;
                prog_strain = true;
            }
		} 
	}
    return;
}

void strainControlCompression(System &sy)
{
    vector <Particle *> particle_active;
	vector <Bond *> bond_active;
	sy.check_active(particle_active, bond_active);
	sy.makeNeighbor();
    sy.calcVolumeFraction();
    // sy.calcModifiedVolumeFraction();
    // double vf_equilibrium = sy.volume_fraction / sy.initial_compaction;
    double vf_equilibrium = sy.volume_fraction * sy.volumefraction_increment;     
    
    cerr << "phi = " << sy.volume_fraction << " " << vf_equilibrium<< endl;
    double lz_checked_output = sy.lz;
    bool compaction = true;
    sy.dt = sy.dt_max;
    sy.outputConfiguration('n');
    sy.stress_z = 0.;
    int counter_relax_for_restructuring = 0;
    int limit_stopnumber = 0;
    sy.checkState(particle_active, bond_active);
    
    if (sy.checkPercolation() == false){
        sy.shiftForPercolation();
        sy.makeNeighbor();
        sy.generateBond( bond_active );
        sy.outputConfiguration('n');
        exit(1);
    }
	while( sy.volume_fraction < sy.max_volume_fraction ){
        ////////////////////////////////////////////////////////////// 
        //// The information of state is kept to be compare 
        //// after time iteration.
        double stress_z_before = sy.stress_z;
        double stress_x_before = sy.stress_x;
		int counterRegenerate_before = sy.counterRegenerate;
        ///////////////////////////////////////////////////////////////
        //// Set up the simulation of this time interval.
        //// 
        double t_next = sy.time + sy.interval_convergence_check*sy.dt_max;
        sy.calc_count = 0;
        
        while ( sy.time < t_next ){
            if (counter_relax_for_restructuring ++ > sy.relax_for_restructuring){
                compaction = true;
                limit_stopnumber = 0;
            } else{
                compaction = false;
                if (limit_stopnumber++ > 10){
                    compaction = true;
                }
            }
            if (sy.volume_fraction > vf_equilibrium){
                compaction = false;
            }
            
            sy.TimeDevStrainControlCompactionEuler(particle_active, bond_active, compaction);
            if (counter_relax_for_restructuring > sy.relax_for_restructuring){
                sy.checkBondFailure(bond_active);
                if (!sy.regeneration_bond.empty()){
                    sy.regeneration();
                    sy.regeneration_bond.clear();
                    counter_relax_for_restructuring = 0;
                }
                if (!sy.rupture_bond.empty()){
                    sy.rupture(bond_active);
                    sy.rupture_bond.clear();
                    counter_relax_for_restructuring = 0;
                }
            }
            
            if ( sy.calc_count % 100 == 0 ){
                sy.makeNeighbor();
            }
            
            sy.generateBond( bond_active );
            sy.wl[0]->addNewContact(particle_active);            
            sy.wl[1]->addNewContact(particle_active);
            sy.calcVolumeFraction();
            sy.time += sy.dt;
            sy.calc_count ++;
        }
        
        ////
        ///////////////////////////////////////////////////////////////
        //sy.calcVolumeFraction();
        sy.calcStress();
        sy.checkState(particle_active, bond_active);
        sy.simuAdjustment();
        
        sy.stress_z_change = abs(sy.stress_z - stress_z_before)/sy.stress_z;
        sy.stress_x_change = abs(sy.stress_x - stress_x_before)/sy.stress_x;
        sy.strain_z = (sy.lz_init - sy.lz)/sy.lz_init;
        
        
		sy.output_log();
		sy.makeNeighbor();
        sy.optimalTimeStep();
        if ( sy.lz < lz_checked_output - 1.0 ){
            sy.outputConfiguration('n');
            lz_checked_output = sy.lz;
        }
        cerr << sy.volume_fraction << " --> " << vf_equilibrium  << endl;
		if ( sy.volume_fraction >= vf_equilibrium ){
            cerr << "check equilibrium" << endl;
            /* Check the equilibrium condition
             * (1) The change of the measured stress is enough small.
             * (2) The difference between the stresses at the two walls is enough small.
             * (3) The maximum velocity of the particles are enough small.
             * (4) The maximum angular velocity of the particles are enough small.
             *  ---------> Bond breakup requires (1)-(4).
             * (5) No restructuring.
             *  ---------> The equilibrium requires (1)-(5).
             */
            if ((sy.stress_z_change < sy.stress_change_convergence
                 && sy.diff_stress_z < sy.diff_stress_convergence)
                || sy.stress_z_change < 1e-8 ){
                if( sy.max_velocity < sy.max_velocity_convergence
                   && sy.max_ang_velocity < sy.max_ang_velocity_convergence
                   && counterRegenerate_before == sy.counterRegenerate){
                    cerr <<"Equilibrium!"<< endl;
                    sy.outputData();
                    sy.outputConfiguration('e');
                    sy.monitorDeformation('e');
                    counter_relax_for_restructuring = 0;
                    vf_equilibrium = vf_equilibrium* sy.volumefraction_increment;     
                    compaction = true;
                }
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
    
    if (sy.simulation == 'c'){
        strainControlCompression(sy);
    } else if (sy.simulation == 's'){
        //        sy.max_volume_fraction = 0.3;
        //        strainControlCompression(sy);
        strainControlShear(sy);
    }
    //stressControlCompression(sy);
    
    
	//relaxation_process(sy);
	//shear_process(sy);
	/*****************************/
    return;
}
