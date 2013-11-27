//
//  system.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_system_h
#define ColloidalGelRheology_system_h
#include "common.h"
#include "particle.h"
#include "grid.h"
#include "bond.h"
#include "wall.h"
#include "contable.h"
#include <vector>
#include <map>
using namespace std;
class Grid;
class Bond;
class Wall;

class System {
protected:
	/*
	 * Given parameters for simulation
	 */
	double dt_max;
	double max_move_step;
	double wall_velocity;
	double stress_change_convergence; // The change of stress during the interval_convergence_check.
	double stress_minimum;
	double initial_compaction;
	double max_volume_fraction;
	double step_strain_x;
	double max_strain_x;
	double dist_generate;
	double max_velocity_convergence;
	double max_ang_velocity_convergence;
	double diff_stress_convergence;
	double max_displacement;
	double eq_max_displacement;
	int interval_convergence_check;
	int interval_makeNeighbor;
	int interval_output_config;
	int relax_for_restructuring;
	int cnt_output_config;
	/*
	 * State parameters of system
	 */
	int n_particle; // Total number of particles
	vec3d force_wall;
	double time;
	double volume_fraction;
	double volume_fraction_init;
//	double area_fraction;
	double strain_x;
	double stress_x;
	double strain_z;
	double stress_z;
	double stress_z_top;
	double stress_z_bot;
	double stress_y;
	double kinetic_energy;
	int counterBreak;
	int counterRegenerate;
	double ave_force;
	double ave_bondforce;
	double max_force;
	double max_velocity;
	double max_ang_velocity;
	double r_min;
	double r_max;
	double sliding_disp_max;
	double bending_angle_max;
	double torsional_angle_max;
	double force_max;
	double average_contact_number;
	/*
	 * Simulation parameters
	 */
	double eta_factor; // relative viscosity used in simulation
	double diff_stress_x;
	double diff_stress_z;
	double stress_z_before;
	double stress_x_before;
	double stress_x_change;
	double stress_z_change;
	int counterRegenerate_before;
	/*
	 * variables for simulations
	 */
	bool prog_strain;
	int counter_relax_for_restructuring;
	int min_relaxation_loop;
	int max_relaxation_loop;
	int cnt_loop;

	double volumefraction_increment;
	double strain_target; // next equilibrium for shear
	double vf_target;// next equilibrium for compaction
	double lz_last_output;
	double strain_x_last_output;
	string bond0_file;
	string bond1_file;
	char fn_common[64];   // 192 // 128 + 64
	char parameters_file[128];
	char parameters[32];
	char init_cluster_file[256]; //256
	char init_cluster[256]; //256
	vector <vec3d> pos_previous;
	string version;
	ofstream fout_data;
	ofstream fout_log;
	ofstream fout_yap;
	ofstream fout_conf;
	ofstream fout_deform;
	vec3d f_ex;
private:
	void preProcesses();
	void timeEvolution();
	void middleProcedures();
	void shiftCenterOfMass(vector<vec3d> &p);
	void setWall();
	void setBondGenerationDistance(double);
	void generateBond();
	void generateBondAll();
	void makeInitialBond(double generation_distance);
	void checkState();
	void checkBondFailure();
	void regeneration();
//	void regeneration_onebyone();
	void rupture();
	void TimeDevStrainControlCompactionEuler();
	void TimeDevStrainControlShearEuler();
	void TimeDevBendingTest();

	bool mechanicalEquilibrium();
	bool reachStrainTarget();
	void setTarget();
	void setFirstTarget();
	bool checkEndCondition();
	void checkActiveElements();
	bool checkPercolation();
	void shiftForPercolation();
	void calcStress();
	void calcShearStress();
	void checkStressDiff();
	bool checkOutputStrain(double);
	void simuAdjustment();
	void setSimulationViscosity();
	void preparationOutput();
	void calcVolumeFraction();
	void calcLocalStrains();
	void setUnitForce(double);
//	void setWallStress(double, double, double);
	void setGravityConstant(double);
	void renew_Lz();
	void optimalTimeStep();
	void initGrid();
	void reset_strain_x();
	void makeNeighbor();
	void monitorDeformation(char equilibrium);
	void outputData();
	void output_log();
	void outputYaplot();
	void outputConfiguration(char equilibrium);
	void outputRestructuring();
	void outlog();

public:
	System();
	~System();
	void importPositions();
	void initDEM();
	void readParameterFile();
	void readParameter(const string &codeword, const string &value);
	void readBondParameter();
	void setParameterFile(char *);
	void setInitClusterFile(string);
	void setVersion(string );
	void strainControlSimulation();
	void bendingSimulation();
	void setExternalforce(double fex){
		f_ex.set(fex,0,0);
	}

	void setLinearChain(int number_of_particles);

	/*
	 * Objects
	 */
	vector<Particle *> particle;
	vector <Particle *> particle_active;
	vector<Bond *> bond;
	vector <Bond *> bond_active;
	vector <Wall *> wl;
	/*
	 * Objects for simulation
	 */
	Grid *grid;
	ConTable *ct;
	vector<int> regeneration_bond;
	vector<int> rupture_bond;
	/*
	 * System parameters
	 */
	double lx, ly, lz;
	double lz_init;

	double dt; // Time step to integrate equatino of motion.
	double eta;
	double eta_rot;
	/*
	 * State of system
	 */
	int n_bond; // The total number of bond including breakup bond
	bool percolation;
	int rup_normal; // Counters for rupture events
	int rup_shear;
	int rup_bend;	
	int rup_torsion;
	/*
	 * Simulation parameters
	 */
	char simulation;
	bool initialprocess; // To prepare initial bonds and so on.
	double lx0, ly0, lz0; // = (lx/2, ly/2, lz/2)
	int n_top;
	int n_bot;
	BondParameter bond0;
	BondParameter bond1;
	double sq_dist_generate;
};

#endif
