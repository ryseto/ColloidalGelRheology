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
#include "contable.h"
#include "StressTensor.h"
#include <vector>
#include <map>
using namespace std;
class Grid;
class Bond;
class Particle;

class System {
protected:
	/*
	 * Given parameters for simulation
	 */
	double dt_max;
	double max_move_step;
	double compaction_speed;
	double shear_strain; 
	double stress_change_convergence; // The change of stress during the interval_convergence_check.
	double stress_minimum;
	double max_volume_fraction;
	double step_shear_strain;
	double max_shear_strain;
	double dist_generate;
	double max_velocity_convergence;
	double max_ang_velocity_convergence;
	double diff_stress_convergence;
	double max_displacement;
	double eq_max_displacement;
	int interval_convergence_check;
	int interval_makeNeighbor;
	int interval_output_config;
	int relax_for_restructuring; //@@@ I should chekc the effect @@@
	int cnt_output_config;
	/*
	 * State parameters of system
	 */
	int n_particle; // Total number of particles
	double time;
	double volume_fraction;
	double volume_fraction_init;
	double strain_z;
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
	int counterRegenerate_before;
	/*
	 * variables for simulations
	 */
	bool prog_strain;
	/*
	 * deformation_type
	 * 0: uniaxial compression
	 * 1: biaxial compression
	 */
	int deformation_type;
	int counter_relax_for_restructuring;
	int min_relaxation_loop;
	int max_relaxation_loop;
	int cnt_loop;
	double first_volumefraction;
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
	ofstream fout_particle;
	ofstream fout_bond;
	vec3d f_ex;
private:
	void (System::*enforceStrain)();
	void affineUniaxialCompression();
	void affineBiaxialCompression();
	void bounadyUniaxialCompression();
	void bounadyBiaxialCompression();
	void preProcesses();
	void timeEvolution();
	void middleProcedures();
	void shiftCenterOfMass(vector<vec3d> &p);
	void setBondGenerationDistance(double);
	void generateBond();
	void makeInitialBond(double generation_distance);
	void checkState();
	void checkBondFailure();
	void regeneration();
	void rupture();
	void TimeDevPeriodicBoundaryCompactionEuler();
	void TimeDevBendingTest();
	bool mechanicalEquilibrium();
	void setTarget();
	void setFirstTarget();
	bool checkEndCondition();
	bool checkPercolation();
	void shiftForPercolation();
	void calcStress();
	void calcShearStress();
	bool checkOutputStrain(double);
	void simuAdjustment();
	void setSimulationViscosity();
	void preparationOutput();
	void calcVolumeFraction();
	void calcLocalStrains();
	void setUnitForce(double);
	void setGravityConstant(double);
	void renew_Lz();
	void optimalTimeStep();
	void initGrid();
	void reset_strain_x();
	void makeNeighbor();
	void makeNeighborPB();
	void monitorDeformation(char equilibrium);
	void outputData();
	void output_log();
	void outputYaplot();
	void outputConfiguration(char equilibrium);
	void outputRestructuring();
	void outlog();
	double system_volume();
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
	vector<Bond *> bond;
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
	double lx;
	double ly;
	double lz;
	double lx_half;
	double ly_half;
	double lz_half;
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
	string simulation;
	bool initialprocess; // To prepare initial bonds and so on.
	int n_top;
	int n_bot;
	BondParameter bond0;
	BondParameter bond1;
	double sq_dist_generate;
	StressTensor total_contact_stressXF;
};

#endif
