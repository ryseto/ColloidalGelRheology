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
class ConTable;


class System {
protected:
	void shiftCenterOfMass(vector<vec3d> &p);
	ofstream fout_data;
	ofstream fout_log;
	ofstream fout_yap;
	ofstream fout_conf;
    ofstream fout_deform;
    string bond0_file;
    string bond1_file;
    double eta_factor;
    double kinetic_energy;
public:
	System();
	~System();
 	void importPositions();
	void setRod(int m);
    void setWall();
    void initDEM();
    
	void readParameterFile();
    void readParameter(const string &codeword, const string &value);
    void readBondParameter();
	void setParameterFile(char *);
	void setInitClusterFile(char *);
	void setVersion(char *);
	void setBondGenerationDistance(double);
    void generateBond(vector<Bond *> &bond_active);
    void generateBond();
    void makeInitialBond(double generation_distance);
    
    void checkState(vector <Particle *> &particle_active, vector <Bond *> &bond_active);
    void checkBondFailure(vector<Bond *> &bond_active);
    
    void regeneration();
    void regeneration_onebyone();
    void rupture(vector<Bond *> &bond_active);
    void TimeDevStrainControlCompactionEuler(vector<Particle *> &particle_active, 
                                             vector<Bond *> &bond_active,
                                             bool wall_move);
    
    void TimeDevStrainControlShearEuler(vector<Particle *> &particle_active, 
                                        vector<Bond *> &bond_active,
                                        bool wall_move);
    void check_active(vector<Particle *> &particle_active, vector<Bond *> &bond_active);
    bool checkPercolation();
    void shiftForPercolation();
    void calcStress();
    void calcShearStress();
    void simuAdjustment();
	void setSimulationViscosity();
	void preparationOutput();
    void calcVolumeFraction();
    void calcLocalStrains();
	void setUnitForce(double);
	void setWallStress(double, double, double);
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

    
    //	void recordEquibiruimVolumeFraction(){
    //		volume_fraction_equilibrium.push_back(volume_fraction);
    //	}
	/*
	 * System parameter 
	 *	
	 *
	 */ 
    ConTable *ct;
    vector<Bond *> bond;
    vector<Particle *> particle;
    vector<Wall *> wl;
    Wall *wl_top;
    Wall *wl_bot;
    Grid *grid;
    vector<int> regeneration_bond;
	vector<int> rupture_bond;
	vector <double> volume_fraction_equilibrium;// It should be removed.    
	// number of particle
	char simulation;
	int n_particle; //4
    BondParameter bond0;
    BondParameter bond1;
    /* 
     * For distinguish initial bond generation.
     */
    bool initialprocess;  
	char fn_common[64];   // 192 // 128 + 64
	char parameters_file[128];
	char parameters[32];
	char init_cluster_file[128]; //256
	char init_cluster[32]; //256
	char version[3];
	/*
	 * State variables 
	 *	
	 *
	 */ 
	double time;
    double volume_fraction;
    double area_fraction;
    //    double modified_volume_fraction;
    int calc_count;
	int n_bond; // The total number of bond including breakup bond
	int counterBreak;
	int counterRegenerate;   // 4
	int rup_normal;
	int rup_shear;
	int rup_bend;
	int rup_torsion;
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
    char string_L[3];
	double lx;
	double ly;
	double lz;
	double lx0;
	double ly0;
	double lz0;
	
	//double init_aggregate_radius;
	/*
	 * calucuration parameters
	 */
	double dt;
	double dt_max;
    double max_move_step;
    
    double eta;
	double eta_rot;
	double dist_generate;
	double sq_dist_generate;
	//double cp_aggregate_radius;
	/*
	 * compression 
	 */
	int n_top;
	int n_bot;
	double lz_init;
    bool percolation;
	vec3d force_wall;
    double wall_velocity;
    
    double strain_x;
    double stress_x;
    double strain_z;
	double stress_z;
    
    
	double stress_z_initial;
	double stress_z_increment_ratio;
    double stress_x_initial;
	double stress_x_increment_ratio;
	double stress_y;
    // The diffrece of stresses between the top and bottom.
    double diff_stress_x;
    double diff_stress_z;
    // The change of stress during the interval_convergence_check.
    double stress_x_change; 
    double stress_z_change;
    double stress_change_convergence;
    double stress_minimum;
    double initial_compaction;
	double max_volume_fraction;
    double step_strain_x;
    double max_strain_x;
    
	int relax_for_restructuring;
    double volumefraction_increment;
    double max_velocity_convergence;
    double max_ang_velocity_convergence;
    double diff_stress_convergence;
	int interval_convergence_check;
	double strain_x_change;
	double strain_z_change;
    double dLz_outputconfig;    
};

#endif
