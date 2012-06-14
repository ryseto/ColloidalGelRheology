//
//  bond.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_bond_h
#define ColloidalGelRheology_bond_h

#include "common.h"
#include <cstdlib>
#include <cmath>
#include "contable.h"
#include "particle.h"
#include "system.h"
class Particle;
class System;

class Bond {
	int bond_number;
    BondParameter para;
    System *sy;
	double force_normal;
	vec3d r_vec;
    /*
     * force1 is always -force0
     */
    bool central_force;
	vec3d force0; 
	vec3d torque0;
	vec3d torque1;
	vec3d force_sliding;
	vec3d moment_sliding;
	vec3d moment_bending;
	double moment_torsion;
	vec3d u01;
    /*
     * This vector indicates the contact point
     * in the original coordinate of particle
     */
    vec3d u_vector[2]; 
    vec3d u_vector_initial[2]; 
    double q_old; // r-2
    vec3d d_slid;
    vec3d d_slid_old;
    vec3d sum_sliding;
    vec3d ang_bend;
    vec3d ang_bend_old;
    
    double ang_tort;
    double ang_tort_old;
    //	double damper_tortion;
    //    double damper_normal;
    //    vec3d damper_sliding;
    //    vec3d damper_bending;
    double kn_compression;
    double q;
    vec3d torque_tmp;
    
    int d[2];
	double *p_tor[2];
	vec3d *pp[2];
	vec3d *pu[2];
	Particle *p_particle0;
	Particle *p_particle1;
    
	double des[4];
    double sq_fsc;
    double sq_mbc;
    double sq_mtc;
    
private:
	void setting();
	void checkRuptureType(vec3d &, vec3d &, double);
	void periodicBoundary_rvec(vec3d &);
    
public:
	Bond(){}
	~Bond();
	Bond(int d0, int d1, System &sy);
    /* status; 
     * 0: The bond is ruptured.
     * 1: Only elastic deformation after the generation.
     * 2: Regenerations are experienced. 
     */
    int status; 
    /* initial_bond
     * ture: The bond is generated before the compression process.
     * false: The bond is generated during the compression process.
     */
    bool initial_bond;
    /* 
     * Distance between particle. 
     * The gap between particles or particle and wall are give by r - 2a.
     */
	double r;
    /*
     *  Evaluation of the rupture condition
     */
	double D_function;
    /*
     * "bond_active" vector is used to calculate efficiently.
     * This number memorize where it is in the vector.
     */
    int number_activebond;
    vec3d e_normal;
    
    
	void addContactForce();
	void calcForce();
    vec3d forceToParticle(int paritlce_number ){
        if (paritlce_number == d[0]){
            return force0;
        } else {
            return -force0;
        }        
    }
	void outputCompression(ofstream &out, double f_max);
	void outputTraction(ofstream &out, double f_max);
    
	void rupture();	
	void chPointer(int i, int particle_num);	
	void regeneration();
    
	inline double val_F_norm(){ return abs(force_normal);} 
    inline double val_F_slid(){ return force_sliding.norm();}
    inline double val_M_bend(){ return moment_bending.norm();}
    inline double val_M_tors(){ return abs(moment_torsion);}
    inline double valBendingAngle(){return ang_bend.norm();}
    inline double valTorsionalAngle(){return abs(ang_tort);}
    inline double valSlidingDisplacement(){return d_slid.norm();}
	inline double forceNorm(){return force0.norm();}
    
    void whichparticle(int &i, int &j);
    void cheackBondStress();
    
    vec3d position(int i){
        return *pp[i];
    }
	void monitor_state(ofstream &out);
    int val_bond_number(){
        return bond_number;
    }
};

#endif
