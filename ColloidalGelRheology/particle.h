//
//  particle.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_particle_h
#define ColloidalGelRheology_particle_h

#include "common.h"
//#include "grid.h"
#include "bond.h"
#include "vec3d.h"
#include "quaternion.h"
//class Grid;
class Bond;
class System;
//extern Grid grid;

class Particle{	
private:
	int cn_size;
	vec3d a_velocity;
	vec3d a_omega;
	vec3d a_velocity_1st;
	vec3d a_omega_1st;
	vec3d velocity;
	vec3d omega;
	vec3d velocity_1st;
	vec3d omega_1st;
	vec3d force;
	vec3d torque;
	vec3d v_old;
	vec3d omega_old;
	vec3d p_mem;
	vec3d d_rotation;
	vec3d v_relative;
	vec3d omega_relative;
	vec3d p_pdcopy;
	vector<int> neighbor;
	ConnectPoint cn[16];
	System *sy;
	double z0;
protected:
	void (vec3d::*p_change)(double, double, double);
	void remove_list_neighbor_all();
public:
	Particle(int particle_number_, const vec3d &position, const int initial_cluster,
             System &sy);
	Particle(int paritcle_number_);
	~Particle(){
	}
	
	vec3d p;
    quaternion orientation;
    int init_cluster;
    int wall_connected;    
	bool wall;
	bool near_boundary;
    bool after_rupture;
	double sq_force;
	int particle_number;
    
	void makeNeighbor();
    void checkNearBoundary();
	void addNeighbor( int neighbor_particle ){
		neighbor.push_back( neighbor_particle );
	}
	inline void setBond(int bond_number, int next_particle){
		cn[cn_size].next = next_particle;
		cn[cn_size].bond = bond_number;
		cn[cn_size].tor_angle = 0.;
		++ cn_size;
	}
    
	void addGravityForce();
    void move_Euler();
    
	inline vec3d vectorForce(){
		return force;
	}
	inline void resetForce(){
		force.reset();
		torque.reset();
	}        
	inline void resetFz(){
		force.z = 0;
	}
	
	inline void resetTorque(){
		torque.reset();
	}
	inline double valForceZ(){
		double fz = force.z;
		force.reset();
		return  fz;
	}
    
	inline void stackForce(const vec3d &force_, const vec3d &torque_){
		force += force_;
		torque += torque_;
	}
    void calc_stack_Force();
    
    void setRotate(vec3d axis, const double angle);
    
	void setNorm_u(){
		for (int i = 0; i < cn_size ; ++i){
			cn[i].u.unitvector();
		}
	}
    
	void setVelocityZero(){
		velocity.reset();
		omega.reset();
	}
	void delConnectPoint(int bond_number);	
	inline double x(){return p.x;}
	inline double y(){return p.y;}
	inline double z(){return p.z;}	
    
	void z_shift( double dz ){ p.z += dz; }
    
    void x_shift( double dx );
    
	inline vec3d pos(){return p;}	
	void setInitial(int disk_number_);
	void setPosition(const vec3d &position);
	inline double distOverlap(const vec3d &pp){
		return dist(p, pp) - 2.0; // ro = 2a = 1.0
	}
	vec3d *p_pos(){return &p; }
	vec3d *pu_back(){return &( cn[cn_size-1].u );}
	vec3d *pu(int i){return &( cn[i].u ); }
	double *p_tor_angle_back(){return &( cn[cn_size-1].tor_angle ); }
	double *p_tor_angle(int i){return &( cn[i].tor_angle ); }
	
	void generateBond();
    //	void generateBond_All();
    
	//void wallInclude(){ wall = true; };
	void output(ofstream &fout);
    
	void cerr(){
		std::cerr << "c " << p.x << ' ' << p.y << ' ' << p.z << std::endl;
		//fout << particle_number << ' ' << p.x << ' ' << p.y << ' ' << p.z << endl;
	}	
	double valForce(){ return force.norm(); }
	double valTorque(){ return torque.norm(); }
	double valVelocity(){ return velocity.norm(); }
	double valOmega(){ return omega.norm(); }
	double kineticEnergy(){
        return 0.5*velocity.sq_norm() + 0.2*omega.sq_norm();
    }
	int valCn_size(){ return cn_size; }
	void zero_velocity(){
		velocity.reset();
		velocity_1st.reset();
	}
	void addForce(double fx, double fy, double fz){
		force.x += fx;
		force.y += fy;
		force.z += fz;
	}
	bool percolate(vector<int> perco_path);
    
    void markWallConnected(int wt, vector<int> &wall_group);
    void outputBond();

    
};



#endif
