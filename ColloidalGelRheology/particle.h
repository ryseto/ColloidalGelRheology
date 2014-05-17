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
#include "bond.h"
#include "vec3d.h"
#include "quaternion.h"
class Bond;
class System;
class Particle{
private:
	int cn_size;
	vec3d a_velocity;
	vec3d velocity;
#ifndef TWODIMENSION
	vec3d omega;
	vec3d a_omega;
	vec3d torque;
	vec3d d_rotation;
#else
	double omega;
	double a_omega;
	double torque;
	double d_rotation;
#endif
	vec3d force;
	vec3d p_mem;
	vec3d v_relative;
	vec3d omega_relative;
	vec3d p_pdcopy;
	vector<int> neighbor;
	ConnectPoint *cn;
	System *sy;
protected:
	void (vec3d::*p_change)(double, double, double);
	void remove_list_neighbor_all();
public:
	Particle(int particle_number_, const vec3d &position, const int initial_cluster,
			 System &sy);
	Particle(int paritcle_number_);
	~Particle();
	
	int particle_number;
	vec3d p;
	quaternion orientation;
	double orientation_angle = 0;
	int init_cluster;
	int wall_connected; //
	bool wall;
	bool near_boundary;
	double sq_force;
	
	void makeNeighbor();
//	void checkNearBoundary();
	void addNeighbor( int neighbor_particle ){
		neighbor.push_back( neighbor_particle );
	}
	inline void setBond(int bond_number, int next_particle){
		cn[cn_size].next = next_particle;
		cn[cn_size].bond = bond_number;
		cn[cn_size].tor_angle = 0.;
		cn_size++;
	}
	void addGravityForce();
	void move_Euler();
	inline vec3d vectorForce(){
		return force;
	}
	inline void resetForce(){
		force.reset();
#ifndef TWODIMENSION
		torque.reset();
#else
		torque = 0;
#endif
	}
	inline void resetFz(){
		force.z = 0;
	}
	inline void resetTorque(){
#ifndef TWODIMENSION
		torque.reset();
#else
		torque = 0;
#endif
	}
	inline double valForceZ(){
		double fz = force.z;
		force.reset();
		return  fz;
	}
#ifndef TWODIMENSION
	inline void stackForce(const vec3d &force_, const vec3d &torque_)
#else
	inline void stackForce(const vec3d &force_, const double &torque_)
#endif
	{
		force += force_;
		torque += torque_;
	}
	void calc_stack_Force();
	void setRotate(vec3d axis, const double angle);
	void setNorm_u(){
		for (int i = 0; i < cn_size ; i++){
			cn[i].u.unitvector();
		}
	}
	void setVelocityZero(){
		velocity.reset();
#ifndef TWODIMENSION
		omega.reset();
#else
		omega = 0;
#endif
	}
	void delConnectPoint(int bond_number);
	void z_shift( double dz ){ p.z += dz; }
	void x_shift( double dx );
	inline vec3d pos(){return p;}
	void setInitial(int disk_number_);
	inline double distOverlap(const vec3d &pp){
		return dist(p, pp) - 2.0; // ro = 2a = 1.0
	}
	vec3d *p_pos(){return &p; }
	vec3d *pu_back(){return &( cn[cn_size-1].u );}
	vec3d *pu(int i){return &( cn[i].u ); }
	double *p_tor_angle_back(){return &( cn[cn_size-1].tor_angle ); }
	double *p_tor_angle(int i){return &( cn[i].tor_angle ); }
	void generateBond();
	void output(ofstream &fout);
	void cerr(){
		std::cerr << "c " << p.x << ' ' << p.y << ' ' << p.z << std::endl;
	}
	double valForce(){ return force.norm(); }
	double valVelocity(){ return velocity.norm(); }
#ifndef TWODIMENSION
	double valOmega(){ return omega.norm(); }
	double kineticEnergy(){
		return 0.5*velocity.sq_norm() + 0.2*omega.sq_norm();
	}
#else
	double valOmega(){ return abs(omega);}
	double kineticEnergy(){
		return 0.5*velocity.sq_norm() + 0.2*omega*omega;
	}
#endif
	int valCn_size(){ return cn_size; }
	void zero_velocity(){
		velocity.reset();
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
