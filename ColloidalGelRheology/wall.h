//
//  wall.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_wall_h
#define ColloidalGelRheology_wall_h

#define DELETE(x) if(x){delete [] x; x = NULL;}
#include "common.h"
#include "bond.h"
#include "particle.h"
#include "vec3d.h"
#include "system.h"
class Particle;
class System;

struct ConnectPointWall {
	int bond;
	int next;
	vec3d p;
};

enum WallType {top, bot};
class Wall{
private:
	System *sy;
	int objectID; 
	WallType walltype;	
	vec3d u;
	vector<int> neighbor;
    vector <Particle*> wall_particle;
    
    
	vec3d force;
    double stress_sc;
	vec3d a_velocity;
	vec3d a_velocity_1st;
	vec3d velocity_1st;
	int cn_size;
public:
	Wall(int, double zInit, System &sy);
	~Wall(){
	}
	double x;
	double y;
	double z;
	bool z_movable;
	bool x_movable;
	bool y_movable;
	vec3d velocity;
    vector <int> wall_group;
	int numConnectParticle(){
		return wall_particle.size();
	}
    void getParticles(vector<int> &particles_list);
    bool markWallConnected();
    
	void addNewContact(vector<Particle *> &particle_active);
    void initWallParticle(int i);
    
	void makeNeighbor();
	void addNeighbor( int neighbor_particle ){
		neighbor.push_back( neighbor_particle );
	}
	int valCn_size(){ return cn_size; }
	inline void resetForce(){force.reset() ;}
	void resetVelocity(){velocity.reset();}
	inline void stackForce(vec3d force_){
		force += force_;
	}
	void move();
	void output(ofstream &fout);
	inline vec3d *pu(){return &u;}
	void velocityLimit(){ velocity = velocity_1st;}
    void compactionStrainControl(double velocity);
    void shearingStrainControl(double velocity);
    double stressSensor_z();
    double stressSensor_x();
    
};
#endif
