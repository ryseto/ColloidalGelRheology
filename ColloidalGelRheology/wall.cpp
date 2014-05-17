//
//  wall.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

//#include <iostream>

#include "wall.h"

Wall::Wall(int topbot, double _z, System &sy_): z(_z){
	sy = &sy_;
	if ( topbot == sy->n_top){
		walltype = top;
		u.set(0., 0., -1.);
	} else if ( topbot ==	sy->n_bot){
		walltype = bot;
		u.set(0., 0., 1.);
	}
	velocity.reset();
	x = sy->lx_half;
	y = sy->ly_half;
}

void Wall::makeNeighbor(){
	vector<int> *neighbor_cell;
	if (walltype == top){
		neighbor_cell = sy->grid->get_neighbor_list_pointer_top();
	} else {
		neighbor_cell = sy->grid->get_neighbor_list_pointer_bot();
	}
	neighbor.clear();
	foreach( vector<int> , *neighbor_cell, i){
		if ( sy->particle[ *i ]->wall == false )
			neighbor.push_back( *i );
	}
}

void Wall::addNewContact(vector<Particle *> &particle_active){
	bool generated = false;
	unsigned long n_neighbor = neighbor.size();
	for(int j=0; j < n_neighbor; j++){
		if (  u.z*(sy->particle[neighbor[j]]->p.z - z) <= 1.0 ){
			if ( sy->particle[ neighbor[j] ]->wall == false){

				wall_particle.push_back( sy->particle[ neighbor[j] ] );
				
				sy->particle[ neighbor[j] ]->wall = true;
				neighbor[j] = neighbor.back();
				neighbor.pop_back();
				n_neighbor --;
				generated = true;
				j--;
			}
		}
	}
	if (generated){
		/*
		 * I should use STL algorithm.
		 */
		
		unsigned long n_particle_active = particle_active.size();
		for (int i=0; i< n_particle_active;  i++){
			if ( particle_active[i]->wall ){
				particle_active[i] = particle_active.back();
				particle_active.pop_back();
				n_particle_active --;
				i --;
			}
		}
	}
}

void Wall::initWallParticle(int i){
	if ( sy->particle[i]->wall == false){
		wall_particle.push_back( sy->particle[i] );
		sy->particle[i]->wall = true;
	}
}

void Wall::compactionStrainControl(double velocity){
	double dz = velocity*sy->dt;
	z += dz;
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->p.z += dz;
	}
	if (walltype == top){
		sy->lz += dz;
	} else {
		sy->lz -= dz;
	}
}

void Wall::moveWallGroup(vec3d shift){
	for (int i = 0; i < wall_group.size(); i++){
		int k = wall_group[i];
		sy->particle[k]->p += shift;
		if (sy->particle[k]->p.x > sy->lx){
			sy->particle[k]->p.x -= sy->lx;
		}
	}
	x += shift.x;
	y += shift.y;
	z += shift.z;
	if ( walltype == top){
		sy->lz += shift.z;
	} else {
		sy->lz -= shift.z;
	}
}

void Wall::shearingStrainControl(double velocity){
	double dx = velocity*sy->dt;
	x += dx;
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->p.x += dx;
		if ( (*iter_p)->p.x < 0 )
			(*iter_p)->p.x +=  sy->lx;
		else if ( (*iter_p)->p.x >=  sy->lx )
			(*iter_p)->p.x -=  sy->lx;
	}
}


void Wall::stressSensor(double &stress_x, double &stress_z){
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->resetForce();
	}
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->calc_stack_Force();
	}
	double force_z = 0;
	double force_x = 0;
	foreach( vector <Particle* >, wall_particle, iter_p){
		force_z += ((*iter_p)->vectorForce()).z ;
		force_x += ((*iter_p)->vectorForce()).x ;
	}
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->resetForce();
	}
	
#ifdef TWODIMENSION
	// 2D
	stress_x = force_x / (sy->lx*2.0);
	stress_z = force_z / (sy->lx*2.0);
#else
	// 3D
	stress_x = force_x / (sy->lx*sy->ly);
	stress_z = force_z / (sy->lx*sy->ly);
#endif
}


double Wall::stressSensor_z(){
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->resetForce();
	}
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->calc_stack_Force();
	}
	double force_z = 0;
	foreach( vector <Particle* >, wall_particle, iter_p){
		force_z += ((*iter_p)->vectorForce()).z ;
	}
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->resetForce();
	}
#ifdef TWODIMENSION
	// 2D
	stress_sc = force_z / (sy->lx*2.0);
#else
	// 3D
	stress_sc = force_z / (sy->lx*sy->ly);
#endif
	return stress_sc;
}

double Wall::stressSensor_x(){
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->resetForce();
	}
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->calc_stack_Force();
	}
	double force_x = 0;
	foreach( vector <Particle* >, wall_particle, iter_p){
		force_x += ((*iter_p)->vectorForce()).x ;
	}
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->resetForce();
	}
	//cerr << force_x << endl;
#ifdef TWODIMENSION
	// 2D
	stress_sc = force_x / (sy->lx*2.0);
#else
	// 3D
	stress_sc = force_x / (sy->lx*sy->ly);
#endif
	return stress_sc;
}

void Wall::getParticles(vector<int> &particles_list){
	foreach( vector <Particle* >, wall_particle, iter_p){
		//wall_particle.push_back();
		particles_list.push_back((*iter_p)->particle_number);
		
	}
}

void Wall::markWallConnected(){
	int wt;
	if (walltype == bot){
		wt = 1;  // bt
	} else {
		wt = 2; // top
	}
	//	bool percolation = false;
	foreach( vector <Particle* >, wall_particle, iter_p){
		(*iter_p)->markWallConnected(wt, wall_group);
		//if ( sy->percolation == true){
		//	percolation = true;
		//}
	}
	cerr << "wall" << wt << " : "  << wall_group.size() << endl;
	//return false;
}


void Wall::output(ofstream &fout){
	double x_wall = x;
	while (x_wall > sy->lx){
		x_wall += - sy->lx;
	}
	while ( x_wall < 0 ){
		x_wall +=  sy->lx;
	}
	
	double zz = z - sy->lz_half;
	fout << "l " <<   -sy->lx_half << ' ' <<-sy->ly_half << ' ' << zz;
	fout <<  ' ' << sy->lx-sy->lx_half << ' ' <<-sy->ly_half << ' ' << zz << endl;
	fout << "l " << sy->lx-sy->lx_half << ' ' <<-sy->ly_half << ' ' << zz ;
	fout <<  ' ' << sy->lx-sy->lx_half << ' ' << sy->ly -sy->ly_half << ' ' << zz << endl;
	fout << "l " << sy->lx-sy->lx_half << ' ' << sy->ly -sy->ly_half << ' ' << zz ;
	fout <<  ' ' <<   -sy->lx_half     << ' ' << sy->ly -sy->ly_half << ' ' << zz << endl;
	fout << "l " <<   -sy->lx_half     << ' ' << sy->ly -sy->ly_half << ' ' << zz ;
	fout <<  ' ' <<   -sy->lx_half     << ' ' <<-sy->ly_half << ' ' << zz << endl;
	fout << "l " << x_wall-sy->lx_half << ' ' <<-sy->ly_half << ' ' << zz;
	fout << ' '  <<	x_wall-sy->lx_half << ' ' << sy->ly -sy->ly_half << ' ' << zz << endl;
	fout << "l " <<   -sy->lx_half     << ' ' << y -sy->ly_half << ' ' << zz;
	fout << ' '  <<	sy->lx-sy->lx_half << ' ' << y -sy->ly_half << ' ' << zz << endl;
}

