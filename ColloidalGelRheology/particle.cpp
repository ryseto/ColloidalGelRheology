//
//  particle.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "particle.h"

Particle::Particle(int particle_number_,
				   const vec3d &position,
				   const int _i_cluster,
				   System &sy_){
	sy = &sy_;
	setInitial(particle_number_);
	p = position;
	init_cluster = _i_cluster;
	if (init_cluster != 0)
		exit(1);
	setVelocityZero();
	orientation.set(1., 0., 0., 0.);
	resetForce();
	wall = false;
#ifdef TWODIMENSION
	// 2D
	cn = new ConnectPoint [6];
#else
	// 3D
	cn = new ConnectPoint [12];
#endif
	wall_connected = 0;
}

Particle::Particle(int particle_number_){
	setInitial(particle_number_);
}

Particle::~Particle(){
	delete [] cn;
}

void Particle::setInitial(int particle_number_){
	particle_number = particle_number_;
	cn_size = 0;
}

void  Particle::makeNeighbor(){
	vector< vector<int> *> neighbor_cells;
	sy->grid->get_neighbor_list_pointer( p, neighbor_cells);
	near_boundary = sy->grid->near_boundary_xy(p);
	neighbor.clear();
	foreach(vector< vector<int> *>, neighbor_cells, iter){
		foreach( vector<int>, *(*iter), i){
			if ( sy->ct->connect(particle_number, *i ) == false
				&& *i < particle_number ){
				neighbor.push_back( *i );
			}
		}
	}
}

void Particle::generateBond(){
	if (near_boundary){
		unsigned long n_neighbor = neighbor.size();
		for (int i=0; i < n_neighbor; i++) {
			p_pdcopy = sy->particle[neighbor[i]]->p;
			if( abs(p.x - p_pdcopy.x) > 10) {
				if (p.x > p_pdcopy.x)
					p_pdcopy.x += sy->lx;
				else
					p_pdcopy.x -= sy->lx;
			}
#ifndef TWODIMENSION
			// 3D
			if( abs(p.y - p_pdcopy.y) > 10. ){
				if (p.y > p_pdcopy.y )
					p_pdcopy.y += sy->ly;
				else
					p_pdcopy.y -= sy->ly;
			}
#endif
			
#ifdef TWODIMENSION
			// 2D
			bool contact_distance = (sq_dist_2d(p, p_pdcopy) <= sy->sq_dist_generate);
#else
			// 3D
			bool contact_distance = (sq_dist(p, p_pdcopy) <= sy->sq_dist_generate );
#endif
			if ( contact_distance ) {
				sy->bond.push_back(new Bond (neighbor[i], particle_number, *sy));
				neighbor[i] = neighbor.back();
				neighbor.pop_back();
				n_neighbor --;
				i--; // because neighbor[i] must be checked again.
			}
		}
	} else {
		unsigned long n_neighbor = neighbor.size();
		for (int i=0; i < n_neighbor; i++){
			if (sq_dist(p, sy->particle[neighbor[i]]->p) <= sy->sq_dist_generate ){
				sy->bond.push_back(new Bond (neighbor[i], particle_number, *sy));
				neighbor[i] = neighbor.back();
				neighbor.pop_back();
				n_neighbor --;
				i --; // because neighbor[i] must be checked again.
			}
		}
	}
}

void Particle::delConnectPoint(int bond_number){
	cn_size --;
	for (int i=0; i < cn_size ; i++){
		if ( cn[i].bond == bond_number ){
			cn[i] = cn[ cn_size ];
			sy->bond[ cn[i].bond ]->chPointer(i, particle_number);
			break;
		}
	}
	return;
}

void Particle::checkNearBoundary()
{
	near_boundary = sy->grid->near_boundary_xy(p);
}

void Particle::move_Euler()
{
	/*
	 * The effect of eta should be checked.
	 */
	if (true){
		a_velocity = force-sy->eta*velocity;
		a_omega = 2.5*(torque-sy->eta_rot*omega);
		velocity += a_velocity*sy->dt;
		omega += a_omega*sy->dt;
	} else {
		velocity = force*sy->dt;
		omega = (torque-sy->eta_rot*omega)*sy->dt;
	}
	
	p += velocity*sy->dt;
	d_rotation = omega*sy->dt;
	orientation.infinitesimalRotation(d_rotation);
	for (int i = 0; i < cn_size ; i++) {
		cn[i].u.rotateInfinitesimal(d_rotation);
#ifndef TWODIMENSION
		// 3D
		cn[i].tor_angle += dot(d_rotation, cn[i].u);
#endif
	}
	resetForce();
	if (near_boundary) {
		if (p.x < 0) {
			p.x += sy->lx;
		} else if (p.x > sy->lx) {
			p.x -= sy->lx;
		}
#ifndef TWODIMENSION
		// 3D
		if (p.y < 0) {
			p.y += sy->ly;
		} else if (p.y > sy->ly) {
			p.y -= sy->ly;
		}
#endif
	}
}

void Particle::output(ofstream &fout){
	fout << "y 9" << endl;
	fout << "@ 2" << endl;
	fout << "c " << p.x - sy->lx0 ;
	fout << ' '  << p.y - sy->ly0 ;
	fout << ' '  << p.z - sy->lz0 << endl;
	vec3d vec[6];
	vec[0].set( 0.5,   0,   0);
	vec[1].set( 0  , 0.5,   0);
	vec[2].set( 0  ,   0, 0.5);
	vec[3].set(-0.5,   0,   0);
	vec[4].set( 0  ,-0.5,   0);
	vec[5].set( 0  ,   0,-0.5);
	fout << "y 10" << endl;
	for(int i=0; i<6; i ++){
		vec[i] = orientation.ori_forward(vec[i]);
		fout << "l " << p.x - sy->lx0 ;
		fout << ' '  << p.y - sy->ly0 ;
		fout << ' '  << p.z - sy->lz0 ;
		fout << ' '  << p.x - sy->lx0 + vec[i].x;
		fout << ' '  << p.y - sy->ly0 + vec[i].y;
		fout << ' '  << p.z - sy->lz0 + vec[i].z;
		fout << endl;
	}
	fout << "y 11" << endl;
	fout << "@ 3" << endl;
	for ( int i =0; i < cn_size ; i++){
		fout << "l " << p.x - sy->lx0 ;
		fout << ' '  << p.y - sy->ly0 ;
		fout << ' '  << p.z - sy->lz0 ;
		fout << ' '  << p.x - sy->lx0 + cn[i].u.x;
		fout << ' '  << p.y - sy->ly0 + cn[i].u.y;
		fout << ' '  << p.z - sy->lz0 + cn[i].u.z;
		fout << endl;
	}
}

void Particle::calc_stack_Force(){
	for (int i = 0; i < cn_size ; i++){
		sy->bond[ cn[i].bond ]->calcForce();
		force += sy->bond[cn[i].bond ]->forceToParticle(particle_number);
	}
}

void Particle::setRotate(vec3d axis, const double angle){
	double tmp_angle = 0;
	double d_ang = 0.001;
	while(tmp_angle < angle){
		tmp_angle += d_ang;
#ifdef TWODIMENSION
		// 2D
		d_rotation = d_ang;
		orientation.infinitesimalRotation( d_rotation );
		for (int i = 0; i < cn_size ; i++){
			cn[i].u.rotateInfinitesimal( d_rotation);
		}
#else
		// 3D
		d_rotation = d_ang*axis;
		orientation.infinitesimalRotation( d_rotation );
		for (int i = 0; i < cn_size ; i++){
			cn[i].u.rotateInfinitesimal( d_rotation );
			cn[i].tor_angle += dot( d_rotation, cn[i].u );
		}
#endif
	}
}

bool Particle::percolate(vector<int> perco_path){
	if (wall == true){
		if (p.z > 5){
			return true;
		} else {
			return false;
		}
	}
	
	for (int i = 0; i < cn_size ; i++){
		bool passed = false;
		for (int j=0; j< perco_path.size(); j++){
			if ( perco_path[i] == cn[i].next){
				passed = false;
				break;
			}
		}
		if (passed == false){
			perco_path.push_back(particle_number);
			if (sy->particle[cn[i].next]->percolate(perco_path) == true){
				return true;
			}
		}
	}
	return false;
}

void Particle::markWallConnected(int wt, vector<int> &wall_group){
	// wt : 2 // top
	// wt : 1 // bot
	if ( wall_connected == -1 ){
		wall_group.push_back(particle_number);
		wall_connected = wt;
		for (int i = 0; i < cn_size ; i++){
			sy->particle[cn[i].next]->markWallConnected(wt, wall_group);
		}
	} else {
		if ( wall_connected + wt == 3){
			sy->percolation = true;
			return;
		}
	}
}

void Particle::x_shift( double dx ){
	p.x += dx;
	if (p.x > sy->lx){
		p.x -= sy->lx;
	} else if (p.x < 0){
		p.x += sy->lx;
	}
}

void Particle::outputBond(){
	cout << "@ 2" << endl;
	cout << "t " << p.x - sy->lx0 << ' ' << p.y- sy->ly0  << ' ' << p.z- sy->lz0 << ' ' << particle_number << endl;
	for (int i = 0; i < cn_size ; i++){
		cout << "l " << p.x - sy->lx0 << ' ' << p.y- sy->ly0  << ' ' << p.z- sy->lz0  << ' ';
		cout << p.x + cn[i].u.x - sy->lx0 << ' ' ;
		cout << p.y + cn[i].u.y - sy->ly0  << ' ';
		cout << p.z + cn[i].u.z - sy->lz0  << endl;
	}
	if (cn_size >= 5){
		cout << "r 1 " << endl;
		cout << "@ 3 " << endl;
		cout << "c ";
		cout << p.x - sy->lx0 << ' ' ;
		cout << p.y - sy->ly0  << ' ';
		cout << p.z - sy->lz0  << endl;
	}
}
