//
//  wall.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

//#include <iostream>

#include "wall.h"

Wall::Wall(int _objectID, double _z, System &sy_) : objectID (_objectID), z(_z){
	sy = &sy_;
	if ( objectID == sy->n_top){
		walltype = top;
		u.set(0., 0., -1.); 
	} else if ( objectID ==	sy->n_bot){
		walltype = bot;		
		u.set(0., 0., 1.); 	
	}
	velocity.reset();
	x = sy->lx0;
	y = sy->ly0;
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
    int n_neighbor = neighbor.size();
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
        int n_particle_active = particle_active.size();
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
    stress_sc = force_z / (sy->lx*sy->ly);
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
    stress_sc = force_x / (sy->lx*sy->ly);
    return stress_sc;
}

void Wall::getParticles(vector<int> &particles_list){    
    foreach( vector <Particle* >, wall_particle, iter_p){
        //wall_particle.push_back();
        particles_list.push_back((*iter_p)->particle_number);
        
    }
}

bool Wall::markWallConnected(){
    int wt;
    if (walltype == bot){
        wt = 1;
    } else {
        wt = 2;        
    }
    foreach( vector <Particle* >, wall_particle, iter_p){
        (*iter_p)->markWallConnected(wt, wall_group);
        if ( sy->percolation == true){
            return true;
        }
    }
    cerr << "wall" << wt << " : "  << wall_group.size() << endl;
    return false;
}


void Wall::output(ofstream &fout){
	double x_wall = x;
	while (x_wall > sy->lx){
		x_wall += - sy->lx;
	}
	while ( x_wall < 0 ){
		x_wall +=  sy->lx;
	} 
	
	double zz = z - sy->lz0;
	fout << "l " <<       -sy->lx0 << ' ' <<        -sy->ly0 << ' ' << zz;
	fout <<  ' ' << sy->lx-sy->lx0 << ' ' <<        -sy->ly0 << ' ' << zz << endl;
	fout << "l " << sy->lx-sy->lx0 << ' ' <<        -sy->ly0 << ' ' << zz ;
	fout <<  ' ' << sy->lx-sy->lx0 << ' ' << sy->ly -sy->ly0 << ' ' << zz << endl;
	fout << "l " << sy->lx-sy->lx0 << ' ' << sy->ly -sy->ly0 << ' ' << zz ;
	fout <<  ' ' <<       -sy->lx0 << ' ' << sy->ly -sy->ly0 << ' ' << zz << endl;
	fout << "l " <<       -sy->lx0 << ' ' << sy->ly -sy->ly0 << ' ' << zz ;
	fout <<  ' ' <<       -sy->lx0 << ' ' <<        -sy->ly0 << ' ' << zz << endl;
	fout << "l " << x_wall-sy->lx0 << ' ' <<        -sy->ly0 << ' ' << zz;
	fout << ' '  <<	x_wall-sy->lx0 << ' ' << sy->ly -sy->ly0 << ' ' << zz << endl;
	fout << "l " <<       -sy->lx0 << ' ' << y -sy->ly0 << ' ' << zz;
	fout << ' '  <<	sy->lx-sy->lx0 << ' ' << y -sy->ly0 << ' ' << zz << endl;
}

