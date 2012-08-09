//
//  grid.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>

#include "grid.h"
#include <algorithm>
//extern vector<Particle *> particle;
#define DELETE(x) if(x){delete [] x; x = NULL;}

Grid::Grid(void){
	std::cerr << "grid object is created" << std::endl;
}

Grid::~Grid(void){
	for (int gx = 0; gx < gx_max; gx++){
		for (int gy = 0; gy < gy_max ; gy++){
			DELETE(vcell[gx][gy]);
		}
		DELETE(vcell[gx]);
	}
	DELETE(vcell);
}

void Grid::init (const int num_of_particle, const double lx_, const double ly_, const double lz_,
				 const double grid_size)
{
	numberOfParticle = num_of_particle;
	h = grid_size;
	gx_max = (int)( lx_/h );
#ifndef TWODIMENSION
	/* 3D */
	gy_max = (int)( ly_/h );
#else
	/* 2D */
	gy_max = 1;
#endif
	gz_max = (int)( lz_/h );
	vcell = new vector<int> ** [gx_max];
	neighbor_cell = new	vector<GridPoint> ** [gx_max];
	for (int gx = 0; gx < gx_max; gx++){
		vcell[gx] = new vector<int> * [gy_max];
		neighbor_cell[gx] = new	vector<GridPoint> * [gy_max];
		for (int gy = 0; gy < gy_max; gy++){
			vcell[gx][gy] = new vector<int> [gz_max];
			neighbor_cell[gx][gy] = new	vector<GridPoint> [gz_max];
		}
	}
	gl.resize(numberOfParticle);
#ifndef TWODIMENSION
	/* 3D */
	for (int gx = 0; gx < gx_max; gx++){
		for (int gy = 0; gy < gy_max; gy++){
			for (int gz = 0; gz < gz_max; gz++){
				////////////////////////////////////////
				GridPoint gp;
				for(gp.x = gx-1; gp.x <= gx+1; gp.x++){
					for(gp.y = gy-1; gp.y <= gy+1; gp.y++){
						for(gp.z = gz-1; gp.z <= gz+1; gp.z++){
							if ( gp.z >= 0 && gp.z < gz_max){
								GridPoint gp_tmp = gp;
								if (gp_tmp.x == -1)
									gp_tmp.x = gx_max - 1;
								else if (gp_tmp.x == gx_max)
									gp_tmp.x = 0;
								if (gp_tmp.y == -1)
									gp_tmp.y = gy_max - 1;
								else if (gp_tmp.y == gy_max)
									gp_tmp.y = 0;
								
								neighbor_cell[gx][gy][gz].push_back( gp_tmp );
							}
						}
					}
				}
			}
		}
	}
#else
	/* 2D */
	int gy = 0;
	for (int gx = 0; gx < gx_max; gx++){
		for (int gz = 0; gz < gz_max; gz++){
			////////////////////////////////////////
			GridPoint gp;
			for(gp.x = gx-1; gp.x <= gx+1; gp.x++){
				for(gp.z = gz-1; gp.z <= gz+1; gp.z++){
					if ( gp.z >= 0 && gp.z < gz_max){
						GridPoint gp_tmp = gp;
						if (gp_tmp.x == -1)
							gp_tmp.x = gx_max - 1;
						else if (gp_tmp.x == gx_max)
							gp_tmp.x = 0;
						gp_tmp.y = 0;
						neighbor_cell[gx][gy][gz].push_back( gp_tmp );
					}
				}
			}
		}
	}
#endif
}

void Grid::remake(vector<Particle *> &particle){
#ifndef TWODIMENSION
	for (int i = 0; i < numberOfParticle; i++ ){
		vcell[gl[i].x][gl[i].y][gl[i].z].clear();
	}
	for (int i = 0; i < numberOfParticle; i++ ){
		gl[i] = p_to_grid( *particle[i]->p_pos() );
		vcell[gl[i].x][gl[i].y][gl[i].z].push_back(i);
	}
#else
	for (int i = 0; i < numberOfParticle; i++ ){
		vcell[gl[i].x][0][gl[i].z].clear();
	}
	for (int i = 0; i < numberOfParticle; i++ ){
		gl[i] = p_to_grid( *particle[i]->p_pos() );
		vcell[gl[i].x][0][gl[i].z].push_back(i);
	}
#endif
}

void Grid::remake_with_walls(double zbot, double ztop, vector<Particle *> &particle){
	/*
	 * zbot and ztop are not position of wall.
	 * Particles locating at z < zbot and z > ztop are
	 * marked as near-wall particles.
	 */
	for (int i = 0; i < numberOfParticle; i++ ){
		vcell[gl[i].x][gl[i].y][gl[i].z].clear();
	}
	vcell_wall[0].clear();
	vcell_wall[1].clear();
	for (int i = 0; i < numberOfParticle; i++ ){
		gl[i] = p_to_grid( *particle[i]->p_pos() );
		vcell[gl[i].x][gl[i].y][gl[i].z].push_back(i);
		if (particle[i]->wall == false) {
			if ( particle[i]->p.z < zbot ) {
				vcell_wall[0].push_back(i);
			} else if ( particle[i]->p.z > ztop ) {
				vcell_wall[1].push_back(i);
			}
		}
	}
}

void Grid::remake_with_bottom(double zbot, vector<Particle *> &particle){
	for (int i = 0; i < numberOfParticle; i++ ){
		vcell[gl[i].x][gl[i].y][gl[i].z].clear();
	}
	vcell_wall[0].clear();
	for (int i = 0; i < numberOfParticle; i++ ){
		gl[i] = p_to_grid( *particle[i]->p_pos() );
		vcell[gl[i].x][gl[i].y][gl[i].z].push_back(i);
		if (particle[i]->wall ==false) {
			if ( particle[i]->p.z < zbot ) {
				vcell_wall[0].push_back(i);
			}
		}
	}
}


GridPoint Grid::p_to_grid(const vec3d &p){
	/* 3D */
	gp_tmp.x = (int)(p.x/h);
	if (gp_tmp.x == gx_max ){
		gp_tmp.x --;
	}
#ifndef TWODIMENSION
	gp_tmp.y = (int)(p.y/h);
	if (gp_tmp.y == gy_max ){
		gp_tmp.y --;
	}
#else
	gp_tmp.y = 0;
#endif
	gp_tmp.z = (int)(p.z/h);
	if (gp_tmp.z == gz_max ){
		gp_tmp.z --;
	}
	return gp_tmp;
}

void Grid::entry(vec3d &p, int i){
	GridPoint gp = p_to_grid(p);
#ifndef TWODIMENSION
	vcell[gp.x][gp.y][gp.z].push_back(i);
#else
	vcell[gp.x][0][gp.z].push_back(i);
#endif
	gl[i] = gp;
}


void Grid::get_neighbor_list(const vec3d &p, vector<int> &neighbor){
	//GridPoint gp;
	gp_tmp = p_to_grid(p);
	foreach(vector<GridPoint>, neighbor_cell[gp_tmp.x][gp_tmp.y][gp_tmp.z], gp){
		foreach(vector<int>, vcell[(*gp).x][(*gp).y][(*gp).z], iter){
			neighbor.push_back( *iter );
		}
	}
}

void Grid::get_neighbor_list_pointer(vec3d &p, vector< vector<int>* > &p_neighbor){
	p_neighbor.clear();
	gp_tmp = p_to_grid(p);
	foreach(vector<GridPoint>, neighbor_cell[gp_tmp.x][gp_tmp.y][gp_tmp.z], gp){
		p_neighbor.push_back( &(vcell[(*gp).x][(*gp).y][(*gp).z]));
	}
}

vector<int> * Grid::get_neighbor_list_pointer_top(){
	return &(vcell_wall[1]);
}

vector<int> * Grid::get_neighbor_list_pointer_bot(){
	return &(vcell_wall[0]);
}

bool Grid::near_boundary_xy(const vec3d &p){
#ifndef TWODIMENSION
	int gx = (int)(p.x/h); // gx = gx_max is possible
	if ( gx == 0 || gx >= gx_max - 1)
		return true;
	
	int gy = (int)(p.y/h); // gy = gy_max is possible
	if (gy == 0 || gy >=  gy_max - 1)
		return true;
	// gx = 0, 1 and gx = gx_max-2 and gx_max-1
#else
	int gx = (int)(p.x/h);
	if ( gx == 0 || gx >= gx_max - 1 ){
		return true;
	}
#endif
	return false;
}
