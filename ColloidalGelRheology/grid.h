//
//  grid.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_grid_h
#define ColloidalGelRheology_grid_h

#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "common.h"
#include "vec3d.h"
#include "particle.h"
const int max_number_point = 20000;
class Particle;

class Grid{
	double h;
	std::vector<GridPoint> gl;
	std::vector<int> ***vcell;
	std::vector<int> vcell_wall[2];
	std::vector<GridPoint> ***neighbor_cell;
	int numberOfParticle;
	GridPoint gp_tmp;
	int gx_max;
	int gy_max;
	int gz_max;
public:
	Grid(void);
	~Grid(void);
	
	void init(const int num_of_particle, const double lx_, const double ly_, const double lz_,
			  const double grid_size);
	void init_z(const int num_of_particle, const double lx_, const double ly_, const double lz_,
				 const double grid_size);
	void remake(std::vector<Particle *> &particle);
	void remake_with_walls(double, double, std::vector<Particle *> &particle);
	void remake_with_bottom(double, std::vector<Particle *> &particle);
	
	void entry(vec3d &p, int i);
	GridPoint p_to_grid(const vec3d &p);
#ifndef TWODIMENSION
	inline std::vector<int>* particle_in_cell(int x, int y, int z){ return &vcell[x][y][z];}
	inline std::vector<int>* particle_in_toplayer(int x, int y) { return &vcell[x][y][gz_max]; }
	inline std::vector<int>* particle_in_botlayer(int x, int y) { return &vcell[x][y][0]; }
#else
	inline std::vector<int>* particle_in_cell(int x, int y, int z){ return &vcell[x][0][z];}
	inline std::vector<int>* particle_in_toplayer(int x, int y) { return &vcell[x][0][gz_max]; }
	inline std::vector<int>* particle_in_botlayer(int x, int y) { return &vcell[x][0][0]; }
#endif
	void reset();
	inline int val_gx_max(){ return gx_max; }
	inline int val_gy_max(){ return gy_max; }
	inline int val_gz_max(){ return gz_max; }
	inline int min_gy(int gy){ return ( gy >= 0 ? gy : 0 ) ; }
	inline int max_gy(int gy){ return ( gy <= gy_max ? gy : gy_max ); }
	int gy(int i){ return gl[i].y;}
	unsigned long size(int x, int y, int z) { return vcell[x][y][z].size(); }
	void gl_resize(int n){ gl.resize(n); }
	void get_neighbor_list(const vec3d &p, std::vector<int> &neighbor);
	void get_neighbor_list_pointer(vec3d &p, std::vector<std::vector<int>* > &p_neighbor);
	std::vector<int> * get_neighbor_list_pointer_top();
	std::vector<int> * get_neighbor_list_pointer_bot();
	bool near_boundary_xy(const vec3d &p);
	
};

#endif
