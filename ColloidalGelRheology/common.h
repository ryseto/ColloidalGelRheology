//
//  common.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_common_h
#define ColloidalGelRheology_common_h
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "vec3d.h"
using namespace std;

#define L_SIZE 5

inline double sq(double x){ return x*x;}

#define foreach(a,b,c) \
for(a::iterator c=(b).begin();(c)!=(b).end(); ++ (c))

#define ForAllBond \
for(vector<Bond *>::iterator bond_iter=(bond).begin();(bond_iter)!=(bond).end(); ++(bond_iter))

#define ForAllParticle \
for(vector<Particle *>::iterator p_iter=(particle).begin();(p_iter)!=(particle).end(); ++(p_iter))


#define fout_pos_c(X,Y,Z) \
fout << "c " << X << ' ' << Y << ' ' << Z  << endl;

#define fout_pos(X,Y,Z) \
fout << ' ' << X << ' ' << Y << ' ' << Z << ' ';

#define out_pos_c(X,Y,Z) \
out << "c " << X << ' ' << Y << ' ' << Z  << endl;

#define out_pos(X,Y,Z) \
out << ' ' << X << ' ' << Y  << ' ' << Z << ' ';

#define o_pos(O,X,Y,Z) \
O << ' ' << X << ' ' << Y << ' ' << Z << ' ';

#define cout_pos(X,Y,Z) \
cout << ' ' << X << ' ' << Y  << ' ' << Z << ' ';

#define cout_pos_c(X,Y,Z) \
cout << "c " << X << ' ' << Y << ' ' << Z << endl;


#define cout_pos_o(X,Y,Z) \
cout << "o " << X << ' ' << Y << ' ' << Z  << endl;


#define SIN fast_sin // @@@@@@@@@@@@@@@@@@@@@@@@@@@@
#define COS fast_cos
inline double fast_sin(double x){return x*(1.0 - 0.1666666667*x*x );}
inline double fast_cos(double x){return 1.0 - 0.5*x*x;}
inline double sqrt_1(double x){return 1.0 + (x - 1.0) * ( 0.625 - 0.125*x);}

struct ConnectPoint {
	int bond;
	int next;
	double tor_angle;
	vec3d u;
};

struct GridPoint {
	int x;
	int y;
	int z;
};

struct BondParameter {
	double kn;
	double ks;
	double kb;
	double kt;
	double kn3;
	
	double fnc;
	double fsc;
	double mbc;
	double mtc;
	double n_max;
	double s_max;
	double b_max;
	double t_max;
	
};


inline int ipow(int p, int q){
	int x = 1;
	for (int i=0; i < q; i++){
		x *= p;
	}
	return x;
}


#endif
