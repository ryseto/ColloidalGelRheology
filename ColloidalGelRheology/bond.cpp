//
//  bond.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "bond.h"

Bond::Bond(const int d0, const int d1, System &sy_)
{
	sy = &sy_;
	status = 1;
	initial_bond = sy->initialprocess;
	cnt_regeneration = 0;
	bond_number = sy->n_bond++;
	d[0] = d0;
	d[1] = d1;
	sy->ct->on_connect( d[0], d[1] );
	p_particle0 = sy->particle[d[0]];
	p_particle1 = sy->particle[d[1]];
	if (initial_bond) {
		para = sy->bond0;
		bondtype = 0;
	} else {
		para = sy->bond1;
		bondtype = 1;
	}
	if (para.fsc == 0) {
		central_force = true;
	} else {
		central_force = false;
	}
	(*p_particle1).setBond(bond_number, d[0]);
	(*p_particle0).setBond(bond_number, d[1]);
	pp[0] = (*p_particle0).p_pos();
	pp[1] = (*p_particle1).p_pos();
	update_rvec();
#ifdef TWODIMENSION
	// 2D
	r = r_vec.norm_2d();
	e_normal = r_vec.division_2d(r);
#else
	// 3D
	r = r_vec.norm();
	e_normal = r_vec/r;
#endif
	if (r < 1.5) {
		cerr << "new bond: r = " << r << " < 1.5" << endl;
		pp[0]->cerr();
		pp[1]->cerr();
		exit(1);
	}
	pu[0] = (*p_particle0).pu_back();
	pu[1] = (*p_particle1).pu_back();
	(*pu[0]) = e_normal;
	(*pu[1]) = - e_normal;
#ifdef TWODIMENSION
	contactangle[0] = atan2( e_normal.z, e_normal.x)-(*p_particle0).orientation_angle;
	contactangle[1] = atan2( -e_normal.z, -e_normal.x)-(*p_particle1).orientation_angle;
#endif
	u_vector[0] = (*p_particle0).orientation.ori_backward(*pu[0]);
	u_vector[1] = (*p_particle1).orientation.ori_backward(*pu[1]);
	if (initial_bond) {
		u_vector_initial[0] = u_vector[0];
		u_vector_initial[1] = u_vector[1];
	}
#ifndef TWODIMENSION
	// 3D
	p_tor[0] = (*p_particle0).p_tor_angle_back();
	p_tor[1] = (*p_particle1).p_tor_angle_back();
	(*p_tor[0]) = 0.;
	(*p_tor[1]) = 0.;
#endif
#ifdef TWODIMENSION
	// 2D
	sq_fsc = sq(para.fsc);
#else
	// 3D
	this need to be written.
	sq_fsc = sq(para.fsc);
	sq_mbc = sq(para.mbc);
	sq_mtc = sq(para.mtc);
#endif
	force_normal = 0;
	force_sliding.reset();
#ifdef TWODIMENSION
	// 2D
	moment_bending = 0;
#else
	// 3D
	moment_bending.reset()
	moment_torsion = 0;
#endif
}

Bond::~Bond()
{
	cerr << "deleted bond" << endl;
}

void Bond::getParticleNumbers(int &i, int &j)
{
	i = (*p_particle0).particle_number;
	j = (*p_particle1).particle_number;
}

void Bond::addContactForce()
{
	if (status == 0) return;
	calcForce();
	(*p_particle0).stackForce(force0, torque0);
	(*p_particle1).stackForce(-force0, torque1);
}

void Bond::calcForce()
{
	update_rvec();
#ifdef TWODIMENSION
	// 2D
	r = r_vec.norm_2d();
	e_normal = r_vec.division_2d(r);
#else
	// 3D
	r = r_vec.norm();
	e_normal = r_vec/r;
#endif
	q = r-2;
	u01 = *pu[1]-*pu[0];
#ifdef TWODIMENSION
	// 2D
	d_slid = u01-dot_2d(u01, e_normal)*e_normal;
	ang_bend = -cross_2d(*pu[0], *pu[1]);
#else
	// 3D
	d_slid = u01-dot(u01, e_normal)*e_normal;
	ang_bend = -cross(*pu[0], *pu[1]);
	ang_tort = -(*p_tor[0])-(*p_tor[1]);
#endif
	force_normal = para.kn*q;
	force_sliding  = para.ks*d_slid;
	moment_bending = para.kb*ang_bend;
#ifndef TWODIMENSION
	// 3D
	moment_torsion = para.kt*ang_tort;
#endif
	force0 = force_normal*e_normal+force_sliding;
#ifdef TWODIMENSION
	// 2D
	moment_sliding = cross_2d(e_normal, force_sliding);
	torque0 = moment_sliding+moment_bending;
	torque1 = moment_sliding-moment_bending;
#else
	// 3D
	moment_sliding = cross(e_normal, force_sliding);
	torque_tmp = moment_bending+moment_torsion *e_normal;
	torque0 = moment_sliding+torque_tmp;
	torque1 = moment_sliding-torque_tmp;
#endif
	return;
}

void Bond::chPointer(int i, int particle_num)
{
	if (d[0] == particle_num) {
		pu[0] = (*p_particle0).pu(i);
#ifndef TWODIMENSION
		// 3D
		p_tor[0] = (*p_particle0).p_tor_angle(i);
#endif
	} else {
		pu[1] = (*p_particle1).pu(i);
#ifndef TWODIMENSION
		// 3D
		p_tor[1] = (*p_particle1).p_tor_angle(i);
#endif
	}
}

void Bond::update_rvec()
{
	r_vec = (*pp[1])-(*pp[0]);
	if (r_vec.x > sy->lx_half) {
		r_vec.x -= sy->lx;
	} else if (r_vec.x < -sy->lx_half) {
		r_vec.x += sy->lx;
	}
	if (r_vec.z > sy->lz_half) {
		r_vec.z -= sy->lz;
	} else if (r_vec.z < -sy->lz_half) {
		r_vec.z += sy->lz;
	}
#ifndef TWODIMENSION
	if (r_vec.y > sy->ly_half) {
		r_vec.y -= sy->ly;
	} else if (r_vec.y < -sy->ly_half) {
		r_vec.y += sy->ly;
	}
#endif
}

void Bond::rupture()
{
	status = 0;
	sy->ct->off_connect(d[0], d[1]);
	(*p_particle0).delConnectPoint(bond_number);
	(*p_particle1).delConnectPoint(bond_number);
	if (d[1] < d[0]) {
		(*p_particle0).addNeighbor(d[1]);
	} else {
		(*p_particle1).addNeighbor(d[0]);
	}
}

void Bond::regeneration()
{
	if (status == 1) {
		status = 2;
		para = sy->bond1;
		sq_fsc = sq(para.fsc);
		sq_mbc = sq(para.mbc);
		sq_mtc = sq(para.mtc);
		if (para.fnc == 0) {
			central_force = true;
		}
	}
	cnt_regeneration++;
	//  check_boundary();
	update_rvec();
#ifdef TWODIMENSION
	// 2D
	r = r_vec.norm_2d();
	e_normal = r_vec.division_2d(r);
#else
	// 3D
	r = r_vec.norm();
	e_normal = r_vec/r;
#endif
	(*pu[0]) = e_normal;
	(*pu[1]) = -e_normal;
#ifdef TWODIMENSION
	contactangle[0] = atan2(e_normal.z, e_normal.x)-(*p_particle0).orientation_angle;
	contactangle[1] = atan2(-e_normal.z, -e_normal.x)-(*p_particle1).orientation_angle;
#endif
#ifndef TWODIMENSION
	// 3D
	(*p_tor[0]) = 0;
	(*p_tor[1]) = 0;
#endif
	u_vector[0] = (*p_particle0).orientation.ori_backward(*pu[0]);
	u_vector[1] = (*p_particle1).orientation.ori_backward(*pu[1]);
	//////////////////////////////////////////////////////////////////////
	u01 = *pu[1] - *pu[0];
}

void Bond::cheackBondStress()
{
	if (status == 0) {
		return;
	}
	/* Condition of bond failure
	 * The bond can be robuster by compression?
	 * The normalization for des[0] should not be cp_f_n_max;
	 * Now, this effect is neglected.
	 */
#ifdef TWODIMENSION
	// 2D
	if (force_normal > para.fnc) {
		sy->rupture_bond.push_back(bond_number);
		sy->rup_normal ++;
	} else if (force_sliding.sq_norm() > sq_fsc) {
		sy->regeneration_bond.push_back(bond_number);
		sy->rup_shear ++;
	} else if ( abs(moment_bending) > para.mbc) {
		sy->regeneration_bond.push_back(bond_number);
		sy->rup_bend ++;
	}
#else
	// 3D
	not yet implemented
	if (force_normal > para.fnc) {
		sy->rupture_bond.push_back(bond_number);
		sy->rup_normal ++;
	} else if (force_sliding.sq_norm() > sq_fsc) {
		sy->regeneration_bond.push_back(bond_number);
		sy->rup_shear ++;
	} else if (abs(moment_bending) > para.mbc) {
		sy->regeneration_bond.push_back(bond_number);
		sy->rup_bend ++;
	} else if (.....) {
		sy->regeneration_bond.push_back(bond_number);
		sy->rup_torsion ++;
	}
#endif
	return;
}

void Bond::monitor_state(ofstream &out)
{
	/* To be checked.
	 *
	 */
	vec3d u_init[2];
	vec3d e_normal_init;
	u_init[0] = (*p_particle0).orientation.ori_forward(u_vector[0]);
	u_init[1] = (*p_particle1).orientation.ori_forward(u_vector[1]);
	update_rvec();
	double normal_distance;
#ifdef TWODIMENSION
	// 2D
	normal_distance = r_vec.norm_2d();
	e_normal_init = r_vec.division_2d(r);
#else
	// 3D
	normal_distance = r_vec.norm();
	e_normal_init = r_vec/normal_distance;
#endif
	u01 = u_init[1]-u_init[0];
	double slide_distance = (u01-dot(u01, e_normal_init)*e_normal_init).norm();
	double dot_prod = dot(u_init[0], -u_init[1]);
	double rolling_angle;
	if (dot_prod > 0.9) {
#ifdef TWODIMENSION
		// 2D
		double cross_prod = abs( cross_2d( u_init[0], - u_init[1]));
#else
		// 3D
		double cross_prod = (cross( u_init[0], - u_init[1])).norm();
#endif
		rolling_angle = asin(cross_prod);
	} else {
		rolling_angle = acos(dot_prod);
	}
	double torsional_angle = 0; 
	out << bond_number << ' ' << status  << ' ';
	out << normal_distance << ' ' << slide_distance << ' ' << rolling_angle << ' ' << torsional_angle << endl ;
}
