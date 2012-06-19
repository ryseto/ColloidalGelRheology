//
//  bond.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>


#include "bond.h"

Bond::Bond(int d0, int d1, System &sy_){
	sy = &sy_;	
    status = 1;    
    initial_bond = sy->initialprocess;
    cnt_regeneration = 0;
	bond_number = sy->n_bond ++; 
	d[0] = d0, d[1] = d1;  
	sy->ct->on_connect( d[0], d[1] );
	p_particle0 = sy->particle[d[0]];
	p_particle1 = sy->particle[d[1]];
    if (initial_bond)
        para = sy->bond0;
    else
        para = sy->bond1;
    
    if (para.fnc == 0){
        central_force = true;
    } else {
        central_force = false;
    }
    (*p_particle1).setBond(bond_number, d[0]);
    (*p_particle0).setBond(bond_number, d[1]);
	
    pp[0] = (*p_particle0).p_pos();
	pp[1] = (*p_particle1).p_pos();
    
    periodicBoundary_rvec(r_vec);    
#ifndef TWODIMENSION
    r = r_vec.norm();
	e_normal = r_vec/r;
#else
	// 2D
	r = r_vec.norm_2d();
	e_normal = r_vec.division_2d(r);
#endif
    if (r < 1.5){
        cerr << "new bond: r = " << r << " < 1.5" << endl;
        pp[0]->cerr();
        pp[1]->cerr();
        exit(1);
    }  
    pu[0] = (*p_particle0).pu_back();
    pu[1] = (*p_particle1).pu_back();
    (*pu[0]) = e_normal;
    (*pu[1]) = - e_normal;
    u_vector[0] = (*p_particle0).orientation.ori_backward(*pu[0]);
    u_vector[1] = (*p_particle1).orientation.ori_backward(*pu[1]);
    if (initial_bond){
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
    sq_fsc = sq(para.fsc);
    sq_mbc = sq(para.mbc);
    sq_mtc = sq(para.mtc);
    
    q_old = 0;
    d_slid_old.reset();
    ang_bend_old.reset();
    ang_tort_old = 0;
}

Bond::~Bond(){
    cerr << "deleted bond" << endl;
}

void Bond::whichparticle(int &i, int &j){
    i = (*p_particle0).particle_number;
    j = (*p_particle1).particle_number;
}

void Bond::addContactForce(){
	calcForce();
	(*p_particle0).stackForce( force0, torque0 );
	(*p_particle1).stackForce(-force0, torque1 );
}

void Bond::calcForce(){
	periodicBoundary_rvec(r_vec);
#ifndef TWODIMENSION
	/* 3D */
	r = r_vec.norm();
	e_normal = r_vec/r;
#else
    /* 2D */
	r = r_vec.norm_2d();
	e_normal = r_vec.division_2d(r);
#endif
    q = r - 2.;    
    if ( central_force ){
        if ( q < 0.){
            /* compression --> negative */
            kn_compression = para.kn + para.kn3*q*q;
            force_normal = kn_compression*q + (2.*sqrt(kn_compression)/sy->dt)*(q-q_old);
        } else {
            force_normal = 0;
        }
        force0 = force_normal*e_normal;
        torque0.reset();
        torque1.reset();
        q_old = q;
        return;
    }        
    /* 3D and 2D */
	u01 = *pu[1] - *pu[0];
#ifndef TWODIMENSION
    /* 3D */
    d_slid = u01 - dot(u01, e_normal)*e_normal;
    ang_bend = - cross(*pu[0], *pu[1]);
    ang_tort = - (*p_tor[0]) - (*p_tor[1]);
#else
    /* 2D */
	d_slid = u01 - dot_2d(u01, e_normal)*e_normal;
    ang_bend = - cross_2d(*pu[0], *pu[1]);    
#endif
	if (q < 0.){
		/* compression --> negative */
        kn_compression = para.kn + para.kn3*q*q;
		//force_normal = kn_compression*q + (2.*sqrt(kn_compression)/sy->dt)*(q-q_old);
        force_normal = kn_compression*q + (para.c_norm/sy->dt)*(q-q_old);
	} else {
		/* traction --> positive */
		force_normal =  para.kn*q + (para.c_norm/sy->dt)*(q-q_old);
	}
    /* 3D and 2D */
    force_sliding = para.ks*d_slid + (para.c_slid/sy->dt)*(d_slid-d_slid_old);
    moment_bending = para.kb*ang_bend + (para.c_bend/sy->dt)*(ang_bend-ang_bend_old);
#ifndef TWODIMENSION
    /* 3D */
	moment_torsion = para.kt*ang_tort + (para.c_tort/sy->dt)*(ang_tort-ang_tort_old);
#endif
    /* 3D and 2D */
	force0 = force_normal*e_normal + force_sliding;
#ifndef TWODIMENSION
    /* 3D */
    moment_sliding = cross(e_normal, force_sliding);
    torque_tmp = moment_bending + moment_torsion *e_normal;
    torque0 = moment_sliding + torque_tmp;
    torque1 = moment_sliding - torque_tmp;
#else
    /* 2D */
    moment_sliding = cross_2d(e_normal, force_sliding);
    //torque_tmp = moment_bending;
    torque0 = moment_sliding + moment_bending;
    torque1 = moment_sliding - moment_bending;
#endif    
    /* This is for dampers*/
    q_old = q;
    d_slid_old = d_slid;
    ang_bend_old = ang_bend;
#ifndef TWODIMENSION
    /* 3D */
    ang_tort_old = ang_tort;
#endif
	return;
}

void Bond::chPointer(int i, int particle_num){
	if (d[0] == particle_num){
		pu[0] = (*p_particle0).pu(i);
#ifndef TWODIMENSION
		p_tor[0] = (*p_particle0).p_tor_angle(i);
#endif
	} else {
		pu[1] = (*p_particle1).pu(i);
#ifndef TWODIMENSION
		p_tor[1] = (*p_particle1).p_tor_angle(i);
#endif
	}
}

void Bond::periodicBoundary_rvec(vec3d & r_vec_tmp){
	r_vec_tmp = (*pp[1]) - (*pp[0]);
	// use boundary flag for particles and bonds.
	if (p_particle0->near_boundary || p_particle1->near_boundary){
        if (abs(r_vec_tmp.x) > 5.){
            if (r_vec_tmp.x > 0.) 
                r_vec_tmp.x -= sy->lx;
            else 
                r_vec_tmp.x += sy->lx;
        }
#ifndef TWODIMENSION
        if (abs(r_vec_tmp.y) > 5.){
            if (r_vec_tmp.y > 0.) 
                r_vec_tmp.y -= sy->ly;
            else 
                r_vec_tmp.y += sy->ly;
        }
#endif
    }
}

void Bond::rupture(){
    /*
     cout << "@ 5" << endl;
     cout << "r 0.3" << endl;
     cout << "c " << (*p_particle0).p.x - sy->lx0 << ' ';
     cout <<         (*p_particle0).p.y - sy->ly0 << ' ';
     cout <<         (*p_particle0).p.z - sy->lz0 << endl;
     cout << "c " << (*p_particle1).p.x - sy->lx0 << ' ';
     cout <<         (*p_particle1).p.y - sy->ly0 << ' ';
     cout <<         (*p_particle1).p.z - sy->lz0 << endl;
     */
    status = 0;
	sy->ct->off_connect( d[0], d[1] );
    (*p_particle0).delConnectPoint(bond_number);
    (*p_particle1).delConnectPoint(bond_number);
	if ( d[1] < d[0] ){
		(*p_particle0).addNeighbor( d[1] );
    }else{
		(*p_particle1).addNeighbor( d[0] );
    }
}

void Bond::outputRuptureMark(){
    /* Output the rupture event for debug.
     *
     */
     cerr << (*p_particle0).particle_number << "--x--" << (*p_particle1).particle_number << endl;
     cerr << des[0] << ' ' << des[1] << ' ' << des[2] << endl;
     cerr << force_normal <<' ' << force_sliding.norm() << ' ' << moment_bending.norm()  << endl;
     cout << "@ 4" << endl;
     cout << "r 0.3" << endl;
     cout << "c " << (*p_particle0).p.x - sy->lx0 << ' ';
     cout <<         (*p_particle0).p.y - sy->ly0 << ' ';
     cout <<         (*p_particle0).p.z - sy->lz0 << endl;
     cout << "c " << (*p_particle1).p.x - sy->lx0 << ' ';
     cout <<         (*p_particle1).p.y - sy->ly0 << ' ';
     cout <<         (*p_particle1).p.z - sy->lz0 << endl;
}

void Bond::regeneration(){
    if (status == 1){
        status = 2;        
        para = sy->bond1;
        sq_fsc = sq(para.fsc);
        sq_mbc = sq(para.mbc);
        sq_mtc = sq(para.mtc);
        if (para.fnc == 0){
            central_force = true;
        }
    }
    cnt_regeneration++;
    //  check_boundary();
    periodicBoundary_rvec(r_vec);
#ifndef TWODIMENSION
	// 3D
	r = r_vec.norm();
	e_normal = r_vec/r;
#else
	// 2D
	r = r_vec.norm_2d();
	e_normal = r_vec.division_2d(r);
#endif
    (*pu[0]) = e_normal;
    (*pu[1]) = - e_normal;
#ifndef TWODIMENSION
    (*p_tor[0]) = 0.;
    (*p_tor[1]) = 0.;    
#endif
    u_vector[0] = (*p_particle0).orientation.ori_backward(*pu[0]);
    u_vector[1] = (*p_particle1).orientation.ori_backward(*pu[1]);    
//////////////////////////////////////////////////////////////////////
    u01 = *pu[1] - *pu[0];
#ifndef TWODIMENSION
    /* 3D */
    d_slid_old.reset();
    ang_bend_old.reset();
    ang_tort_old = 0;
#else
    /* 2D */
	d_slid_old.reset();
    ang_bend_old.reset();
#endif
    D_function = 0;
}

void Bond::cheackBondStress(){
	/* Condition of bond failure
     * The bond can be robuster by compression?
     * The normalization for des[0] should not be cp_f_n_max;
     * Now, this effect is neglected.
     */
    if (central_force){
        if ( q > 0.1 ){
            sy->rupture_bond.push_back(bond_number);
        }
        return;
    }
    if ( q < 0 ){
        // force_normal < 0;
        des[0] = force_normal * para.reinforce_factor;
    } else {
        if (para.fnc == 0){
            des[0] = 1;
        } else {
            des[0] = sq(force_normal/para.fnc);
        }
    }
	des[1] = force_sliding.sq_norm()/sq_fsc;
	des[2] = moment_bending.sq_norm()/sq_mbc;
#ifndef TWODIMENSION
	des[3] = sq(moment_torsion)/sq_mtc;
    D_function = des[0] + des[1] + des[2] + des[3] ;
#else
    D_function = des[0] + des[1] + des[2];
#endif
	if ( D_function >= 1 ){	
        if ( q >= 0 ){
            sy->rupture_bond.push_back(bond_number);                            
            sy->rup_normal ++;
        } else {            
            sy->regeneration_bond.push_back(bond_number);
#ifndef TWODIMENSION
            if ( des[1] > des[2] ){
                if (des[1] > des[3]){
                    sy->rup_shear ++; // 1 > 2 && 1 > 3
                } else {
                    sy->rup_torsion ++; // 3 > 1 > 2 
                }
            } else {
                if (des[2] > des[3]){
                    sy->rup_bend ++; // 2 > 1 && 2 > 3
                } else {                    
                    sy->rup_torsion ++;// 3 > 2 > 1 
                } 
            }
#else
            if (des[1] > des[2]){
                sy->rup_shear ++;                 
            } else {
                sy->rup_bend ++;
            }
#endif
        }
    }
    return;
}

void Bond::monitor_state(ofstream &out){
    vec3d u_init[2];
    vec3d r_vec_init;
    vec3d e_normal_init;
    u_init[0] = (*p_particle0).orientation.ori_forward(u_vector[0]);
    u_init[1] = (*p_particle1).orientation.ori_forward(u_vector[1]);
    //r_vec = *pp[1] - *pp[0];
	periodicBoundary_rvec(r_vec_init);
    double normal_distance;
#ifndef TWODIMENSION
	// 3D
    normal_distance = r_vec_init.norm();
	e_normal_init = r_vec_init/normal_distance;
#else
	// 2D
	normal_distance = r_vec_init.norm_2d();
	e_normal_init = r_vec_init.division_2d(r);
#endif
	u01 = u_init[1] - u_init[0];
    double slide_distance = (u01 - dot(u01, e_normal_init)*e_normal_init).norm();
    double dot_prod = dot( u_init[0], - u_init[1]);
    double rolling_angle;
    if ( dot_prod > 0.9 ){
        double cross_prod = (cross( u_init[0], - u_init[1])).norm();
        rolling_angle = asin(cross_prod);
    }else{
        rolling_angle = acos(dot_prod);
    }
    double torsional_angle = 0; 
    out << bond_number << ' ' << status  << ' ';
    out << normal_distance << ' ' << slide_distance << ' ' << rolling_angle << ' ' << torsional_angle << endl ;
}
