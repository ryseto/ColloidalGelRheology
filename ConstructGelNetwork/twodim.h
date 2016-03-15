//
//  twodim.h
//  ColloidalGelRheology
//
//  Created by Ryohei Seto on 2012/08/17.
//
//

#ifndef ColloidalGelRheology_twodim_h
#define ColloidalGelRheology_twodim_h
#include "tools.h"
#include "vec2d.h"
#include "Floc.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
extern double r;
extern int n;
extern double Lx;
extern double Lz;
extern double Lz_pd;
double sq_max_cluster_distance;
double max_cluster_size;

int ipow(int p, int q);
void twodim_prepare_fout(int filling, ofstream &fout,
						 int rank, int rsd, double phi);
double twodim_collision_value(vec2d const &po, vec2d const &p, double theta);
void twodim_make_fractal_object(int m, vector<vec2d> &p);
double twodim_collision_height(vec2d po, vec2d p);
void twodim_put_first_fo(vector<vec2d> &p, int fo_size);

void twodim_drop_fo(vector<vec2d> &p, int i_premier, int fo_size);
void twodim_fo_random_fall( vector<vec2d> &p, int fo_size);
void twodim_set_initial_pairs(vector<vec2d> &p);
double twodim_calc_y_max(vector<vec2d> &p);
void twodim_output_result(ofstream &fout, vector<vec2d> &p,
						  int number_of_floc_for_vf, int fo_size);
int twodim_fo_random_put( vector<vec2d> &p, int fo_size, int total_number);
int twodim_put_fo(vector<vec2d> &p, int i_first, int fo_size);

void monitor(Floc *floc,
			 int floc_size, int number_of_floc_for_vf)
{
	cout << "@ 0 " << endl;
	for (int i = 0; i < number_of_floc_for_vf; i++){
		for (int k = 0; k < floc_size ; k++){
			vec2d p = floc[i].pos(k);
			p.shift(0, 1);
			cout << "c " << p.x() - Lx/2<< ' ' << 0 <<  ' ' <<  p.z() - Lz/2 << endl;
			if ( p.x() < 1){
				cout << "@ 3 " << endl;
				cout << "c " << p.x() + Lx - Lx/2<< ' ' << 0 <<  ' ' <<  p.z() - Lz/2 << endl;
				cout << "@ 0 " << endl;
			} else if ( p.x() > Lx - 1){
				cout << "@ 3 " << endl;
				cout << "c " << p.x() - Lx - Lx/2<< ' ' << 0 <<  ' ' <<  p.z() - Lz/2 << endl;
				cout << "@ 0 " << endl;
			}
			if ( p.z() < 2){
				cout << "@ 3 " << endl;
				cout << "c " << p.x() - Lx/2<< ' ' << 0 <<  ' ' <<  p.z() + Lz_pd - Lz/2 << endl;
				cout << "@ 0 " << endl;
			} else if ( p.z() > Lz_pd - 1){
				cout << "@ 3 " << endl;
				cout << "c " << p.x() - Lx/2<< ' ' << 0 <<  ' ' <<  p.z() - Lz_pd - Lz/2 << endl;
				cout << "@ 0 " << endl;
			}
		}
	}
	
    cout << "l " << -Lx/2 << ' ' << 0 << ' ' << -Lz/2 << ' ';
    cout << Lx/2 << ' ' << 0 << ' ' << -Lz/2 << endl;
    cout << "l " << -Lx/2 << ' ' << 0 << ' ' << Lz/2 << ' ';
    cout << Lx/2 << ' ' << 0 << ' ' << Lz/2 << endl;
	//	cerr << "N = " << floc_size*number_of_floc_for_vf << endl;
	cout << endl;
}

void twodim_output_flocs(ofstream &fout,
						 Floc *floc,
						 int number_of_floc_for_vf,
						 int floc_size)
{
	for (int i=0; i<number_of_floc_for_vf; i++) {
		for (int k=0; k<floc_size ; k++) {
			vec2d p1 = floc[i].pos(k);
			for (int j=0; j<number_of_floc_for_vf; j++) {
				for (int l=0; l<floc_size ; l++) {
					if (i != j || k != l) {
						vec2d p2 = floc[j].pos(l);
						vec2d dp = p2-p1;
						if (sq_dist_pd(p1, p2) < 3.99) {
							cerr << i << ": " << k << endl;
							cerr << j << ": " << l << endl;
							cerr <<  dp.sq_norm() << endl;
							dp.cerr();
							exit(1);
						}
					}
				}
			}
		}
	}
	
	for (int i=0; i<number_of_floc_for_vf; i++) {
		for (int k=0; k<floc_size ; k++) {
			vec2d p = floc[i].pos(k);
			//p[k].shift(0, -y_max/2);
			//p.shift(0, 1);
			//p[i + j*fo_size].shift(-L/2, 0);
			//p[i + j*fo_size].cout();
			//p[i + j*fo_size].out_data();
			fout << setprecision(15);
			fout << p.x() << ' ' <<  p.z() << ' ' << 0 << endl;
		}
	}
}

bool checkCollision(Floc &floc1, Floc &floc2)
{
	for (int i = 0; i < floc1.size(); i++){
		for (int j = 0; j < floc2.size(); j++){
			if (sq_dist_pd(floc1.pos(i), floc2.pos(j) ) < 4.05)
				return true;
		}
	}
	return false;
}

bool checkCollision(Floc *floc,
					int f1, int f2,
					vector <int> &f1_col,
					vector <int> &f2_col)
{
	unsigned long floc_size = floc[0].size();
	for (int i1_connect = 0; i1_connect < floc[f1].connect.size(); i1_connect++){
		int i1_floc = floc[f1].connect[i1_connect];
		for (int i=0; i<floc_size; i++) {
			vec2d p1 = floc[i1_floc].pos(i);
			for (int i2_connect = 0; i2_connect < floc[f2].connect.size(); i2_connect++) {
				int i2_floc = floc[f2].connect[i2_connect];
				for (int j=0; j<floc_size; j++) {
					vec2d p2 = floc[i2_floc].pos(j);
					if (sq_dist_pd(p1, p2) <= 4.3) { // @@@ why 4.3?
						f1_col.push_back(i1_floc);
						f2_col.push_back(i2_floc);
					}
				}
			}
		}
	}
	return false;
}

void solveCollision(Floc *floc,
					int moving_floc,
					int collided_floc,
					vec2d direction,
					double &lamda_min,
					int &col_f1,
					int &col_f2)
{
	unsigned long floc_size = floc[0].size();
		
	for (int k=0; k < floc[moving_floc].connect.size(); k++) {
		int ii = floc[moving_floc].connect[k];
		for (int m=0; m < floc_size ; m++) {
			vec2d p1 = floc[ii].pos(m);
			for (int kk = 0; kk < floc[collided_floc].connect.size(); kk ++) {
				int jj = floc[collided_floc].connect[kk];
				for (int mm=0; mm<floc_size ; mm++) {
					vec2d dp = diff_pd(p1, floc[jj].pos(mm));
					double dpdv = dot(dp, direction);
					double disc = dpdv*dpdv - (dp.sq_norm()-4);
					if (disc > 0) {
						double sq_disc = sqrt(disc);
						double lamda1 = -dpdv-sq_disc;
						double lamda2 = -dpdv+sq_disc;
						if (lamda1 > 0
							&& lamda_min > lamda1) {
							lamda_min = lamda1;
							if (floc[ii].parent != -1) {
								col_f1 = floc[ii].parent;
							} else {
								col_f1 = ii;
							}
							if (floc[jj].parent != -1) {
								col_f2 = floc[jj].parent;
							} else {
								col_f2 = jj;
							}
						}
						if (lamda2 > 0
							&& lamda_min > lamda2) {
							lamda_min = lamda2;
							if (floc[ii].parent != -1) {
								col_f1 = floc[ii].parent;
							} else {
								col_f1 = ii;
							}
							if (floc[jj].parent != -1) {
								col_f2 = floc[jj].parent;
							} else {
								col_f2 = jj;
							}
						}
					}
				}
			}
		}
	}
}

void moveRandom(Floc *floc, int floc_size, int num_of_floc)
{
	vector<int> cluster_list;
	for (int j = 0; j < num_of_floc; j++) {
		floc[j].parent = -1;
		cluster_list.push_back(j);
	}
	vector <int> f1_col;
	vector <int> f2_col;
	while (true) {
		int num_cluster = (int)cluster_list.size();
		int i = cluster_list[lrand48() % num_cluster];
		double theta = 2*M_PI*drand48();
		vec2d direction(cos(theta),sin(theta));
		for (int k=0; k < floc[i].connect.size(); k++) {
			int ii = floc[i].connect[k];
			floc[ii].move(direction);
		}
		for (int jj = 0; jj < num_cluster; jj++) {
			int j = cluster_list[jj];
			if (i != j) {
				checkCollision(floc, i, j,
							   f1_col, f2_col);
			}
		}

		if (f1_col.size() > 0) {
			int col_f1, col_f2;
			for (int k=0; k < floc[i].connect.size(); k++) {
				int ii = floc[i].connect[k];
				floc[ii].move(-direction);
				cerr << ii  << ":" << ' ';
			}
			cerr << endl;
			double lamda_min = 1;
			for (int ii=0	; ii < f1_col.size() ; ii++) {
				solveCollision(floc, f1_col[ii] , f2_col[ii] ,
							   direction, lamda_min,
							   col_f1, col_f2 );
			}
			f1_col.clear();
			f2_col.clear();

			for (int k=0; k < floc[i].connect.size(); k++) {
				int ii = floc[i].connect[k];
				floc[ii].move(lamda_min*direction);
			}
			if (lamda_min < 1) {
				if (col_f1 < col_f2) {
					for (int k = 0; k < floc[col_f2].connect.size(); k++) {
						floc[col_f1].connect.push_back(floc[col_f2].connect[k]);
						floc[floc[col_f2].connect[k]].parent = col_f1;
					}
					cluster_list.erase(remove(cluster_list.begin(), cluster_list.end(), col_f2),
									   cluster_list.end());
				} else {
					for (int k = 0; k < floc[col_f1].connect.size(); k++) {
						floc[col_f2].connect.push_back(floc[col_f1].connect[k]);
						floc[floc[col_f1].connect[k]].parent = col_f2;
					}
					cluster_list.erase(remove(cluster_list.begin(), cluster_list.end(), col_f1),
									   cluster_list.end());
				}
				cerr <<  "num_cluster = " << num_cluster << endl;
			}
		}
		if (floc[0].connect.size() == num_of_floc) {
			break;
		}
	}
}

void putRandom(Floc *floc, int floc_size, int num_of_floc)
{
	for (int i = 0; i < num_of_floc; i++) {
		bool collide;
		int cnt = 0;
		do {
			collide = false;
			vec2d p_center_trial(Lx*drand48(), Lz_pd*drand48());
			floc[i].setCenterPosition(p_center_trial);
			for (int j = 0; j < i ; j ++) {
				if (checkCollision(floc[i], floc[j])) {
					collide = true;
					break;
				}
			}
			cnt ++;
		} while(collide);
	}
	for (int i = 0; i < num_of_floc; i++) {
		floc[i].connect.push_back(i);
	}
}

void buildSpaceFillingNetwork2D(int filling, int rank, int rsd, double vol_frac)
{
	r = 1;
	int k = 16; // number of particle 2^k
	srand48(rsd);
	int floc_size = ipow(2, rank);
	cerr << "fractal size = " << floc_size << endl;
	int number_of_particle_for_vf = vol_frac*Lx*Lz / M_PI;
	double dbl_number_of_floc_for_vf = (double)number_of_particle_for_vf / floc_size;
	int number_of_floc_for_vf;
	cerr << "number_of_particle_for_vf = " << number_of_particle_for_vf << endl;
	if (dbl_number_of_floc_for_vf - (int)dbl_number_of_floc_for_vf > 0.5) {
		number_of_floc_for_vf = (int)dbl_number_of_floc_for_vf  + 1;
	} else {
		number_of_floc_for_vf = (int)dbl_number_of_floc_for_vf ;
	}
	
	n = ipow(2, k);
	cerr << n << endl;
	if (n < number_of_floc_for_vf*floc_size) {
		cerr << "n= " << n << endl;
		cerr << "number_of_floc_for_vf = " << number_of_floc_for_vf << endl;
		cerr << " k should be larger" << endl;
		exit(1);
	}
	cerr << "n = " << n << endl;
	vector<vec2d> p;
	p.resize(n);
	// Making n/2 pairs with random angles
	twodim_set_initial_pairs(p);
	// Build cluster with the rank
	if (rank >= 2) {
		for (int i=2; i <= rank; i++) {
			int floc_size_tmp = ipow(2, i);
			cerr << floc_size_tmp << endl;
			twodim_make_fractal_object(floc_size_tmp, p);
		}
	} else if (rank == 1) {
		int fo_size_tmp = 2;
		twodim_make_fractal_object(fo_size_tmp, p);
	}
	
	Floc *floc;
	floc = new Floc [number_of_floc_for_vf ];
	double sq_r_max = 0;
	for (int i=0; i < number_of_floc_for_vf ; i++) {
		vector < vec2d > floc_pos;
		for (int j = 0; j < floc_size; j++) {
			floc_pos.push_back(p[j+floc_size*i]);
			if (floc_pos.back().sq_norm() > sq_r_max) {
				sq_r_max = floc_pos.back().sq_norm() ;
			}
		}
		floc[i].setSystemSize( Lx, Lz_pd);
		floc[i].setPosition(floc_pos);
	}
	cerr << sq_r_max << endl;
	max_cluster_size = sqrt(sq_r_max);
	sq_max_cluster_distance = 4*sq_r_max;
	//	int total_number;

	ofstream fout;
    if (filling == 1) {
		//fo_random_fall(p, fo_size);
    } else if (filling == 2) {
		//n = twodim_fo_random_put(p, floc_size, number_of_floc_for_vf);
		putRandom(floc, floc_size, number_of_floc_for_vf);
		//twodim_prepare_fout(filling, fout, rank, rsd, vol_frac, "d");
		twodim_output_result(fout, p, number_of_floc_for_vf, floc_size);
		
		twodim_output_flocs(fout, floc,
							number_of_floc_for_vf,
							floc_size);
		fout.close();
		moveRandom(floc, floc_size, number_of_floc_for_vf);
	}
	
	twodim_prepare_fout(filling, fout, rank, rsd, vol_frac);
	//	twodim_output_result(fout, p, number_of_floc_for_vf, floc_size);
	twodim_output_flocs(fout, floc,
						number_of_floc_for_vf,
						floc_size);
	fout.close();
	
	return ;
}

void twodim_prepare_fout(int filling, ofstream &fout,
						 int rank, int rsd, double phi)
{
	char name_string[64];
    int floc_size = ipow(2, rank);
	sprintf(name_string, "W%dH%d_cca%d_phi%1.2f_rand%d",
			(int)(Lx), (int)(Lz), floc_size, phi, rsd);
    string fn_fout;
    if (filling == 1) {
        fn_fout = "sediment_2D_" + (string)name_string + ".dat";
    } else if (filling == 2) {
        fn_fout = "rw_2D_" + (string)name_string + ".dat"; //cluster random walk
    }
	fout.open(fn_fout.c_str());
}

/////////////////////////////////////////////////////////////////
double twodim_collision_value(vec2d const &po, vec2d const &p,
							  double theta)
{
	double l1, l2;
	double A, B;
	A = (p.x()-po.x())*cos(theta)+(p.z()-po.z())*sin(theta);
	//B = sq_dist_pd(p, po) - 4*r*r;
	B = sq_dist(p, po) - 4*r*r;
	if (A*A-B < 0) {
		return 0;
	}
	l1 = - A - sqrt( A*A-B );
	l2 = - A + sqrt( A*A-B );
    
	if (l1 > l2) {//@@@@@@@@@@ l1 > l2
		return l1;
	} else {
		return l2;
	}
}

void twodim_make_fractal_object(int m, vector<vec2d> &p)
{
	// fo_size_tmp, p
	// m means the size of fractal cluster
	for (int k=0; k < n; k+= m) {
		double theta = 2*M_PI*drand48();
		double l;
		
		double l_max = 0;
		for (int i = k; i < k+m/2; i++){ //fixed particles
			for (int j = k+m/2; j < k+m; j++){ // particles shifted
				l = twodim_collision_value(p[i], p[j], theta);
				if (l >= l_max){
					l_max = l;
				}
			}
		}
		for (int j = k+m/2; j < k+m; j++) {
			p[j] = p[j].plus(l_max*cos(theta), l_max*sin(theta));
		}
		
		vec2d cm(0, 0); // center of mass
		for(int i=k; i<k+m; i++) {
			cm += p[i]/m;
		}
		
		for(int i=k; i<k+m; i++) {
			p[i] = p[i]-cm;
		}
	}
}

double twodim_collision_height(vec2d po, vec2d p)
{
	double h, A, B;
	A = p.z()-po.z();
	//B = sq_dist(p, po) - 4*r*r;
	B = sq_dist_pd(p, po)-4*r*r;
	if (A*A-B < 0) {
		return -1;
	}
	h = -A+sqrt(A*A-B);
    
	if (h < -p.z()) {
		return -1;
	}
	return h;
}

void twodim_put_first_fo(vector<vec2d> &p, int fo_size)
{
	double x_mid = Lx*drand48(); //@@@@@@@@@@@ one fractal
	double z_mid = Lz_pd*drand48();
	for (int i=0; i < fo_size ; i++) {
		p[i].shift(x_mid, z_mid);
	}
    for (int i=0; i < fo_size ; i++) {
        p[i].periodic_range_xz(Lx, Lz_pd);
    }
}

void twodim_drop_fo(vector<vec2d> &p, int i_premier, int fo_size)
{
	int i_dernier = i_premier + fo_size;
	double x_rand = Lx*drand48();
	double h_max = -1;
	double y_min = 0;
	for (int i= i_premier ; i < i_dernier; i++) {
		p[i].shift( x_rand, 0);
		p[i].periodic_range(Lx);
		if (y_min > p[i].z()) {
			y_min = p[i].z();
		}
		for (int i_stacked = 0; i_stacked<i_premier; i_stacked++) {
			double h = twodim_collision_height(p[i_stacked], p[i]);
			if (h > h_max) {
				h_max = h;
			}
		}
	}
	
	if (h_max == -1) {
		for (int i=i_premier; i < i_dernier; i++) {
			p[i].shift(0, -y_min);
		}
	} else {
		double y_min = 0;
		for (int i=i_premier; i < i_dernier; i++) {
			p[i].shift(0, h_max);
			if (p[i].z() < y_min){
				y_min = p[i].z();
			}
		}
		if (y_min < 0){
			for (int i=i_premier; i < i_dernier; i++) {
				p[i].shift(0, -y_min);
			}
		}
	}
}

void twodim_fo_random_fall(vector<vec2d> &p, int fo_size)
{
	//    unsigned long n =
    unsigned long fo_num = p.size()/fo_size;
    twodim_put_first_fo(p, fo_size);
	for (int j = 1; j < fo_num; j++) {
		twodim_drop_fo(p, j*fo_size, fo_size);
	}
	
}

void twodim_set_initial_pairs(vector<vec2d> &p)
{
    unsigned long num = p.size();
	double theta;
	for (int i=0; i < num/2; i++) {
		theta = 2*M_PI*drand48();
		p[2*i].set(r*cos(theta), r*sin(theta));
		p[2*i+1].set(-r*cos(theta), -r*sin(theta));
	}
}

double twodim_calc_y_max(vector<vec2d> &p)
{
	double y_max = 0;
	for (int j=0; j < n; j++) {
		if (p[j].z() > y_max) {
			y_max = p[j].z();
		}
	}
	return y_max;
}

void twodim_output_result(ofstream &fout, vector<vec2d> &p,
						  int number_of_floc_for_vf, int fo_size)
{
	for (int k = 0; k < fo_size*number_of_floc_for_vf; k++) {
		//p[k].shift(0, -y_max/2);
		p[k].shift(0, 1);
		//p[i + j*fo_size].shift(-L/2, 0);
		//p[i + j*fo_size].cout();
		//p[i + j*fo_size].out_data();
		cout << "@ 0 " << endl;
		fout << setprecision(15);
		//		fout << p[k].x() << ' ' <<  p[k].z() << ' ' << 0 << endl;
		//            cout << "@ " << col << endl;
		cout << "c " << p[k].x() - Lx/2<< ' ' << 0 <<  ' ' <<  p[k].z() - Lz/2 << endl;
		cout << "@ 3 " << endl;
		
		if ( p[k].x() < 1)
			cout << "c " << p[k].x() + Lx - Lx/2<< ' ' << 0 <<  ' ' <<  p[k].z() - Lz/2 << endl;
		else if ( p[k].x() > Lx - 1)
			cout << "c " << p[k].x() - Lx - Lx/2<< ' ' << 0 <<  ' ' <<  p[k].z() - Lz/2 << endl;
		if ( p[k].z() < 2)
			cout << "c " << p[k].x() - Lx/2<< ' ' << 0 <<  ' ' <<  p[k].z() + Lz_pd - Lz/2 << endl;
		else if ( p[k].z() > Lz_pd - 1)
			cout << "c " << p[k].x() - Lx/2<< ' ' << 0 <<  ' ' <<  p[k].z() - Lz_pd - Lz/2 << endl;
		
	}
    cout << "l " << -Lx/2 << ' ' << 0 << ' ' << -Lz/2 << ' ';
    cout << Lx/2 << ' ' << 0 << ' ' << -Lz/2 << endl;
    cout << "l " << -Lx/2 << ' ' << 0 << ' ' << Lz/2 << ' ';
    cout << Lx/2 << ' ' << 0 << ' ' << Lz/2 << endl;
	cerr << "N = " << fo_size*number_of_floc_for_vf << endl;
}

int twodim_put_fo(vector<vec2d> &p, int i_first, int fo_size)
{
    int max_trial = 10000;
	int i_last = i_first + fo_size;
	vec2d p_tmp;
	int count = 0;
	
try_again1:
	double x_rand = Lx*drand48();
	double y_rand = Lz_pd*drand48();
	for (int i = i_first ; i < i_last; i++) {
		p_tmp = p[i].plus(x_rand, y_rand);
		for( int j = 0; j < i_first ; j++) {
            for (int pb = 0; pb < 4 ; pb++) {
                vec2d p_depl;
                p_depl = p[j];
                switch (pb) {
                    case 0:
                        break;
                    case 1:
						if (p_depl.x() > 0.5*Lx) {
                            p_depl.shift(-Lx,0);
						} else {
                            p_depl.shift( Lx,0);
						}
                        break;
                    case 2:
						if (p_depl.z() > 0.5*Lz_pd) {
                            p_depl.shift(0, -Lz_pd);
						} else {
                            p_depl.shift(0,  Lz_pd);
						}
                        break;
                    case 3:
						if (p_depl.x() > 0.5*Lx) {
                            p_depl.shift(-Lx,0);
						} else {
                            p_depl.shift( Lx,0);
						}
						
						if (p_depl.z() > 0.5*Lz_pd) {
                            p_depl.shift(0, -Lz_pd);
						} else {
                            p_depl.shift(0,  Lz_pd);
						}
                        break;
                }
                if (sq_dist_pd(p_tmp, p_depl) <= 4) {
                    count ++;
                    if (count > max_trial) {
                        return 1;
                    }
                    goto try_again1;
                }
            }
        }
    }
	
	//try_again_move:
	
    double angle = 2*M_PI*drand48();
	vec2d v( cos(angle), sin(angle));
	double lambda = 100;
	for (int i = i_first ; i < i_last; i++) {
        p_tmp = p[i].plus(x_rand, y_rand);
        //p_tmp.periodic_range_xy(L, Lz_pd);
        for ( int j = 0; j < i_first ; j++) {
            double A = v.sq_norm();
            for (int pb = 0; pb < 4; pb++) {
                vec2d p_depl;
                p_depl = p[j];
                switch (pb) {
                    case 0:
                        break;
                    case 1:
						if (p_depl.x() > 0.5*Lx) {
                            p_depl.shift(-Lx,0);
						} else {
                            p_depl.shift( Lx,0);
						}
                        break;
                    case 2:
						if (p_depl.z() > 0.5*Lz_pd) {
                            p_depl.shift(0, -Lz_pd);
						} else {
                            p_depl.shift(0,  Lz_pd);
						}
                        break;
                    case 3:
						if (p_depl.x() > 0.5*Lx) {
                            p_depl.shift(-Lx,0);
						} else {
							p_depl.shift( Lx,0);
						}
						if (p_depl.z() > 0.5*Lz_pd) {
                            p_depl.shift(0, -Lz_pd);
						} else {
                            p_depl.shift(0,  Lz_pd);
						}
                        break;
                }
                double B = dot(v, p_tmp - p_depl);
                vec2d tmp = (p_tmp - p_depl);
                double C = tmp.sq_norm() - 4;
                double DD = B*B - A*C;
                if (DD > 0) {
                    double sqrt_DD = sqrt(DD);
                    double lambda1 = (-B + sqrt_DD)/A;
                    double lambda2 = (-B - sqrt_DD)/A;
					if (lambda1 > 0 && lambda1 < lambda) {
                        lambda = lambda1;
					}
					if (lambda2 > 0 && lambda2 < lambda) {
                        lambda = lambda2;
					}
                }
            }
        }
	}

	if ( lambda == 100){
		goto try_again1;
    }
	
	for (int i = i_first; i < i_last; i++){
		p[i].shift(x_rand, y_rand);
		//p[i] += lambda*v;
        p[i].periodic_range_xz(Lx, Lz_pd);
    }
	return 0;
}

int twodim_fo_random_put(vector<vec2d> &p,
						 int floc_size,
						 int number_of_particle_for_vf)
{
	//	unsigned long fo_num = p.size()/fo_size;
	twodim_put_first_fo(p, floc_size);
	//int j_max = (int)(total_number/fo_size);
	for (int j = 1; j < number_of_particle_for_vf; j++) {
		twodim_put_fo(p, j*floc_size, floc_size);
	}
	return (int)n;
}

#endif
