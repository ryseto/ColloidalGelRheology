//
//  threedim.h
//  ColloidalGelRheology
//
//  Created by Ryohei Seto on 2012/08/17.
//
//

#ifndef ColloidalGelRheology_threedim_h
#define ColloidalGelRheology_threedim_h
#include "tools.h"
#include "vec3d.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
using namespace std;
extern double r;
extern int n;
extern double L;


void sediment_make3dNetwork(int rank, int m, int rsd);
//int ipow(int p, int q);
void make_fractal_object(int m, vector<vec3d> &p);
void fo_random_fall( vector<vec3d> &p, int fo_size);
void make_fractal_object(int m, vector<vec3d> &p);
void set_initial_pairs(vector<vec3d> &p);
//double calc_z_max(vector<vec3d> &p);
void output_result(ofstream &fout, vector<vec3d> &p, int fo_size);
void prepare_fout(ofstream &fout, int rank, int rsd);

double collision_value(vec3d const &po, vec3d const &p, vec3d v);
double collision_height(vec3d po, vec3d p);
void put_first_fo(vector<vec3d> &p, int fo_size);
void drop_fo(vector<vec3d> &p, int i_premier, int fo_size);


//void make3dNetwork(int argc, char ** argv)
void sediment_make3dNetwork(int rank, int m, int rsd){
	if (rank < 1){
		cerr << "rank < 1 is not allowed." << endl;
		exit(1);
	}
    int k = m - rank; // 粒子数 2^k
    //	int rsd = atoi(argv[7]);
	srand48(rsd);
	int fo_size = ipow(2, rank);
	n = ipow(2, rank + k);
	vector<vec3d> p;
	p.resize(n);
	
	cerr << "fractal size = " << fo_size << endl;
	cerr << "number of fractal object = " << n/fo_size << endl;
	cerr << "number of particle = " << n << endl;
	
	// 出力ファイルの準備
	ofstream fout;
    prepare_fout(fout, rank, rsd);
	
	// 粒子をランダムな角で対にする。
	// n/2個の対を生成
	set_initial_pairs(p);
	
	// 与えられたランクのクラスターを作る。
	if ( rank >= 2){
		for( int i = 2; i <= rank; i++){
			int fo_size_tmp = ipow(2, i);
			make_fractal_object(fo_size_tmp, p);
		}
	} else if ( rank == 1 ){
		int fo_size_tmp = 2;
		make_fractal_object(fo_size_tmp, p);
		
	}
	
	// ランダムに落として積み上げる。
	fo_random_fall(p, fo_size);
	
	// 結果を fout に出力。
	
	output_result(fout, p, fo_size);
	fout.close();
	return ;
}

inline vec3d randUniformSphere(double r){
	double z = 2.*drand48() - 1.;
	double phi = 2.*M_PI*drand48();
	double sin_theta = sqrt(1.-z*z);
	return vec3d( r*sin_theta*cos(phi),  r*sin_theta*sin(phi), r*z);
}

void prepare_fout(ofstream &fout, int rank, int rsd){
	char name_string[64];
	sprintf(name_string, "L%d_rank%d_rand%d", (int)(L), rank, rsd) ;
	string fn_fout = "sfn3D_" + (string)name_string + ".dat";
	fout.open(fn_fout.c_str());
}

/////////////////////////////////////////////////////////////////
double collision_value(vec3d const &po, vec3d const &p, vec3d v){
	double l1, l2;
	double A, B;
	vec3d dp = po - p;
	// A = (p.x()-po.x())*cos(theta)+(p.y()-po.y())*sin(theta);
	A = - dot(dp,v);
	B =   dp.sq_norm() - 4*r*r;
	//B = sq_dist_pd(p, po) - 4*r*r;
	//B = sq_dist(p, po) - 4*r*r;
	if ( A*A - B < 0){
		return 0;
	}
	l1 = - A - sqrt( A*A-B );
	l2 = - A + sqrt( A*A-B );
	
	if (l1 > l2) //@@@@@@@@@@ l1 > l2
		return l1;
	else
		return l2;
}

void make_fractal_object(int m, vector<vec3d> &p){
	// fo_size_tmp, p
	// m はフラクタルのサイズ
	for (int k=0; k < n; k+= m){
		vec3d v = randUniformSphere(1.0);
		double l;
		double l_max = 0;
		for (int i = k; i < k+m/2; i++){ //固定する粒子
			for (int j = k+m/2; j < k+m; j++){//平行移動する粒子
				l = collision_value(p[i], p[j], v);
				if( l >= l_max ){
					l_max = l;
				}
			}
		}
		// 最大の移動パラメータ l の最大値で平行移動
		for (int j = k+m/2; j < k+m; j++){
			p[j] = p[j] + l_max*v; //.plus( l_max*, l_max*sin(theta), 0 );
		}
		vec3d cm(0,0,0);//重心
		for(int i=k; i<k+m; i++){
			cm += p[i]/m;
		}
		for(int i=k; i<k+m; i++){
			p[i] = p[i] - cm;
		}
	}
}

//int ipow(int p, int q){
//	int x = 1;
//	for (int i=0; i < q; i++){
//		x *= p;
//	}
//	return x;
//}

double collision_height(vec3d po, vec3d p){
	double h, A, B;
	vec3d p_po = p - po;
	if (abs(p_po.x) > L/2){
		if (p_po.x > 0)
			p_po.x += -L;
		else
			p_po.x += +L;
	}
	if (abs(p_po.y) > L/2){
		if (p_po.y > 0)
			p_po.y += -L;
		else
			p_po.y += +L;
	}
	
	A = p.z - po.z;
	B = p_po.sq_norm() - 4*r*r;
	//B = sq_dist_pd(p, po) - 4*r*r;
	if (A*A-B < 0){
		return -1;
	}
	h = - A + sqrt( A*A-B );
	if (h < -p.z){
		return -1;
	}
	return h;
}

void put_first_fo(vector<vec3d> &p, int fo_size){
	//double x_rand = L*drand48();
	double x_rand = 0; //@@@@@@@@@@@ one fractal
	double y_rand = 0;
	double z_min = 0;
	for (int i=0; i < fo_size ; i++){
		//p[i].change(x_rand, y_rand , 0);
		p[i].x += x_rand, p[i].y += y_rand;
		if (p[i].z < z_min){
			z_min = p[i].z;
		}
	}
    
	for (int i=0; i < fo_size ; i++){
		//p[i].periodic_range(L);
		//p[i].change(0,0, -z_min);
		p[i].z += -z_min;
	}
    
}

void drop_fo(vector<vec3d> &p, int i_premier, int fo_size){
	int i_dernier = i_premier + fo_size;
	double x_rand = L*drand48();
	double y_rand = L*drand48();
	double h_max = -1;
	double z_min = 0;
	for (int i= i_premier ; i < i_dernier; i++){
		//p[i].change( x_rand, y_rand, 0);
		p[i].x += x_rand, p[i].y += y_rand;
		//p[i].periodic_range(L);
		if (z_min > p[i].z){
			z_min = p[i].z;
		}
		for(int i_stacked = 0; i_stacked<i_premier; i_stacked++){
			double h = collision_height(p[i_stacked], p[i]);
			if(h > h_max)
				h_max = h;
		}
	}
	
	if (h_max == -1){
		for (int i=i_premier; i < i_dernier ; i++){
			//p[i].change(0, 0, -z_min);
			p[i].z += -z_min;
		}
	} else {
		double z_min = 0;
		for (int i=i_premier; i < i_dernier ; i++){
			//p[i].change(0, 0, h_max);
			p[i].z += h_max;
			if (p[i].z < z_min){
				z_min = p[i].z;
			}
			
		}
		if (z_min < 0 ){
			for (int i=i_premier; i < i_dernier ; i++){
				//p[i].change(0, 0, -z_min);
				p[i].z += -z_min;
			}
		}
	}
	for (int i=i_premier; i < i_dernier ; i++){
		if (p[i].x > L)
			p[i].x -= L;
		else if (p[i].x < 0)
			p[i].x += L;
		if (p[i].y > L)
			p[i].y -= L;
		else if (p[i].y < 0)
			p[i].y += L;
	}
}

void fo_random_fall( vector<vec3d> &p, int fo_size){
	int n = (int)p.size();
	int fo_num = n/fo_size;
	put_first_fo(p, fo_size);
	
	for (int j = 1; j < fo_num; j++){
		drop_fo(p, j*fo_size, fo_size);
	}
	
}

void set_initial_pairs(vector<vec3d> &p){
	int num = (int)p.size();
	
    //	double theta, phi;
	for (int i=0; i < num/2; i++){
		//theta = M_PI*drand48();
		//phi = 2*M_PI*drand48();
		//vec3d dp(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
		vec3d dp = randUniformSphere(r);
		p[2*i] = dp;
		p[2*i+1] = -dp;
		
		//p[2*i + 1] = -p[2*i];
	}
}

//double calc_z_max(vector<vec3d> &p){
//	double z_max = 0;
//	for (int j=0; j < n; j++){
//		if( p[j].z() > z_max ){
//			z_max = p[j].z();
//		}
//	}
//	return z_max;
//}

void output_result(ofstream &fout, vector<vec3d> &p, int fo_size){
	//double y_max = calc_y_max(p);
	double lowestZ = 1.0;
	//fout << "r 0.5" << endl;
	for (int j=0; j < n/fo_size; j++){
		for (int i=0; i < fo_size ; i++){
			int k = i + j*fo_size;
			//p[k].shift(0, -y_max/2);
			//p[k].change(0, 0, lowestZ);
			p[k].z += lowestZ;
			//p[i + j*fo_size].shift(-L/2, 0);
			//			p[i + j*fo_size].cout();
			//p[i + j*fo_size].out_data();
			fout << setprecision(15);
			//fout << "c " << p[k].x() << ' ' <<  p[k].y()  << ' ' <<  p[k].z() - 50 << endl;
			while (p[k].x < 0){
				p[k].x += L;
			}
			
			while (p[k].x > L){
				p[k].x -= L;
			}
			
			while (p[k].y < 0){
				p[k].y += L;
			}
			
			while (p[k].x > L){
				p[k].y -= L;
			}
			fout << p[k].x << ' ' <<  p[k].y  << ' ' <<  p[k].z << endl;
		}
	}
}

#endif
