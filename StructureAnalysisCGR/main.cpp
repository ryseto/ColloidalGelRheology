//
//  main.cpp
//  StructureAnalysisCGR
//
//  Created by Ryohei SETO on 6/21/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <sstream>

#include <math.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
//#include "CoreFoundation/CoreFoundation.h"
#include <sys/stat.h> // for stat and mkdir
#include "vec3d.h"
#include "grid2.h"
using namespace std;

double lx, lz;
double lx_half, lz_half;

void outputZaverageVolumeFraction(vector <vec3d> &p, ofstream &fout_phi_z);


bool import_compaction(vector <vec3d> &p, ifstream &fin,
					   bool back_bone, double f_threshold){
	unsigned long num_particle;
	char buf[10];
	char buf2[10];
	double info[13];
	char buf_line[512];
	fin >> buf >> buf2;
	while (buf[0]  != '#'){
		cerr  << "** " << buf_line << endl;
		cerr << "stop: " << buf << endl;
		exit(1);
	}
	bool equilibrium = false;
	if (buf2[0] == 'e'){
		equilibrium  = true;
		cerr << "eq" << endl;
	} else {
		cerr << "not eq" << endl;
	}
	int num_infomation = 10;
	for(int i=0;i < num_infomation; i++){
		fin >> info[i];
	}

	double vf = info[0];
	cerr << "vf = " << vf << endl;
	//	stress_z = force_z / (sy->lx*2.0);
	double stress = info[1];
	double wall_force = stress*lx*2.0;
	fin.getline(buf_line, 1000);

	double z_bot = info[5];
	double z_top = info[6];
	
	lz = z_top-z_bot;
	
	fin >> buf >> num_particle;
	cerr << "num_particle = " << num_particle << endl;
	double x, y, z;
	vector< vec3d > p_import;
	p_import.resize(num_particle);
	//int tmp;
	for (int i = 0; i < num_particle ; i++){
		fin >> x >> y >> z ;
		//		cerr << x << ' ' <<  y << ' ' << z << endl;
		p_import[i].set(x+lx/2, y, z - z_bot);
	
		if (  z - z_bot < 0){
			cerr << "t" << z_top << endl;
			cerr << "b" << z_bot << endl;
			cerr << z << endl;
			cerr << z + z_bot << endl;
			exit(1);
		}
		fin.getline(buf_line, 512);
	}
	int num_bond;
	fin >> buf >> num_bond;
	cerr << "num_bond " << num_bond << endl;

//	fout_conf << i_p0 << ' ' << i_p1 << ' '; //1,2
//	fout_conf << initial_bond << ' ' << (*bond_iter)->status << ' '; //3,4
//	fout_conf << (*bond_iter)->val_F_norm() << ' ';//5
//	fout_conf << (*bond_iter)->val_F_slid() << ' ';//6
//	fout_conf << (*bond_iter)->val_M_bend() << ' ';//7
//#ifndef TWODIMENSION
//	fout_conf << (*bond_iter)->val_M_tors() << ' ';//8
//#else
//	fout_conf << 0 << ' ';                         //8
//#endif
//	fout_conf << (*bond_iter)->cnt_regeneration     // 9

	int p1, p2;
	int init_bond;
	int bond_status;
	double f_norm;
	double f_slid;
	double M_bend;
	double M_tors;
	int cnt_regeneration;
	vector<int> bb_check;
	bb_check.resize(num_particle);
	vector <int> bond_p1;
	vector <int> bond_p2;
	vector <double> bond_stress;
	bond_p1.resize(num_bond);
	bond_p2.resize(num_bond);
	bond_stress.resize(num_bond);
	//	double sum_stress = 0;
	//	int count = 0;
	vector <double> particle_stress;
	particle_stress.resize(num_particle);
	for (int i = 0; i < num_bond ; i++){
		fin >> p1 >> p2 ;
		fin >> init_bond >> bond_status;
		fin >> f_norm >> f_slid >> M_bend >> M_tors;
		fin >> cnt_regeneration;
		bond_p1[i] = p1;
		bond_p2[i] = p2;
//		fin.getline(buf_line, 20);
//		cerr << buf_line << endl;
		particle_stress[p1] += 0.5*(abs(f_norm)+abs(f_slid)+abs(M_bend));
		particle_stress[p2] += 0.5*(abs(f_norm)+abs(f_slid)+abs(M_bend));
		//		bond_stress[i] = f_norm;
		//		if (abs(f_norm) > 1e-6 ){
		//			sum_stress += abs(f_norm);
		//			count ++;
		//		}
	}
	fin.getline(buf_line, 512);
	cerr << buf_line << endl;
	//	double average_stress = sum_stress / count;
	//	if (count > 0){
	//			cout << vf << ' '<< average_stress << endl;
	//	}
	//	double lo = 4.0*pow(vf,-0.7);
	if (back_bone){
		for (int i = 0; i < num_particle ; i++){
			if (particle_stress[i] > f_threshold*wall_force ){
				bb_check[i] = 1;
			}
		}
	}

	if (equilibrium){
		p.clear();
		if (back_bone){
			for (int i=0; i< num_particle; i++){
				if (bb_check[i] == 1){
					p.push_back(p_import[i]);
					cout << "@ 3\n";
				} else {
					cout << "@ 0\n";
				}
				cout << "c " << p_import[i].x - lx/2 << " 0 " << p_import[i].z - lz/2 << endl;
			}
			cout << endl;
			//		cerr << "BB"  << num_particle << endl;
		} else {
			for (int i=0; i< num_particle; i++){
				p.push_back(p_import[i]);
			}
		}
//		for (int i = 0; i < num_bond ; i++){
//			cout << bond_stress[i]  << ' ';
//		}
//		cout << endl;
//		cerr << "N : Nbb = " << num_particle << ' ' << p.size() << endl;
//		cerr << "average_stress =" << average_stress << endl;
		
	}
	
	num_particle = p.size();
//	cerr << "num_particle = " << num_particle << endl;
//	if ( num_particle > 0 ){
//		cerr << num_particle << endl;
//		exit(1);
//	}

	return equilibrium;
}

bool import_shear(vector <vec3d> &p, ifstream &fin, double &strain){
	int num_particle;
	char buf[10];
	double info[13];
	char buf_line[512];
	fin >> buf ;
	if (buf[0] != '#'){
		cerr << "stop: " << buf << endl;
		exit(1);
	}
	fin >> buf;
	bool equilibrium = false;
	if (buf[0] == 'e')
		equilibrium  = true;
	int num_infomation = 13;
	for(int i=0;i<num_infomation; i++){
		fin >> info[i];
		cerr << info[i] << ' ';
	}
	cerr << endl;
	strain = info[5];
	double z_bot = info[6];
	double z_top = info[7];
	lz = z_top-z_bot;
	cerr << "lz="<< lz << endl;
	
	fin >> buf >> num_particle;
	
	cerr << "N = " << num_particle << endl;
	double x, y, z;
	p.resize(num_particle);
	//int tmp;
	for (int i = 0; i < num_particle ; i++){
		fin >> x >> y >> z ;
		p[i].set(x+lx/2, y, z - z_bot);
		if (  z - z_bot < 0){
			cerr << "t" << z_top << endl;
			cerr << "b" << z_bot << endl;
			cerr << z << endl;
			cerr << z + z_bot << endl;
			exit(1);
		}
		fin.getline(buf_line, 512);
		
	}
	int num_bond;
	fin >> buf >> num_bond;
	for (int i = 0; i < num_bond ; i++){
		fin.getline(buf_line, 512);
		//cerr << buf_line << endl;
	}
	fin.getline(buf_line, 512);
	return equilibrium;
}

double sq_distancePeriodicBoundaries(const vec3d &p1, vec3d p2, bool z_pd ){
	if (abs(p1.x-p2.x) > lx_half) {
		if (p2.x > p1.x) {
			p2.x -= lx;
		} else {
			p2.x += lx;
		}
	}
	if (z_pd) {
		if (abs(p1.z - p2.z) > lz_half){
			if (p2.z > p1.z) {
				p2.z -= lz;
			} else {
				p2.z += lz;
			}
		}
	}
	return (p1-p2).sq_norm_xz();
}

void densityDensityCorrelationFunction(double r_min,
									   int resolution,
									   double **ddc,
									   vector <vec3d> &p,
									   bool z_pd){

	Grid2 grid;
	double r_max;
	if (lx < lz){
		r_max = lx/2;
	} else {
		r_max = lz/2;
	}
	double cell_h = 2;
	unsigned long num_particle = p.size();
	double ly = 0;
	grid.init(num_particle, lx, ly, lz, cell_h);
	for (int i=0; i < num_particle; i++){
		grid.entry(p[i], i);
	}
	srand(14235);

	double d_logr = (log(r_max)-log(r_min))/resolution;
	int  mmax = 500000;
	vec3d p1, p0;
	vector<int> neighbor;
	vector<int> neighbor1;
	for (int k=0; k < resolution; k++){
		double r = r_min*exp(k*d_logr);
		int count_p0_inside = 0;
		int count_p1_inside = 0;
		while (count_p0_inside < mmax) {
			p0.set(lx*drand48(), 0, 1+(lz-2)*drand48()); // 1 < z < lz-1
			neighbor.clear();
			grid.get_neighbor_list(p0, neighbor);
			for (int i=0; i<neighbor.size(); i++){
				double sq_distance = sq_distancePeriodicBoundaries(p0, p[neighbor[i]], z_pd);
				if (sq_distance <= 1) {
					bool p1_success = true;
					double theta = 2*M_PI*drand48();
					p1.x = p0.x+r*cos(theta);
					p1.z = p0.z+r*sin(theta);
					if (z_pd) {
						if (p1.z < 0) {
							p1.z += lz;
						} else if (p1.z >= lz) {
							p1.z -= lz;
						}
					} else {
						if (p1.z > lz-1 || p1.z < 1) {
							p1_success = false;
						}
					}
					if (p1.x < 0) {
						p1.x += lx;
					} else if (p1.x >= lx) {
						p1.x -= lx;
					}
					if (p1_success) {
						count_p0_inside ++;
						neighbor1.clear();
						grid.get_neighbor_list(p1, neighbor1);
						for (int j=0; j<neighbor1.size(); j++) {
							double sq_distance1 = sq_distancePeriodicBoundaries(p1, p[neighbor1[j]], z_pd);
							if (sq_distance1 < 1) {
								count_p1_inside ++;
								break;
							}
						}
					}
					break;
				}
			}
		}
		//fout << r << ' ' << (double)count_p1_inside / count_p0_inside << ' ' << (double)count_p0_inside/mmax<< endl;
		ddc[k][0] = r;
		ddc[k][1] = (double)count_p1_inside/count_p0_inside;
	}
}


double circleArea(double r, double dz){
	double circle_area = M_PI*r*r;
	if (dz < r) {
		double theta = acos(dz/r);
		double lack_area = circle_area*theta/M_PI-dz*r*sin(theta);
		return circle_area-lack_area;
	} else {
		return circle_area;
	}
}

void pairCorrelationFunction(double r_min,
							 int resolution,
							 double **gr,
							 vector <vec3d> &p,
							 bool z_pd){
	double r_max;
	if (lx < lz) {
		r_max = lx/2-5;
	} else {
		r_max = lz/2-5;
	}
	cerr << " r_max = " << r_max << endl;
	
	double d_logr = (log(r_max)-log(r_min))/resolution;
	double logr_min = log(r_min);
	for(int k=0; k<resolution; k++){
		gr[k][1] = 0;
	}
	for (int i=0; i<p.size(); i++){
		for (int j=0; j<p.size(); j++){
			if (i != j) {
				vec3d dp = p[i]-p[j];
				if (z_pd) {
					if (dp.z > 0) {
						double periodic_dz = dp.z-lz;
						if (-periodic_dz < dp.z){
							dp.z -= lz;
						}
					} else {
						double periodic_dz = dp.z+lz;
						if (periodic_dz < -dp.z){
							dp.z += lz;
						}
					}
				}
				if (dp.x > 0) {
					double periodic_dx = dp.x-lx;
					if (-periodic_dx < dp.x){
						dp.x -= lx;
					}
				} else {
					double periodic_dx = dp.x+lx;
					if (periodic_dx < -dp.x){
						dp.x += lx;
					}
				}
				
				double distance = dp.norm();
				if (distance < r_max) {
					int k = (log(distance)-logr_min)/d_logr;
					double r_k = exp(logr_min+d_logr*k);
					double r_k1 = exp(logr_min+d_logr*(k+1));
					double d_area;
					if (z_pd == false) {
						double dz_bot = p[i].z-1;
						double dz_top = (lz-1)-p[i].z;
						if (dz_bot>dz_top){
							d_area = circleArea(r_k1, dz_top)-circleArea(r_k, dz_top);
						} else {
							d_area = circleArea(r_k1, dz_bot)-circleArea(r_k, dz_bot);
						}
						if (dz_bot < r_k1 && dz_top < r_k1){
							cerr << dz_bot << " and " << dz_top << " < " << r_k1 << endl;
							exit(1);
						}
					} else {
						d_area = M_PI*(r_k1*r_k1 - r_k*r_k);
					}
					gr[k][1] += (M_PI/d_area);
				}
			}
		}
	}
//	for(int k=0; k<resolution; k++){
//
//		double r_1 = exp(logr_min+(k+1)*d_logr);
//		double area = M_PI*(r_1*r_1-r*r);

//		gr[k][1] = gr[k][1]/area;
//	}
	unsigned long np = p.size();
	for(int k=0; k<resolution; k++){
		double r = exp(logr_min+k*d_logr);
		gr[k][0] = r;
		gr[k][1] = gr[k][1]/np;
	}
	//	double r = r_min*exp(k*d_logr);
	//fout << r << ' ' << (double)count_p1_inside / count_p0_inside << ' ' << (double)count_p0_inside/mmax<< endl;
	//ddc[k][0] = r;
	//ddc[k][1] = (double)count_p1_inside/count_p0_inside;
	//}
}




void importInitConfig(char *simu_name, string path, vector <vec3d> &p){
	ifstream fin;
	fin.open(path.c_str());
	double x, z;
	int tmp;
	while (fin >> x >> z >> tmp) {
		cerr << x << ' ' << z << endl;
		p.push_back(vec3d(x, 0, z));
	}
	cerr << "N = " << p.size() << endl;
}

void evaluateVolumeFractionZ(char *simu_name, string path, vector <vec3d> &p){
	string simu_name_str = simu_name;
	ofstream fout_phi_z;
	string phiz_filename = "phi_z_"+simu_name_str+".dat";
	fout_phi_z.open(phiz_filename.c_str());
	outputZaverageVolumeFraction(p, fout_phi_z);
	fout_phi_z.close();
}

void evaluateVolumeFractionZ_compt(char *simu_name, string path, vector <vec3d> &p){

}

void evaluateDDC_InitConfig(char *simu_name, string path, vector <vec3d> &p,
							int resolution, double r_min,
							bool z_periodic_boundary){
	ofstream fout;
	string simu_name_str = simu_name;
	string outpuffile = "dcf_"+simu_name_str+".dat";
	cerr << outpuffile << endl;
	fout.open(outpuffile.c_str());
	double phi = M_PI*p.size()/(lx*lz);//
	cerr << "phi = " << phi << endl;
	double **ddc;
	ddc = new double* [resolution];
	for (int k=0; k<resolution; k++) {
		ddc[k] = new double [2];
	}
	densityDensityCorrelationFunction(r_min,
									  resolution, ddc,
									  p, z_periodic_boundary);
	if (fout.is_open()) {
		cerr << " opened successfully"  << endl;
		for (int k=0; k< resolution; k++){
			fout << ddc[k][0] << ' ';
			fout << ddc[k][1]/phi << endl;
			cerr << ddc[k][1] << endl;
		}
	} else {
		cerr << " fail to open" <<  outpuffile << endl;
	}
	fout.close();
	for (int k=0; k < resolution; k++){
		delete [] ddc[k];
	}
	delete [] ddc;
	//densityDensityCorrelationFunction(p, lx, lz);
	//densityDensityCorrelationFunction2(p, lx,lz, fout);
	fout.close();
}

void evaluatePairCorrelation_InitConfig(char *simu_name, string path, vector <vec3d> &p,
										int resolution, double r_min,
										bool z_periodic_boundary){

	ofstream fout;
	string simu_name_str = simu_name;
	string outpuffile = "paircol_"+simu_name_str+".dat";
	cerr << outpuffile << endl;
	fout.open(outpuffile.c_str());
	double phi = M_PI*p.size()/(lx*lz);//
	cerr << "phi = " << phi << endl;
	double **gr;
	gr = new double* [resolution];
	for (int k=0; k<resolution; k++) {
		gr[k] = new double [2];
	}
	pairCorrelationFunction(r_min,
							resolution, gr,
							p, true);
	if (fout.is_open()) {
		cerr << " opened successfully"  << endl;
		for (int k=0; k<resolution; k++){
			if (gr[k][1] != 0){
				fout << gr[k][0] << ' ';
				fout << gr[k][1]/phi << endl;
			}
		}
	} else {
		cerr << " fail to open" <<  outpuffile << endl;
	}
	fout.close();
	for (int k=0; k < resolution; k++){
		delete [] gr[k];
	}
	delete [] gr;
	//densityDensityCorrelationFunction(p, lx, lz);
	//densityDensityCorrelationFunction2(p, lx,lz, fout);
	fout.close();
}


void evaluateDDC_Compaction(char *simu_name, string path, vector <vec3d> &p,
							int resolution, double r_min,
							bool z_periodic_boundary,
							bool back_bone,
							int i){
	double **ddc;
	ddc = new double* [resolution];
	for (int k=0; k < resolution; k++){
		ddc[k] = new double [2];
	}
	if (back_bone){
		cerr << "bb" << i << endl;
	}
	if (back_bone == false || i >= 1){
		densityDensityCorrelationFunction(r_min, resolution, ddc,
										  p, z_periodic_boundary);
	}
	double phi = (M_PI*p.size())/(lx*lz);
	ofstream fout;
	ofstream fout_phi_z;
	ostringstream ddc_filename;
	ostringstream phiz_filename;
	if (back_bone == false || i >= 1){
		if (back_bone){
			cerr << "bb" << i << endl;
			ddc_filename << simu_name <<  "/bb_ddc" << i << ".dat";
		} else {
			ddc_filename << simu_name <<  "/ddc" << i << ".dat";
			phiz_filename << simu_name <<  "/phi_z_" << i << ".dat";
		}
		fout.open(ddc_filename.str().c_str());
		fout << "# " << phi << endl;
		for (int k=0; k< resolution; k++){
			fout << ddc[k][0] << ' ';
			fout << ddc[k][1]/phi << endl;
		}
		fout.close();
	}
	for (int k=0; k < resolution; k++){
		delete [] ddc[k];
	}
	delete [] ddc;
	////////////////////////////////////////////////////////////////////////
	fout_phi_z.open( phiz_filename.str().c_str());
	outputZaverageVolumeFraction(p, fout_phi_z);
	fout_phi_z.close();
}


void evaluatePairCol_Compaction(char *simu_name, string path, vector <vec3d> &p,
								int resolution, double r_min,
								bool z_periodic_boundary,
								bool back_bone,
								int i){
	double **gr;
	gr = new double* [resolution];
	for (int k=0; k < resolution; k++){
		gr[k] = new double [2];
	}
	pairCorrelationFunction(r_min, resolution, gr,
							p, false);
	double phi = (M_PI*p.size())/(lx*lz);
	ofstream fout;
	ostringstream gr_filename;


	gr_filename << simu_name <<  "/gr" << i << ".dat";
	fout.open(gr_filename.str().c_str());
	fout << "# " << phi << endl;
	for (int k=0; k< resolution; k++){
		fout << gr[k][0] << ' ';
		fout << gr[k][1]/phi << endl;
	}
	fout.close();
	
	for (int k=0; k < resolution; k++){
		delete [] gr[k];
	}
	delete [] gr;
	////////////////////////////////////////////////////////////////////////
	//	fout_phi_z.open( phiz_filename.str().c_str());
	//	outputZaverageVolumeFraction(p, fout_phi_z);
	//	fout_phi_z.close();
}

void evaluateBoxCounting(char *simu_name, string path, vector <vec3d> &p,
						 int resolution, double r_min,
						 bool z_periodic_boundary,
						 bool back_bone,
						 int i){
	r_min = 2;
	double r_max;
	if (lx < lz) {
		r_max = lx/2-5;
	} else {
		r_max = lz/2-5;
	}
	cerr << " r_max = " << r_max << endl;
	
	double d_logr = (log(r_max)-log(r_min))/resolution;
	double logr_min = log(r_min);
	ofstream fout;
	ostringstream bc_filename;
	bc_filename << simu_name <<  "/bc" << i << ".dat";
	fout.open(bc_filename.str().c_str());
	
	for (int k=0; k<resolution; k++){
		double l = r_min*exp(k*d_logr);
		double slide_step = l/3;
		unsigned long np = p.size();
		vector<double> density;
		for (double zz=1; zz<lz-1-l; zz+=slide_step) {
			for (double xx=0; xx<lx;  xx+=slide_step) {
				int np_in_box = 0;
				for (int j=0; j< np; j++){
					double xmax = xx+l;
					if (xmax > lx) {
						if (p[j].z >= zz &&
							p[j].z < zz+l) {
							if (p[j].x >= xx
								|| p[j].x < xmax-lx){
								np_in_box ++;
							}
						}
					} else {
						if (p[j].x >= xx &&
							p[j].x < xmax &&
							p[j].z >= zz &&
							p[j].z < zz+l) {
							np_in_box ++;
						}
					}
				}
				density.push_back(np_in_box*M_PI/(l*l));
			}
		}
		cerr << "box size" << density.size() << endl;
		double phi = (M_PI*p.size())/(lx*lz);
		double sum = 0;
		for (int k=0; k<density.size(); k++){
			sum += density[k];
		}
		double mean = sum/density.size();
		double sum_sq_deviation = 0;
		for (int k=0; k<density.size(); k++){
			double deviation = density[k] - mean;
			sum_sq_deviation += deviation*deviation;
		}
		double std_dev = sqrt(sum_sq_deviation/(density.size()-1));
		
		fout << phi << ' ' << l << ' ' << mean << ' ' << std_dev/mean << ' ' << density.size()  << endl;

//		for (int k=0; k<density.size(); k++){
//			fout << density[k] << endl;
//		}
	

	}
	fout.close();
	////////////////////////////////////////////////////////////////////////
	//	fout_phi_z.open( phiz_filename.str().c_str());
	//	outputZaverageVolumeFraction(p, fout_phi_z);
	//	fout_phi_z.close();
}

					
					

void outputZaverageVolumeFraction(vector <vec3d> &p, ofstream &fout_phi_z){
	int nz = 20;
	int *z_cnt;
	z_cnt = new int [nz];
	for (int iz=0; iz < nz ; iz++){
		z_cnt[iz] = 0;
	}
	double dz = lz / nz;
	for (int i=0; i < p.size(); i++){
		if ( p[i].z > 0 &&
			p[i].z < lz  ){
			int iz = (int)( p[i].z/dz);
			z_cnt[iz]++;
		}
	}
	for (int iz = 0; iz < nz ; iz++){
		fout_phi_z <<  iz*dz + 0.5*dz - lz/2 << ' ' << z_cnt[iz]*M_PI / (lx * dz) << endl;
	}
	fout_phi_z.close();
	delete [] z_cnt;
}

int main(int argc, const char * argv[])
{
	char datatype = argv[1][0];
	string path = argv[2];
	double f_threshold = 0;
	char analyze_type = argv[3][0];
	double phi;
	double r_min = 1;
	double r_max;
	bool back_bone = false;
	int resolution = 100;
	bool z_periodic_boundary;
	/* import_type
	 *  0: x z clusterid
	 *  1:
	 */
	char simu_name[256];
	if (datatype =='i') {
		z_periodic_boundary = true;
		unsigned long i0 = path.find_first_of("2D_W*H*");
		if (i0 != 1) {
			unsigned long i1 = path.find_first_of("W",i0+1)+1 ;
			unsigned long i2 = path.find_first_of("H",i1);
			unsigned long i3 = path.find_first_of("_",i2);
			char width_str[4];
			char height_str[4];
			sprintf(width_str, "%s", (path.substr(i1,i2-i1)).c_str());
			sprintf(height_str, "%s", (path.substr(i2+1,i3-i2-1)).c_str());
			lx = atof(width_str);
			lz = atof(height_str);
		}
		unsigned long i5 = path.find_last_of(".dat");
		sprintf(simu_name, "%s", path.substr(0,i5-3).c_str());
	} else if (datatype =='c') {
		z_periodic_boundary = false;
		unsigned long i0 = path.find_first_of("D");
		if (i0 != 1){
			unsigned long i1 = path.find_first_of("W",i0)+1 ;
			unsigned long i2 = path.find_first_of("H",i1);
			unsigned long i3 = path.find_first_of("_",i2);
			char width_str[4];
			char height_str[4];
			sprintf(width_str, "%s", (path.substr(i1, i2-i1)).c_str());
			sprintf(height_str, "%s", (path.substr(i2+1,i3-i2-1)).c_str());
			lx = atof(width_str);
			lz = atof(height_str);
		}
		unsigned long i4 = path.find_first_of("conf_");
		unsigned long i5 = path.find_last_of(".dat");
		sprintf(simu_name, "%s", (path.substr(i4+5,(i5-3)-(i4+5))).c_str());
	} else if (datatype =='s') {
		z_periodic_boundary = false;
		//		import_type = 0;
		unsigned long i1 = path.find_first_of("W") + 1;
		unsigned long i2 = path.find_first_of("H",i1);
		unsigned long i3 = path.find_first_of("_",i2);
		char width_str[4];
		char height_str[4];
		sprintf(width_str, "%s", (path.substr(i1,i2-i1)).c_str());
		sprintf(height_str, "%s", (path.substr(i2+1,i3-i2-1)).c_str());
		lx = atof(width_str);
		lz = atof(height_str);
		
		unsigned long i4 = path.find_first_of("conf_");
		unsigned long i5 = path.find_last_of(".dat");
		cerr << i4 << ' ' << i5 << endl;
		sprintf(simu_name, "%s", (path.substr(i4+5,(i5-3)-(i4+5))).c_str());
		cerr << simu_name << endl;
		//conf_simu1_floc16_sed_L200_R1_H300_1.dat
		//cerr << lx << ' ' << lz << endl;
	}
		
	lx_half = lx/2;
	lz_half = lz/2;
	
	ifstream fin;
	fin.open( path.c_str());
	vector <vec3d> p;
	if (datatype == 'i'){
		importInitConfig(simu_name, path, p);
		if (analyze_type == 'g'){
			evaluatePairCorrelation_InitConfig(simu_name, path, p,
											   resolution, r_min,
											   z_periodic_boundary);
		} else if (analyze_type == 'd') {
			evaluateDDC_InitConfig(simu_name, path, p,
								   resolution, r_min,
								   z_periodic_boundary);
			evaluateVolumeFractionZ(simu_name, path, p);
		}
	} else if (datatype == 'c'){
		struct stat st;
		if (stat(simu_name, &st) == 0) {
			//printf("%s already exists\n", simu_name);
		} else {
			//if (mkdir(mydir, S_IRWXU|S_IRWXG|S_IRGRP|S_IXGRP) != 0) {
			if (mkdir(simu_name, S_IRWXU | S_IRWXG | S_IRWXO) != 0) {
				cerr << "Error creating directory: %s\n" << endl;
				exit(1);
			}
		}
		//printf("%s successfully created\n", simu_name);
		//		bool initialdata = true;
		int i = 0;
		while (true){
			bool equilibrium = import_compaction (p, fin, back_bone, f_threshold);
			phi = M_PI*p.size()/(lx*lz);
			cerr << "n = "  << p.size() << endl;
			r_max = lz/2;
			if (equilibrium){
				cerr << "Equilibrium!" << endl;
				if (analyze_type == 'g'){
					evaluatePairCol_Compaction(simu_name, path, p,
											   resolution, r_min,
											   false,
											   back_bone, i);
				} else if (analyze_type == 'd'){
					evaluateDDC_Compaction(simu_name, path, p,
										   resolution, r_min,
										   z_periodic_boundary,
										   back_bone, i);
				} else if (analyze_type == 'b'){
					resolution = 50;
					evaluateBoxCounting(simu_name, path, p,
										resolution, r_min,
										z_periodic_boundary,
										back_bone, i);
				}
				////////////////////////////////////////////////////////////////////////
				i++;
			}
		}
	} else if (datatype == 's'){
		struct stat st;
		if (stat(simu_name, &st) == 0) {
			printf("%s already exists\n", simu_name);
		} else {
			//if (mkdir(mydir, S_IRWXU|S_IRWXG|S_IRGRP|S_IXGRP) != 0) {
			if (mkdir(simu_name, S_IRWXU | S_IRWXG | S_IRWXO) != 0) {
				cerr << "Error creating directory: %s\n" << endl;
				exit(1);
			}
		}
		printf("%s successfully created\n", simu_name);
		int i=0;
		while (true){
			double strain;
			bool equilibrium = import_shear(p, fin, strain);
			if (lz > 299){
				equilibrium = true;
			}
			if (equilibrium){
				double **ddc;
				ddc = new double* [resolution];
				for (int k=0; k < resolution; k++){
					ddc[k] = new double [2];
				}
//				r_max = lz/2;

				double area_fraction = (M_PI*p.size())/(lx*lz);
				densityDensityCorrelationFunction(r_min, resolution,
												  ddc, p, z_periodic_boundary);
				ofstream fout;
				ostringstream ddc_filename;

				ddc_filename << simu_name << "/ddc" << i << ".dat";
				
				fout.open(ddc_filename.str().c_str());
				fout << "# " << area_fraction << ' ' << strain << endl;
				for (int k=0; k< resolution; k++){
					fout << ddc[k][0] << ' ';
					fout << ddc[k][1] << endl;
				}
				fout.close();

				
				for (int k=0; k < resolution; k++){
					delete [] ddc[k];
				}
				delete [] ddc;
				
				i++;
			}
		}
	}
	
	
	
	//fin.close();
//	fout.close();
	
	return 0;
}

