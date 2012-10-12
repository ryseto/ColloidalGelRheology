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

bool import_compaction(vector <vec3d> &p, ifstream &fin, double lx, double &lz,
					   bool back_bone){
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
//	exit(1);
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
		//particle_stress[p1] += 0.5*(abs(f_norm)+abs(f_slid));
		//		particle_stress[p2] += 0.5*(abs(f_norm)+abs(f_slid));
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
			if (particle_stress[i] > 0.05*wall_force ){
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
				//	cout << "@ 3\n";
				} else {
				//	cout << "@ 0\n";
				}
				//cout << "c " << p_import[i].x - lx/2 << " 0 " << p_import[i].z - lz/2 << endl;
			}
			//cout << endl;
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

bool import_shear(vector <vec3d> &p, ifstream &fin, double lx, double lz, double &strain){
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

double sq_distancePeriodicBoundaries(vec3d p1, vec3d p2,double lx, double lz, bool z_pd ){
	if (abs(p1.x - p2.x) > lx/2){
		if (p2.x > p1.x)
			p2.x -= lx;
		else
			p2.x += lx;
	}
	if (z_pd){
		if (abs(p1.z - p2.z) > lz/2){
			if (p2.z > p1.z)
				p2.z -= lz;
			else
				p2.z += lz;
		}
	}
	return (p1 - p2).sq_norm_xz();
}


void densityDensityCorrelationFunction(double r_min,
									   double r_max,
									   int resolution,
									   double **ddc,
									   vector <vec3d> &p,
									   double lx, double lz,
									   bool z_pd
									   ){
	Grid2 grid;
	double cell_h = 2;
	unsigned long num_particle = p.size();

	double ly = 0;
	grid.init(num_particle, lx, ly, lz, cell_h);
	for (int i=0; i < num_particle; i++){
		grid.entry(p[i], i);
	}
	srand(14235);
	if (r_max > lx/2){
		r_max = lx/2;
	}
	if (r_max > lz){
		r_max = lz;
	}
	
	double d_logr = (log(r_max)-log(r_min))/resolution;
	int  mmax = 100000;
	vec3d p1, p0;
	vector<int> neighbor;
	for (int k=0; k < resolution; k++){
		double r = r_min*exp(k*d_logr );
		int count_p0_inside = 0;
		int count_p1_inside = 0;
		while ( count_p0_inside < mmax ){
			p0.set(lx*drand48(), 0, lz*drand48());
			neighbor.clear();
			grid.get_neighbor_list(p0, neighbor);

			for (int i=0; i < neighbor.size(); i++){
				double sq_distance = sq_distancePeriodicBoundaries(p0, p[neighbor[i]], lx, lz, z_pd);

				if (sq_distance <= 1){
					count_p0_inside ++;
					if (z_pd){
						double theta = 2*M_PI*drand48();
						p1.x = p0.x + r*cos(theta);
						p1.z = p0.z + r*sin(theta);
						if (p1.z < 0){
							p1.z += lz;
						} else if (p1.z >= lz){
							p1.z -= lz;
						}
					} else {
						do{
							double theta = 2*M_PI*drand48();
							
							p1.x = p0.x + r*cos(theta);
							p1.z = p0.z + r*sin(theta);
							
						} while (p1.z > lz-1 || p1.z < +1);

					}
					if (p1.x < 0){
						p1.x += lx;
					}
					else if (p1.x >= lx){
						p1.x -= lx;
					}
//					cout << p1.x << ' ' << p1.z << endl;
					
					vector<int> neighbor1;
					grid.get_neighbor_list(p1, neighbor1);
					for (int j=0; j<neighbor1.size(); j++){
						double sq_distance1 = sq_distancePeriodicBoundaries(p1, p[neighbor1[j]], lx, lz, z_pd);
						if (sq_distance1 < 1){
							count_p1_inside ++;
							break;
						}
					}
					break;
				}
			}
		}
		//fout << r << ' ' << (double)count_p1_inside / count_p0_inside << ' ' << (double)count_p0_inside/mmax<< endl;
		ddc[k][0] = r ;
		ddc[k][1] = (double)count_p1_inside / count_p0_inside;
	}
}

int main(int argc, const char * argv[])
{
	char datatype = argv[1][0];
	string path = argv[2];
	double lx, lz;
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
	if (datatype =='i'){
		z_periodic_boundary = true;
		unsigned long i0 = path.find_first_of("2D_W*H*");
		if (i0 != 1){
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
		
	} else if (datatype =='c'){
		z_periodic_boundary = false;

		if (argc == 4 && argv[3][0]=='b'){
			back_bone = true;
		} else {
			back_bone = false ;
		}
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
		
	} else if (datatype =='s'){
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
		
	ifstream fin;
	fin.open( path.c_str());
	
	

	vector <vec3d> p;
	if (datatype == 'i'){
		ofstream fout;
		string simu_name_str = simu_name;
		string outpuffile = "dcf_" + simu_name_str + ".dat";
		cerr << outpuffile << endl;
		
		fout.open( outpuffile.c_str());
		double x, z;
		int tmp;
		while(!fin.eof()){
			fin >> x >> z >> tmp;
			cerr << x << ' ' << z << endl;
			p.push_back(vec3d(x,0,z));
		}
		p.pop_back();
		
		phi = M_PI*p.size()/(lx*lz);
		cerr << "phi = " << phi << endl;
		
		double **ddc;
		ddc = new double* [resolution];
		for (int k=0; k < resolution; k++){
			ddc[k] = new double [2];
		}
		r_max = lz/2;

		densityDensityCorrelationFunction(r_min, r_max,
										  resolution, ddc,
										  p, lx, lz, z_periodic_boundary);
		
		if(fout.is_open()){
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
		// densityDensityCorrelationFunction2(p, lx,lz, fout);
		fout.close();
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
			bool equilibrium = import_compaction(p, fin, lx, lz, back_bone);
			phi = M_PI*p.size()/(lx*lz);
			cerr << "n = "  << p.size() << endl;
			//			if (i == 0){
			//				equilibrium = true;
			//			}
			r_max = lz/2;
			if (equilibrium){
				cerr << "Equilibrium!" << endl;
				double **ddc;
				ddc = new double* [resolution];
				for (int k=0; k < resolution; k++){
					ddc[k] = new double [2];
				}
				densityDensityCorrelationFunction(r_min, r_max, resolution, ddc,
												  p, lx, lz, z_periodic_boundary);

				double area_fraction = (M_PI*p.size())/(lx*lz);
				ofstream fout;
				ostringstream ddc_filename;
				if (back_bone){
					ddc_filename << simu_name <<  "/bb_ddc" << i << ".dat";
				} else {
					ddc_filename << simu_name <<  "/ddc" << i << ".dat";
					
				}
				fout.open(ddc_filename.str().c_str());
				fout << "# " << area_fraction << endl;
				for (int k=0; k< resolution; k++){
					fout << ddc[k][0] << ' ';
					fout << ddc[k][1]/phi << endl;
				}
				fout.close();
				for (int k=0; k < resolution; k++){
					delete [] ddc[k];
				}
				delete [] ddc;
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
			bool equilibrium = import_shear(p, fin, lx, lz, strain);
			if (lz > 299){
				equilibrium = true;
			}
			if (equilibrium){
				double **ddc;
				ddc = new double* [resolution];
				for (int k=0; k < resolution; k++){
					ddc[k] = new double [2];
				}
				r_max = lz/2;
				
				double area_fraction = (M_PI*p.size())/(lx*lz);
				densityDensityCorrelationFunction(r_min, r_max, resolution, 
												  ddc, p, lx,lz, z_periodic_boundary);
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

