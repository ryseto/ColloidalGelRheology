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

bool import_compaction(vector <vec3d> &p, ifstream &fin, double lx, double &lz){
	int num_particle;
	char buf[10];
	double info[13];
	char buf_line[512];
	fin >> buf ;
	cerr << buf << endl;
	if (buf[0] != '#'){
		cerr << "stop: " << buf << endl;
		exit(1);
	}
	fin >> buf;
	cerr << buf << endl;
	bool equilibrium = false;
	if (buf[0] == 'e')
		equilibrium  = true;
	int num_infomation = 10;
	for(int i=0;i<num_infomation; i++){
		fin >> info[i];
		cerr << info[i] << ' ';
	}
	cerr << endl;
	double z_bot = info[3];
	double z_top = info[4];
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

bool import_shear(vector <vec3d> &p, ifstream &fin, double lx, double lz, double &strain){
	int num_particle;
	char buf[10];
	double info[13];
	char buf_line[512];
	fin >> buf ;
	cerr << buf << endl;
	if (buf[0] != '#'){
		cerr << "stop: " << buf << endl;
		exit(1);
	}
	fin >> buf;
	cerr << buf << endl;
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

double sq_distancePeriodicBoundaries(vec3d p1, vec3d p2,double lx, double lz){
	if (abs(p1.x - p2.x) > lx/2){
		if (p2.x > p1.x)
			p2.x -= lx;
		else
			p2.x += lx;
	}
	
	if (abs(p1.z - p2.z) > lz/2){
		if (p2.z > p1.z)
			p2.z -= lz;
		else
			p2.z += lz;
	}
	return (p1 - p2).sq_norm_xz();
}


void densityDensityCorrelationFunction(
									   double r_min,
									   double r_max,
									   int resolution,
									   double **ddc,
									   vector <vec3d> &p,
									   double lx, double lz,
									   ofstream &fout){
	Grid2 grid;
	double cell_h = 2;
	unsigned long num_particle = p.size();
	double ly = 0;
	grid.init(num_particle, lx, ly, lz, cell_h);
	for (int i=0; i < num_particle; i++){
		grid.entry(p[i], i);
	}
	cerr << "------" << endl;
	srand(14235);
	if (r_max > lx/2)
		r_max = lx/2;
	double d_logr = (log(r_max)-log(r_min))/resolution;
	int  mmax = 100000;
	vec3d p1,p0;
	for (int k=0; k<resolution; k++){
		double r = r_min*exp(k*d_logr );
		cerr << "r = " << r << endl;
		int count_p0_inside = 0;
		int count_p1_inside = 0;
		for (int m=0; m < mmax; m++){
			p0.set(lx*drand48(), 0, lz*drand48());
			vector<int> neighbor;
			grid.get_neighbor_list(p0, neighbor);
			for (int i=0; i<neighbor.size(); i++){
				double sq_distance = sq_distancePeriodicBoundaries(p0, p[neighbor[i]], lx,lz);
				if (sq_distance <= 1){
					count_p0_inside ++;
					double theta = 2*M_PI*drand48();
					p1.x = p0.x + r*cos(theta);
					p1.z = p0.z + r*sin(theta);
					if (p1.x < 0){p1.x += lx;}
					else if (p1.x >= lx){p1.x -= lx;}
					if (p1.z < 0){ p1.z += lz;}
					else if (p1.z >= lz){  p1.z -= lz;  }
					vector<int> neighbor1;
					grid.get_neighbor_list(p1, neighbor1);
					for (int j=0; j<neighbor1.size(); j++){
						double sq_distance1 = sq_distancePeriodicBoundaries(p1, p[neighbor1[j]], lx,lz);
						if (sq_distance1 < 1){
							count_p1_inside ++;
							break;
						}
					}
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
	//srand48(1243);
	double lx, lz;
	double r_min = 1;
	double r_max;
	int resolution = 100;
	
	/* import_type
	 *  0: x z clusterid
	 *  1:
	 */
	int import_type;
	char simu_name[128];
	if (datatype =='i'){
		unsigned long i0 = path.find_first_of("2D_W*H*");
		if (i0 != 1){
			import_type = 0;
			unsigned long i1 = path.find_first_of("W") + 1;
			unsigned long i2 = path.find_first_of("H",i1);
			unsigned long i3 = path.find_first_of("_",i2);
			char width_str[4];
			char height_str[4];
			sprintf(width_str, "%s", (path.substr(i1,i2-i1)).c_str());
			sprintf(height_str, "%s", (path.substr(i2+1,i3-i2-1)).c_str());
			lx = atof(width_str);
			lz = atof(height_str);
			cerr << lx << ' ' << lz << endl;
		}
		
	} else if (datatype =='c'){
		lx = atof(argv[3]);
		import_type = 1;
		unsigned long i4 = path.find_first_of("conf_");
		unsigned long i5 = path.find_last_of(".dat");
		sprintf(simu_name, "%s", (path.substr(i4+5,(i5-3)-(i4+5))).c_str());
		
	} else if (datatype =='s'){
		import_type = 0;
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
	
	ofstream fout;
	string simu_name_str = simu_name;
	string outpuffile = "dcf_" + simu_name_str + ".dat";
	cerr << outpuffile << endl;
	fout.open( outpuffile.c_str());
	
	//if (! fin.is_open() ){
	//  cerr << "no initial file" << endl;
	//exit(1);
	//}
	//int num_particle;
	vector <vec3d> p;
	if (datatype == 'i'){
		double x, z;
		int tmp;
		while(!fin.eof()){
			fin >> x >> z >> tmp;
			p.push_back(vec3d(x,0,z));
		}
		p.pop_back();
		
		//densityDensityCorrelationFunction(p, lx, lz);
		//  densityDensityCorrelationFunction2(p, lx,lz, fout);
		
	} else if (datatype == 'c'){
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
		
		bool initialdata = true;
		int i = 0;
		while (true){
			bool equilibrium = import_compaction(p, fin, lx, lz);
			if (i==0){
				equilibrium = true;
				initialdata = false;
			}
			//vector< double[2]> ddc;
			r_max = lz/2;
			cerr << r_min << ' ' << r_max << endl;
			
			
			if (equilibrium){
				double **ddc;
				ddc = new double* [resolution];
				for (int k=0; k < resolution; k++){
					ddc[k] = new double [2];
				}
				densityDensityCorrelationFunction(r_min, r_max, resolution, ddc,
												  p, lx,lz, fout);
				
				double area_fraction = (M_PI*p.size())/(lx*lz);
				ofstream fout;
				ostringstream ddc_filename;
				ddc_filename << simu_name <<  "/ddc" << i << ".dat";
				fout.open(ddc_filename.str().c_str());
				fout << "# " << area_fraction << endl;
				for (int k=0; k< resolution; k++){
					fout << ddc[k][0] << ' ';
					fout << ddc[k][1] << endl;
					
				}
				fout.close();
				for (int k=0; k < resolution; k++){
					delete [] ddc[k];
				}
				delete [] *ddc;
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
												  ddc, 
												  p, lx,lz, fout);
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
	
	
	
	fin.close();
	fout.close();
	
	return 0;
}

