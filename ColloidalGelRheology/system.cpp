//
//  system.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

//#include <iostream>
#include "system.h"

System::System()
{
	n_particle = 0;
	n_bond = 0;
	version = "0";
	ct = new ConTable;
	grid = new Grid;
}

System::~System()
{
	fout_conf.close();
	fout_data.close();
	fout_log.close();
	fout_deform.close();
	delete ct;
	delete grid;
}

void System::affineUniaxialCompression(){
	double lz_previous = lz;
	lz -= compaction_speed*dt;
	double z_compaction = lz/lz_previous;
	lz_half = lz/2;
	foreach (vector<Particle *>, particle, it_particle) {
		(*it_particle)->p.z *= z_compaction;
	}
}

void System::affineBiaxialCompression(){
	double lx_previous = lx;
	double lz_previous = lz;
	lx -= compaction_speed*dt;
	lx_half = lx/2;
	lz -= compaction_speed*dt;
	lz_half = lz/2;
	double z_compaction = lz/lz_previous;
	double x_compaction = lx/lx_previous;
	foreach (vector<Particle *>, particle, it_particle) {
		(*it_particle)->p.x *= x_compaction;
		(*it_particle)->p.z *= z_compaction;
	}
}

void System::bounadyUniaxialCompression(){
	lz -= compaction_speed*dt;
	lz_half = lz/2;
}

void System::bounadyBiaxialCompression(){
	lx -= compaction_speed*dt;
	lx_half = lx/2;
	lz -= compaction_speed*dt;
	lz_half = lz/2;
}

void System::TimeDevPeriodicBoundaryCompactionEuler()
{
	if (prog_strain) {
		(this->*enforceStrain)();
	}
 	foreach (vector<Bond *>, bond, it_bond) {
		(*it_bond)->addContactForce();
	}
	//////////////////////////////
	// calculate
	// x*_{j+1} and v*_{j+1}
	foreach (vector<Particle *>, particle, it_particle) {
		(*it_particle)->move_Euler();
	}
	time += dt;
}

void System::TimeDevBendingTest()
{
	particle[0]->stackForce(-0.5*f_ex, 0);
	particle[5]->stackForce(f_ex, 0);
	particle[10]->stackForce(-0.5*f_ex, 0);
 	foreach (vector<Bond *>, bond, it_bond) {
		(*it_bond)->addContactForce();
	}
	//////////////////////////////
	// calculate
	// x*_{j+1} and v*_{j+1}
	foreach (vector<Particle *>, particle, it_particle) {
		(*it_particle)->move_Euler();
	}
	time += dt;
}

bool System::mechanicalEquilibrium()
{
	static int cnt_check = 0;
	bool mechanical_equilibrium = false;
	if (simulation == "s"){
		if ((max_velocity < max_velocity_convergence // less important
			 && max_ang_velocity < max_ang_velocity_convergence // less important
			 && counterRegenerate_before == counterRegenerate)) {
			mechanical_equilibrium = true;
		}
	} else {
		max_displacement = 0;
		for (int i=0; i<n_particle; i++) {
			vec3d disp = particle[i]->p-pos_previous[i];
			if (disp.x > lx_half) {
				disp.x -= lx;
			} else if (disp.x < -lx_half) {
				disp.x += lx;
			}
#ifndef TWODIMENSION
			if (disp.y > ly_half) {
				disp.y -= ly;
			} else if (disp.y < -ly_half) {
				disp.y += ly;
			}
#endif
			if (disp.z > lz_half) {
				disp.z -= lz;
			} else if (disp.z < -lz_half) {
				disp.z += lz;
			}
			if (disp.norm() > max_displacement) {
				max_displacement = disp.norm();
			}
		}
		cerr << "cnt_check : " << cnt_check << endl;
		cerr << "counterRegenerate : " << counterRegenerate_before << ' ' << counterRegenerate << endl;
		cerr << "max_velocity : " << max_velocity << endl;
		cerr << "max_ang_velocity : " << max_ang_velocity << endl;
		cerr << "max_displacement : " << max_displacement << endl;
		if (max_displacement > 1) {
			mechanical_equilibrium = false;
			max_displacement = 0;
		} else if (cnt_check > 0
				   && counterRegenerate_before == counterRegenerate
				   && max_velocity < max_velocity_convergence
				   && max_ang_velocity < max_ang_velocity_convergence
				   && max_displacement < eq_max_displacement) {
			mechanical_equilibrium = true;
		} else {
			mechanical_equilibrium = false;
		}
	}
	counterRegenerate_before = counterRegenerate;
	for (int i=0; i<n_particle; i++) {
		pos_previous[i] = particle[i]->p;
	}
	cnt_check++;
	if (mechanical_equilibrium == true) {
		cnt_check = 0;
	}
	return mechanical_equilibrium;
}

void System::preProcesses()
{
	initialprocess = true;
	setSimulationViscosity();
	preparationOutput();
	setBondGenerationDistance(2.0);
	initGrid();
	initDEM();
	makeInitialBond(2.0001);
	initialprocess = false;
	calcVolumeFraction();
	dt = dt_max;
	checkState();
	outputData();
	pos_previous.resize(n_particle);
	for (int i=0; i<n_particle; i++) {
		pos_previous[i] = particle[i]->p;
	}
	cnt_output_config = 0;
}

void System::timeEvolution()
{
	double t_next = time+interval_convergence_check*dt_max;
	int calc_count = 0;
	counter_relax_for_restructuring = 0;
	while (time < t_next) {
		prog_strain = (volume_fraction < vf_target);
		/* Main time evolution */
		if (simulation == "s") {
			/* Shear strain is not yet implemented.
			 * TimeDevStrainControlShearEuler();
			 */
		} else {
			TimeDevPeriodicBoundaryCompactionEuler();
		}
		/* Breakup of bonds */
		if (counter_relax_for_restructuring++ >= relax_for_restructuring) {
			checkBondFailure();
			if (!regeneration_bond.empty()) {
				regeneration();
				counter_relax_for_restructuring = 0;
			}
			if (!rupture_bond.empty()) {
				rupture();
				counter_relax_for_restructuring = 0;
			}
		}
		/* Bond generation */
		if (calc_count++ % interval_makeNeighbor == 0) {
			double cell_h = 2.2;
			grid->init_z(n_particle, lx, ly, lz, cell_h);
			//grid->init(n_particle, lx, ly, lz, cell_h);
			makeNeighborPB();
		}
		generateBond();
	}
	if (simulation != "s") {
		calcVolumeFraction();
	}
}

void System::middleProcedures()
{
	simuAdjustment();
	calcStress();
	checkState();
	makeNeighborPB();
	optimalTimeStep();
	output_log();
	cnt_output_config++;
}

void System::strainControlSimulation()
{
	preProcesses();
	cnt_loop = 0;
	max_displacement = 0;
	if (simulation == "c1") {
		enforceStrain = &System::affineUniaxialCompression;
		
	} else if (simulation == "c2") {
		enforceStrain = &System::affineBiaxialCompression;
	} else {
		exit(1);
	}
	if (simulation == "c1" || simulation == "c2") {
		interval_output_config = (int)(2./(compaction_speed*dt))/interval_convergence_check;
	} else {
		interval_output_config = 20;
	}
	int cnt = 0;
	outputConfiguration('e');
	vf_target = volume_fraction;
	setTarget();
	while (true) {
		timeEvolution();
		middleProcedures();
		outlog();
 		if (cnt ++ % 50 == 0) {
			outputConfiguration('n');
		}
		if (volume_fraction >= vf_target) {
			if (mechanicalEquilibrium()) {
				cerr <<"Equilibrium!"<< endl;
				calcStress();
				/*  Output */
				outputData();
				outputConfiguration('e');
				monitorDeformation('e');
				/*  Prepare next step */
				setTarget();
				if (checkEndCondition()) {
					break;
				}
			}
		}
		cnt_loop ++;
	}
	return;
}

void System::bendingSimulation()
{
	initialprocess = true;
	setSimulationViscosity();
	preparationOutput();
	setBondGenerationDistance(2);
	initGrid();
	initDEM();
	makeInitialBond(2.0002);
	initialprocess = false;
	dt = dt_max;
	checkState();
	int cnt = 0 ;
	double time_interval = 300;
	double time_max = time_interval;
	for (int i=0; i<2 ; i++) {
		if (i == 1) {
			f_ex.set(0,0,0);
		}
		while(true) {
			TimeDevBendingTest();
			checkBondFailure();
			if (!regeneration_bond.empty()){
				regeneration();
				counter_relax_for_restructuring = 0;
			}
			if (!rupture_bond.empty()){
				rupture();
				counter_relax_for_restructuring = 0;
			}
			if (cnt++ % 100 ==0){
				
				fout_conf << "# " << f_ex.x << endl;
				foreach( vector<Particle *>, particle, it_particle){
					cout << "c ";
					cout << (*it_particle)->p.x -0.5*lx << ' ';
					cout << (*it_particle)->p.y -0.5*ly << ' ';
					cout << (*it_particle)->p.z -0.5*lz << endl;
					fout_conf << (*it_particle)->p.x -0.5*lx << ' ';
					fout_conf << (*it_particle)->p.y -0.5*ly << ' ';
					fout_conf << (*it_particle)->p.z -0.5*lz << ' ';
					fout_conf << (*it_particle)->orientation.q[0] << ' ';
					fout_conf << (*it_particle)->orientation.q[1] << ' ';
					fout_conf << (*it_particle)->orientation.q[2] << ' ';
					fout_conf << (*it_particle)->orientation.q[3] << endl;
				}
				cout << "l 0 -0.01 -12 0 -0.01 12 \n";
				cout << "l -12 -0.01 0 12 -0.01 0 \n";
				cout << endl;
				double diff_x = particle[5]->p.x-particle[0]->p.x;
				fout_data << time << ' ' << diff_x << ' ' << f_ex.x << ' ' << counterRegenerate << endl;
			}
			if (time > time_max) {
				break;
			}
		}
		time_max +=time_interval;
	}
	return;
}

void System::outlog()
{
	if (simulation == "s") {
		cerr << cnt_output_config << ' ' << shear_strain << ": " << strain_target  << endl;
	} else {
		cerr << cnt_output_config << ' ' << volume_fraction << ": " << vf_target  << endl;
	}
}

void System::setTarget()
{
	if (simulation == "s") {
		strain_target += step_shear_strain;
	} else {
		vf_target *= volumefraction_increment;
	}
}

bool System::checkEndCondition()
{
	bool end_simulation = false;
	if (simulation == "s") {
		if (shear_strain > max_shear_strain) {
			end_simulation = true;
		}
	} else {
		if (volume_fraction > max_volume_fraction) {
			end_simulation = true;
		}
	}
	return end_simulation;
}

bool System::checkOutputStrain(double diff)
{
	if (simulation == "s") {
		if (shear_strain > strain_x_last_output+diff) {
			lz_last_output = lz;
			return true;
		}
	} else {
		if (lz < lz_last_output-diff) {
			lz_last_output = lz;
			return true;
		}
	}
	return false;
}

void System::setParameterFile(char *parameters_file_)
{
	sprintf(parameters_file, "%s", parameters_file_);
	string s_parameters_file = parameters_file;
	unsigned long i_backslash = s_parameters_file.find_last_of( "/") + 1;
	unsigned long i_extention = s_parameters_file.find( ".txt" );
	sprintf(parameters, "%s",
			(s_parameters_file.substr( i_backslash, i_extention-i_backslash)).c_str());
	cerr << parameters_file << ' ' << parameters << endl;
}

void System::initDEM()
{
	ct->set(n_particle);
	rup_normal = rup_shear = rup_bend = rup_torsion = 0;
	rupture_bond.clear();
	regeneration_bond.clear();
	time = 0.;
	counterBreak = 0;
	counterRegenerate = 0;
	shear_strain = 0;
	strain_z = 0;
	counter_relax_for_restructuring = 0;
	lz_last_output = lz_init;
	counterRegenerate_before = 0;
}

void System::preparationOutput()
{
	char fn_snapshot[128];
	char fn_data[128];
	char fn_log[128];
	char fn_conf[128];
	char fn_deform[128];
	char fn_particle[128];
	char fn_bond[128];
	char name_string[128];
	if (simulation == "c1" || simulation == "c2" || simulation == "s") {
		sprintf(name_string, "%s_%s_%s_%s",
				parameters, simulation.c_str(),
				init_cluster, version.c_str());
	} else if (simulation == "bt") {
		sprintf(name_string, "%s_%1.2f", parameters, f_ex.x );
	} else {
		cerr << "Not yet implemented" << endl;
		exit(1);
	}
	string fn_common = (string)name_string;
	cerr << version << endl;
	cerr << fn_common.c_str() << endl;
	sprintf(fn_data, "data_%s.dat", fn_common.c_str());
	sprintf(fn_log, "log_%s.dat", fn_common.c_str());
	sprintf(fn_conf, "conf_%s.dat", fn_common.c_str());
	sprintf(fn_deform, "deform_%s.dat", fn_common.c_str());
	sprintf(fn_particle, "par_%s.dat", fn_common.c_str());
	sprintf(fn_bond, "bond_%s.dat", fn_common.c_str());
	if (simulation == "bt") {
		sprintf(fn_snapshot, "ss_%s.yap", fn_common.c_str());
		fout_yap.open(fn_snapshot);
	}
	fout_data.open(fn_data);
	fout_log.open(fn_log);
	fout_conf.open(fn_conf);
	fout_deform.open(fn_deform);
	fout_particle.open(fn_particle);
	fout_bond.open(fn_bond);
}

void System::setInitClusterFile(string init_cluster_file_)
{
	sprintf(init_cluster_file, "%s", init_cluster_file_.c_str());
	string s_init_cluster_file = init_cluster_file;
	unsigned long i_backslash = s_init_cluster_file.find_last_of( "/") + 1;
	unsigned long i_extention = s_init_cluster_file.find( ".dat" );
	sprintf(init_cluster, "%s",
			(s_init_cluster_file.substr(i_backslash, i_extention-i_backslash)).c_str());
}

/*
 * to distinguish the output files :
 * file name is FileName_***.dat
 * "***" is given by command argument.
 * Default: version = 0
 */
void System::setVersion(string version_)
{
	//sprintf(version, "%s", version_);
	version = version_;
}

void System::setBondGenerationDistance(double distance)
{
	dist_generate = distance;
	sq_dist_generate = sq(distance);
}

void System::generateBond()
{
	for (int i=0; i<n_particle; i++) {
		particle[i]->generateBond();
	}
	return;
}

void System::setSimulationViscosity()
{
	eta = eta_factor*2*sqrt(bond0.kb);
	eta_rot = (4./3)*eta;
}

void System::optimalTimeStep()
{
	double dt_new1 = max_move_step/max_velocity;
	double dt_new2 = max_move_step/max_ang_velocity;
	double dt_new;
	if (dt_new1 < dt_new2) {
		dt_new = dt_new1;
	} else {
		dt_new = dt_new2;
	}
	if (dt_new > dt_max) {
		dt_new = dt_max;
	} else {
		cerr << "small time step"  << dt_new << endl;
	}
	dt = dt_new;
}

void System::calcStress()
{
	total_contact_stressXF.reset();
	foreach(vector<Bond *>, bond, it_bond) {
		(*it_bond)->calcContactStress();
		total_contact_stressXF += (*it_bond)->getContactStressXF();
	}
	total_contact_stressXF /= system_volume();
}

void prepareBond(BondParameter & _bond)
{
	_bond.n_max = _bond.fnc/_bond.kn;
	_bond.s_max = _bond.fsc/_bond.ks;
	_bond.b_max = _bond.mbc/_bond.kb;
	_bond.t_max = _bond.mtc/_bond.kt;
	cerr << "k:" << _bond.kn << ' ' << _bond.ks << ' ' << _bond.kb << ' ' << _bond.kt << endl;
	cerr << "fc:" << _bond.fnc << ' '  << _bond.fsc << ' ' << _bond.mbc << ' ' << _bond.mtc << endl;
	cerr << "critical_strain:" << _bond.n_max ;
	cerr << " * " << _bond.s_max;
	cerr << ' ' << _bond.b_max << ' ' << _bond.t_max << endl;
}

void readBond(const string &codeword,
			  const string &value,
			  BondParameter &bondparameter)
{
	map<string,int> keylist;
	keylist["fnc:"]=1; const int _fnc = 1;
	keylist["fsc:"]=2; const int _fsc = 2;
	keylist["mbc:"]=3; const int _mbc = 3;
	keylist["mtc:"]=4; const int _mtc = 4;
	keylist["kn:"]=5; const int _kn = 5;
	keylist["ks:"]=6; const int _ks = 6;
	keylist["kb:"]=7; const int _kb = 7;
	keylist["kt:"]=8; const int _kt = 8;
	keylist["kn3:"]=9; const int _kn3 = 9;
	cerr << codeword << ' ' << value << endl;
	switch (keylist[codeword]) {
		case _fnc: bondparameter.fnc = atof(value.c_str()); break;
		case _fsc: bondparameter.fsc = atof(value.c_str()); break;
		case _mbc: bondparameter.mbc = atof(value.c_str()); break;
		case _mtc: bondparameter.mtc = atof(value.c_str()); break;
		case _kn: bondparameter.kn = atof(value.c_str()); break;
		case _ks: bondparameter.ks = atof(value.c_str()); break;
		case _kb: bondparameter.kb = atof(value.c_str()); break;
		case _kt: bondparameter.kt = atof(value.c_str()); break;
		case _kn3: bondparameter.kn3 = atof(value.c_str()); break;
		default:
			cerr << "Bond: The codeword " << codeword << " is'nt associated with an parameter" << endl;
			exit(1);
	}
}

void setBondParameter(BondParameter &bondparameter, string &bond_file)
{
	string path_bond = bond_file;
	ifstream fin;
	if (fin.is_open()) {
		fin.close();
	}
	fin.open(path_bond.c_str());
	string codeword;
	string value;
	while (!fin.eof()) {
		fin >> codeword;
		if (codeword == "#") {
			char buf[1024]; fin.get(buf, 1024);
		} else if (codeword == "!") {
			break;
		} else {
			fin >> value ;
			if (fin.eof()) break;
			readBond(codeword, value, bondparameter);
		}
	}
	fin.close();
	prepareBond(bondparameter);
	return;
}

void System::readBondParameter()
{
	setBondParameter(bond0, bond0_file);
	setBondParameter(bond1, bond1_file);
}

void System::readParameterFile()
{
	//////////////////////////////////////////
	// input arguments:
	//////////////////////////////////////////
	// Read parameter file
	ifstream fin_parameter;
	fin_parameter.open(parameters_file);
	if (fin_parameter.is_open()) {
		cerr << parameters_file << " is opend." << endl;
	}
	string codeword;
	string value;
	while (!fin_parameter.eof()){
		fin_parameter >> codeword;
		if ( codeword == "#") {
			char buf[1024]; fin_parameter.get(buf, 1024);
		} else if (codeword == "!") {
			break;
		} else {
			fin_parameter >> value;
			if (fin_parameter.eof()) {
				break;
			}
			//set_value(codeword, value, unit);
			readParameter(codeword, value);
		}
	}
	fin_parameter.close();
	if (!fin_parameter.is_open()) {
		cerr << parameters_file << " is closed." << endl;
	}
}

void System::readParameter(const string &codeword, const string &value)
{
	map<string,int> keylist;
	keylist["bond0_file:"] = 2; const int _bond0_file = 2;
	keylist["bond1_file:"] = 3; const int _bond1_file = 3;
	keylist["volumefraction_increment:"] = 12; const int _volumefraction_increment = 12;
	keylist["max_volume_fraction:"] = 13; const int _max_volume_fraction = 13;
	keylist["step_shear_strain:"] = 14; const int _step_shear_strain = 14;
	keylist["max_shear_strain:"] = 15; const int _max_shear_strain = 15;
	keylist["compaction_speed:"] = 20;const int _compaction_speed = 20;
	keylist["dt_max:"] = 21; const int _dt_max = 21;
	keylist["max_move_step:"] = 22; const int _max_move_step = 22;
	keylist["eta_factor:"] = 23;const int _eta_factor = 23;
	keylist["relax_for_restructuring:"] = 30; const int _relax_for_restructuring = 30;
	keylist["interval_convergence_check:"] = 31; const int _interval_convergence_check = 31;
	keylist["interval_makeNeighbor:"] = 32; const int _interval_makeNeighbor = 32;
	keylist["max_velocity_convergence:"] = 33; const int _max_velocity_convergence = 33;
	keylist["max_ang_velocity_convergence:"] = 34; const int _max_ang_velocity_convergence = 34;
	keylist["diff_stress_convergence:"] = 35;const int _diff_stress_convergence = 35;
	keylist["stress_change_convergence:"] = 36;const int _stress_change_convergence = 36;
	keylist["stress_minimum:"] = 37; const int _stress_minimum = 37;
	keylist["min_relaxation_loop:"] = 38; const int _min_relaxation_loop = 38;
	keylist["max_relaxation_loop:"] = 39; const int _max_relaxation_loop = 39;
	keylist["eq_max_displacement:"] = 40; const int _eq_max_displacement = 40;
	cerr << codeword << ' ' << value << endl;
	switch (keylist[codeword]) {
		case _bond0_file: bond0_file = value; break;
		case _bond1_file: bond1_file = value; break;
		case _interval_convergence_check: interval_convergence_check = atoi(value.c_str()); break;
		case _interval_makeNeighbor: interval_makeNeighbor = atoi(value.c_str()); break;
		case _max_velocity_convergence: max_velocity_convergence = atof(value.c_str()); break;
		case _compaction_speed: compaction_speed = atof(value.c_str()); break;
		case _dt_max: dt_max = atof(value.c_str()); break;
		case _eta_factor: eta_factor = atof(value.c_str()); break;
		case _max_ang_velocity_convergence: max_ang_velocity_convergence = atof(value.c_str()); break;
		case _max_volume_fraction: max_volume_fraction = atof(value.c_str()); break;
		case _step_shear_strain: step_shear_strain = atof(value.c_str()); break;
		case _max_shear_strain: max_shear_strain = atof(value.c_str()); break;
		case _max_move_step: max_move_step = atof(value.c_str()); break;
		case _volumefraction_increment: volumefraction_increment = atof(value.c_str()); break;
		case _diff_stress_convergence: diff_stress_convergence = atof(value.c_str()); break;
		case _stress_change_convergence: stress_change_convergence = atof(value.c_str()); break;
		case _relax_for_restructuring: relax_for_restructuring = atoi(value.c_str()); break;
		case _stress_minimum: stress_minimum = atof(value.c_str()); break;
		case _min_relaxation_loop: min_relaxation_loop = atoi(value.c_str()); break;
		case _max_relaxation_loop: max_relaxation_loop = atoi(value.c_str()); break;
		case _eq_max_displacement: eq_max_displacement = atof(value.c_str()); break;
		default:
			cerr << "The codeword " << codeword << " is'nt associated with an parameter" << endl;
			exit(1);
	}
}

void System::importPositions()
{
	vector<vec3d> init_aggregate;
	vector<int> init_aggregate_cluster;
	ifstream fin;
	string path = init_cluster_file;
	string s_parameters_file = parameters_file;
	unsigned long i0 = path.find_first_of("D");
	unsigned long i1 = path.find_first_of("W",i0)+1;
	unsigned long i2 = path.find_first_of("H",i1);
	unsigned long i3 = path.find_first_of("_",i2);
	char width_str[4];
	char height_str[4];
	sprintf(width_str, "%s", (path.substr(i1,i2-i1)).c_str());
	sprintf(height_str, "%s", (path.substr(i2+1,i3-i2-1)).c_str());
	fin.open(path.c_str());
	if (!fin.is_open()) {
		cerr << "no initial file" << endl;
		exit(1);
	}
	lx = atof(width_str);
	lz_init = atof(height_str);
	lz = lz_init;
#ifdef TWODIMENSION
	ly = 0;
	ly_half = 0;
#else
	ly = lx;
	ly_half = ly/2;
#endif
	lx_half = lx/2;
	lz_half = lz/2;
	double x, y, z;
	int i_cluster;
	if (simulation == "c1"
		|| simulation == "c2"
		|| simulation == "s") {
		cerr << "AA" << endl;
		do {
#ifdef TWODIMENSION
			fin >> x >> z >> i_cluster ;
#else
			fin >> x >> y >> z;
#endif
#ifdef TWODIMENSION
			vec3d new_p(x, y, z);
#else
			vec3d new_p(x, 0, z);
#endif
			init_aggregate.push_back(new_p);
			init_aggregate_cluster.push_back(i_cluster);
			
		} while (!fin.eof());
	} else {
		cerr << simulation << endl;
		exit(1);
	}
	init_aggregate.pop_back();
	int count_particle_number = 0;
	cerr << "N=" << init_aggregate.size() << endl;
	for (int i = 0; i < init_aggregate.size(); i++) {
		vec3d new_p = init_aggregate[i];
		particle.push_back(new Particle(count_particle_number, new_p,
										init_aggregate_cluster[i],
										*this) );
		count_particle_number++;
	}
	n_particle = count_particle_number;
	if (n_particle == 0) {
		cerr << "failed to import the initial configuration.\n";
		exit(1);
	}
	volume_fraction_init = n_particle*M_PI/(lx*lz_init);
	cerr << "volume_fraction_init = " << volume_fraction_init << endl;
}

void System::setLinearChain(int number_of_particles)
{
	int m = number_of_particles;
	int count_particle_number = 0;
	for (int i=0; i<m; i++) {
		vec3d new_p(0.5*lx, 0, 2*i-number_of_particles+1+0.5*lz);
		particle.push_back(new Particle(count_particle_number, new_p, 0, *this));
		++ count_particle_number;
	}
	n_particle = count_particle_number;
	if (n_particle == 0) {
		cerr << "failed to import the initial configuration.\n";
		exit(1);
	}
}

void System::shiftCenterOfMass(vector<vec3d> &p)
{
	vec3d cm(0,0,0);
	foreach(vector<vec3d>, p, p_iter) {
		cm += (*p_iter);
	}
	cm *= 1.0/p.size();
	foreach(vector<vec3d>, p, p_iter) {
		*p_iter -= cm;
		p_iter->x += lx/2;
#ifdef TWODIMENSION
		p_iter->y += ly/2;
#endif
		p_iter->z += lz/2;
	}
}

void System::calcLocalStrains()
{
	/* The minimum and maximum distances between two particles.
	 */
	r_min = 2;
	r_max = 2;
	bending_angle_max = 0;
	torsional_angle_max = 0;
	sliding_disp_max = 0;
	for (int i=0; i < n_bond; i++) {
		if (bond[i]->status) {
			if (bond[i]->r < r_min) {
				r_min = bond[i]->r;
			}
			if (bond[i]->r > r_max) {
				r_max = bond[i]->r;
			}
			if (bond[i]->valSlidingDisplacement() > sliding_disp_max) {
				sliding_disp_max = bond[i]->valSlidingDisplacement();
			}
			if (bond[i]->valBendingAngle() > bending_angle_max) {
				bending_angle_max = bond[i]->valBendingAngle();
			}
#ifndef TWODIMENSION
			if (bond[i]->valTorsionalAngle() > torsional_angle_max) {
				torsional_angle_max = bond[i]->valTorsionalAngle();
			}
#endif
		}
	}
	return;
}

void System::initGrid()
{
	if (n_particle == 0) {
		cerr << "before setting system box, ";
		cerr << "it is required to import configuration of particles.";
		cerr << endl;
		exit(1);
	}
	double cell_h = 2.2;
	grid->init(n_particle, lx, ly, lz, cell_h);
}

void System::outputData()
{
	static bool firsttime = true;
	static int rup_normal_previous = 0;
	static int rup_shear_previous = 0;
	static int rup_bend_previous = 0;
	static int rup_torsion_previous = 0;
	static double strain_previous = 0;
	double del_strain = 0;
	if (firsttime) {
		fout_data << "#1 volume_fraction\n";
		fout_data << "#2 shear_strain\n";
		fout_data << "#3 ParticlePressure\n";
		fout_data << "#4 shear_stress\n";
		fout_data << "#5 compressive_stress(sigma_22)\n";
		fout_data << "#6 lateral_stress\n";
		fout_data << "#7 normal_stress_diff_1\n";
		fout_data << "#8 normal_stress_diff_2\n";
		fout_data << "#9  rup_normal(total number)\n";
		fout_data << "#10 rup_shear(total number)\n";
		fout_data << "#11 rup_bend(total number)\n";
		fout_data << "#12 rup_torsion(total number)\n";
		fout_data << "#13 rate_rup_normal\n";
		fout_data << "#14 rate_rup_shear\n";
		fout_data << "#15 rate_rup_bend\n";
		fout_data << "#16 rate_rup_torsion\n";
		fout_data << "#17 average_contact_number\n";
		fout_data << "#18 del_strain\n";
		fout_data << "# volume_fraction " << volume_fraction << endl;
		fout_data << "# number_of_particle " << particle.size() << endl;
		fout_data << "# lx " << lx << endl;
		fout_data << "# ly " << ly << endl;
		fout_data << "# lz " << lz << endl;
	}
	if (simulation == "c1" || simulation == "c2") {
		del_strain = (strain_z - strain_previous);
		strain_previous = strain_z;
	} else if (simulation == "s") {
		del_strain = (shear_strain - strain_previous);
		strain_previous = shear_strain;
	}
	double rate_rup_normal, rate_rup_shear, rate_rup_bend, rate_rup_torsion;
	if (del_strain == 0) {
		rate_rup_normal = 0;
		rate_rup_shear = 0;
		rate_rup_bend = 0;
		rate_rup_torsion = 0;
	} else {
		rate_rup_normal = (rup_normal-rup_normal_previous)/del_strain;
		rate_rup_shear = (rup_shear-rup_shear_previous)/del_strain;
		rate_rup_bend = (rup_bend-rup_bend_previous)/del_strain;
		rate_rup_torsion = (rup_torsion-rup_torsion_previous)/del_strain;
	}
	if (firsttime) {
		firsttime = false;
		return;
	}
	shear_strain = 0;
	fout_data << volume_fraction << ' '; // 1
	fout_data << shear_strain << ' '; // 2
	fout_data << total_contact_stressXF.getParticlePressure() << ' '; // 3
	fout_data << total_contact_stressXF.getShearStress() << ' '; // 4
	fout_data << total_contact_stressXF.getCompressiveStress() << ' ' ; // 5
	fout_data << total_contact_stressXF.getLateralStress() << ' '; // 6
	fout_data << total_contact_stressXF.getNormalStress1() << ' '; // 7
	fout_data << total_contact_stressXF.getNormalStress2() << ' '; // 8
	///////////////////////////////////////
	fout_data << rup_normal << ' '; //9
	fout_data << rup_shear << ' ';  //10
	fout_data << rup_bend << ' ';   //11
	fout_data << rup_torsion << ' ';//12
	///////////////////////////////////////
	fout_data << rate_rup_normal << ' '; //13
	fout_data << rate_rup_shear << ' ';  //14
	fout_data << rate_rup_bend << ' ';   //15
	fout_data << rate_rup_torsion << ' ';//16
	///////////////////////////////////////
	fout_data << average_contact_number << ' '; //17
	fout_data << del_strain << ' '; // 18
	fout_data << endl;
	rup_normal_previous = rup_normal;
	rup_shear_previous = rup_shear;
	rup_bend_previous = rup_bend;
	rup_torsion_previous = rup_torsion;
}

void System::output_log()
{
	static bool firsttime = true;
	if (firsttime) {
		firsttime = false;
		fout_log << "#1 time\n";
		fout_log << "#2 volume_fraction\n";
		fout_log << "#3 shear_strain\n";
		fout_log << "#4 ParticlePressure\n";
		fout_log << "#5 ShearStress\n";
		fout_log << "#6 CompressiveStress(sigma_22)\n";
		fout_log << "#7 LateralStress\n";
		fout_log << "#8 NormalStress1\n";
		fout_log << "#9 NormalStress2\n";
		fout_log << "#10 lz\n";
		fout_log << "#11 lx\n";
		fout_log << "#12 average_contact_number\n";
		fout_log << "#13 ave_force\n";
		fout_log << "#14 max_force\n";
		fout_log << "#15 r_min\n";
		fout_log << "#16 r_max\n";
		fout_log << "#17 sliding_disp_max\n";
		fout_log << "#18 bending_angle_max\n";
		fout_log << "#19 torsional_angle_max\n";
		fout_log << "#17 counterRegenerate\n";
		fout_log << "#18 counterBreak\n";
		fout_log << "#19 bond number\n";
		fout_log << "#20 max_move_step/max_velocity\n";
		fout_log << "#21 counterRegenerate\n";
		fout_log << "#22 counterBreak\n";
		fout_log << "#23 numBond\n";
		fout_log << "#24 max_velocity\n";
		fout_log << "#25 max_ang_velocity\n";
		fout_log << "#26 ave_bondforce\n";
		fout_log << "#27 cnt_loop\n";
		fout_log << "#28 max_displacement\n";
		fout_log << "#29 prog_strain\n";
	}
	fout_log << time << ' '; // 1
	fout_log << volume_fraction << ' '; // 2
	fout_log << shear_strain << ' '; // 3
	fout_log << total_contact_stressXF.getParticlePressure() << ' '; // 4
	fout_log << total_contact_stressXF.getShearStress() << ' '; // 5
	fout_log << total_contact_stressXF.getCompressiveStress() << ' ' ; // 6
	fout_log << total_contact_stressXF.getLateralStress() << ' '; // 7
	fout_log << total_contact_stressXF.getNormalStress1() << ' '; // 8
	fout_log << total_contact_stressXF.getNormalStress2() << ' '; // 9
	fout_log << lz << ' '; // 10
	fout_log << lx << ' '; // 11
	fout_log << average_contact_number << ' '; //12
	fout_log << ave_force << ' '; //13
	fout_log << max_force << ' '; //14
	fout_log << r_min << ' ';// 15
	fout_log << r_max << ' ';// 16
	fout_log << sliding_disp_max << ' '; // 17
	fout_log << bending_angle_max << ' '; // 18
	fout_log << torsional_angle_max << ' '; // 19
	fout_log << max_move_step/max_velocity << ' '; // 20
	fout_log << counterRegenerate <<' '; //21
	fout_log << counterBreak << ' ';//22
	fout_log << n_bond-counterBreak << ' '; // 23
	fout_log << max_velocity << ' ';//24
	fout_log << max_ang_velocity << ' '; // 25
	fout_log << ave_bondforce << ' '; // 26
	fout_log << cnt_loop << ' '; // 27
	fout_log << max_displacement << ' '; // 28
	fout_log << prog_strain << endl; // 29
}

void System::outputConfiguration(char equilibrium)
{
	int number_of_live_bonds = n_bond-counterBreak;
	fout_particle << "# " << n_particle << ' ' << lx << ' ' << lz << endl;
	ForAllParticle {
		fout_particle << (*p_iter)->particle_number << ' ';
		fout_particle << (*p_iter)->p.x << ' ';
		fout_particle << (*p_iter)->p.z << ' ';
		fout_particle << (*p_iter)->orientation_angle  << endl;
	}
	fout_bond << "# num_bond " << number_of_live_bonds << endl;
	fout_bond << "# box_size " << lx << ' ' << lz << ' ' << endl;
	fout_bond << "# bond0_n " << bond0.fnc << ' ' << bond0.kn << endl;
	fout_bond << "# bond0_s " << bond0.fsc << ' ' << bond0.ks << endl;
	fout_bond << "# bond0_b " << bond0.mbc << ' ' << bond0.kb << endl;
	fout_bond << "# bond1_n " << bond1.fnc << ' ' << bond1.kn << endl;
	fout_bond << "# bond1_s " << bond1.fsc << ' ' << bond1.ks << endl;
	fout_bond << "# bond1_b " << bond1.mbc << ' ' << bond1.kb << endl;
	int cnt_bond = 0;
	ForAllBond {
		if ((*bond_iter)->status) {
			int i_p0, i_p1;
			(*bond_iter)->getParticleNumbers(i_p0, i_p1);
			fout_bond << i_p0 << ' ';
			fout_bond << i_p1 << ' ';
			double angle[2];
			angle[0] = (*bond_iter)->contactangle[0];
			angle[1] = (*bond_iter)->contactangle[1];
			while ( angle[0] < 0) {
				angle[0] += 2*M_PI;
			}
			while ( angle[0] > 2*M_PI) {
				angle[0] -= 2*M_PI;
			}
			while ( angle[1] < 0) {
				angle[1] += 2*M_PI;
			}
			while ( angle[1] > 2*M_PI) {
				angle[1] -= 2*M_PI;
			}
			fout_bond << angle[0] << ' ';
			fout_bond << angle[1] << ' ' ;
			fout_bond << (*bond_iter)->bondtype << endl;
			cnt_bond ++;
		}
	}
	if (cnt_bond != number_of_live_bonds) {
		exit(1);
	}
	fout_conf << "# " << equilibrium;
	fout_conf << ' ' << volume_fraction;
	fout_conf << ' ' << shear_strain;
	fout_conf << ' ' << total_contact_stressXF.getParticlePressure();
	fout_conf << ' ' << total_contact_stressXF.getShearStress();
	fout_conf << ' ' << lx;
	fout_conf << ' ' << lz;
	fout_conf << ' ' << rup_normal;
	fout_conf << ' ' << rup_shear;
	fout_conf << ' ' << rup_bend;
	fout_conf << ' ' << rup_torsion;
	fout_conf << ' ' << average_contact_number;
	fout_conf << endl;
	fout_conf << "P " << particle.size() << endl;
	ForAllParticle {
		fout_conf << (*p_iter)->p.x-lx_half << ' ';
		fout_conf << (*p_iter)->p.y-ly_half << ' ';
		fout_conf << (*p_iter)->p.z-lz_half << ' ';
		fout_conf << (*p_iter)->orientation.q[0] << ' ';
		fout_conf << (*p_iter)->orientation.q[1] << ' ';
		fout_conf << (*p_iter)->orientation.q[2] << ' ';
		fout_conf << (*p_iter)->orientation.q[3] << ' ';
		fout_conf << (*p_iter)->init_cluster << ' ' ;
		fout_conf << (*p_iter)->valCn_size() << endl;
	}
	fout_conf << "B " << number_of_live_bonds << endl;
	ForAllBond {
		if ((*bond_iter)->status) {
			int i_p0, i_p1;
			int initial_bond;
			if ((*bond_iter)->initial_bond) {
				initial_bond = 1;
			} else {
				initial_bond = 0;
			}
			(*bond_iter)->getParticleNumbers(i_p0, i_p1);
			fout_conf << i_p0 << ' ' << i_p1 << ' ';
			fout_conf << initial_bond << ' ' << (*bond_iter)->status << ' ';
			fout_conf << (*bond_iter)->val_F_norm() << ' ';
			fout_conf << (*bond_iter)->val_F_slid() << ' ';
			fout_conf << (*bond_iter)->val_M_bend() << ' ';
#ifndef TWODIMENSION
			fout_conf << (*bond_iter)->val_M_tors() << ' ';
#else
			fout_conf << 0 << ' ';
#endif
			fout_conf << (*bond_iter)->cnt_regeneration << endl;
		}
	}
	lz_last_output = lz;
}

void System::monitorDeformation(char equilibrium)
{
	fout_deform << "# " << equilibrium;
	fout_deform << ' ' << total_contact_stressXF.getParticlePressure();
	fout_deform << ' ' << total_contact_stressXF.getShearStress();
	fout_deform << ' ' << volume_fraction;
	fout_deform << ' ' << shear_strain;
	fout_deform << ' ' << lz;
	fout_deform << ' ' << rup_normal;
	fout_deform << ' ' << rup_shear;
	fout_deform << ' ' << rup_bend;
	fout_deform << ' ' << rup_torsion;
	fout_deform << ' ' << average_contact_number;
	fout_deform << endl;
	int count_bond = 0;
	for (int i=0; i < n_bond; i++) {
		if (bond[i]->initial_bond) {
			count_bond ++;
		}
	}
	fout_deform << "B " << count_bond   <<endl;
	for (int i=0; i < n_bond; i++) {
		if (bond[i]->initial_bond) {
			bond[i]->monitor_state(fout_deform);
		}
	}
}

void System::calcVolumeFraction()
{
#ifndef TWODIMENSION
	double volume = lx*ly*lz;
	static const double volume1 = (4.0/3.0)*M_PI;
	volume_fraction = n_particle*volume1 / volume;
#else
	cerr << "lz = " << lz << endl;
	volume_fraction = n_particle* M_PI/(lx * lz);  // Area fraction
#endif
}

/*
 *	Check the state
 *	- average number of bond par particle
 *  -
 *  - The total force should be zero in the equilibrium
 *  - Because the force from wall is not added here,
 *  - The force is calculated for only active particles (not contact to the walls.)
 */
void System::checkState()
{
	int sum_num_of_contacts_for_particle = 0;
	for (int i=0; i < n_particle; i++) {
		sum_num_of_contacts_for_particle += particle[i]->valCn_size();
	}
	average_contact_number = ((double)sum_num_of_contacts_for_particle)/n_particle;
	for (int i=0; i < bond.size(); i++) {
		bond[i]->addContactForce();
	}
	double bforce_max = 0;
	double bforce_sum = 0;
	for (int i=0; i < bond.size(); i++) {
		double bforce = bond[i]->forceNorm();
		bforce_sum += bforce;
		if (bforce > bforce_max) {
			bforce_max = bforce;
		}
	}
	double sum_force = 0;
	max_force = 0;
	max_velocity = 0;
	max_ang_velocity = 0;
	for (int i=0; i < n_particle; i++) {
		double force = particle[i]->valForce();
		double velocity = particle[i]->valVelocity();
		double ang_velocity = particle[i]->valOmega();
		sum_force += force;
		if (force > max_force) {
			max_force = force ;
		}
		if (velocity > max_velocity) {
			max_velocity = velocity;
		}
		if (ang_velocity > max_ang_velocity) {
			max_ang_velocity = ang_velocity;
		}
	}
	ave_force = sum_force/n_particle;
	for (int i=0; i < n_particle; i++) {
		particle[i]->resetForce();
	}
	calcLocalStrains();
}

/*
 *	Adjust
 *	- center of mass is set to the center
 *  - Bond vector
 */
void System::simuAdjustment()
{
	ForAllParticle {
		(*p_iter)->setNorm_u();
	}
}

void System::makeNeighborPB()
{
	grid->remake(particle);
	ForAllParticle{
		(*p_iter)->makeNeighbor();
	}
}

void System::checkBondFailure()
{
	foreach(vector<Bond *>, bond, it_bond) {
		(*it_bond)->cheackBondStress();
	}
}

void System::regeneration()
{
	foreach (vector<int>, regeneration_bond, b) {
		bond[*b]->regeneration();
		counterRegenerate ++;
	}
	regeneration_bond.clear();
}

void System::rupture()
{
	foreach (vector<int>, rupture_bond, b) {
		bond[*b]->rupture();
		counterBreak ++;
	}
	rupture_bond.clear();
}

void System::makeInitialBond(double generation_distance)
{
	double tmp = sq_dist_generate;
	sq_dist_generate = sq(generation_distance);
	makeNeighborPB();
	generateBond();
	sq_dist_generate = tmp;
}

double System::system_volume()
{
#ifdef TWODIMENSION
	return lx*lz*2;
#else
	return lx*ly*lz;
#endif
}
