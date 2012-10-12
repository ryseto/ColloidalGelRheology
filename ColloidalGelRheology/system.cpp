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
	version="0";
	ct = new ConTable;
	grid = new Grid;
}

System::~System(){
	fout_conf.close();
	fout_data.close();
	fout_log.close();
	fout_deform.close();
	delete ct;
	delete grid;
}

void System::TimeDevStrainControlCompactionEuler(){
	if(prog_strain){
		wl[0]->compactionStrainControl( wall_velocity ); // bot
		wl[1]->compactionStrainControl(-wall_velocity ); // top
	}
 	foreach( vector<Bond *>, bond_active, it_bond){
		(*it_bond)->addContactForce();
	}
	//////////////////////////////
	// calculate
	// x*_{j+1} and v*_{j+1}
	foreach( vector<Particle *>, particle_active, it_particle){
		(*it_particle)->move_Euler();
	}
	time += dt;
}

void System::TimeDevStrainControlShearEuler(){
	if(prog_strain){
		wl[0]->shearingStrainControl(-wall_velocity ); // bot
		wl[1]->shearingStrainControl(+wall_velocity ); // top
	}
 	foreach( vector<Bond *>, bond_active, it_bond){
		(*it_bond)->addContactForce();
	}
	//////////////////////////////
	// calculate
	// x*_{j+1} and v*_{j+1}
	foreach( vector<Particle *>, particle_active, it_particle){
		(*it_particle)->move_Euler();
	}
	time += dt;
}

void System::TimeDevBendingTest(){

	particle_active[0]->stackForce(-f_ex, 0);
	particle_active[5]->stackForce(2*f_ex, 0);
	particle_active[10]->stackForce(-f_ex, 0);

 	foreach( vector<Bond *>, bond_active, it_bond){
		(*it_bond)->addContactForce();
	}
	
	//////////////////////////////
	// calculate
	// x*_{j+1} and v*_{j+1}
	foreach( vector<Particle *>, particle_active, it_particle){
		(*it_particle)->move_Euler();
	}
	time += dt;
}

bool System::mechanicalEquilibrium(){
	bool mechanical_equilibrium = false;
	if (simulation == 's'){
		/* Check the equilibrium condition for shear
		 * (1.1) The difference between the stresses at the two walls is enough small.
		 * (1.2) The stress is less than the minimum stress.
		 * (2) The change of the measured stress is enough small.
		 * (3) The maximum velocity of the particles are enough small.
		 * (4) The maximum angular velocity of the particles are enough small.
		 *  ---------> Bond breakup requires (1)-(4).
		 * (5) No restructuring.
		 *  ---------> The equilibrium requires (1)-(5).
		 */
		if ((diff_stress_x < diff_stress_convergence
			 ||  stress_x < stress_minimum)
			&& stress_x_change < stress_change_convergence
			&& max_velocity < max_velocity_convergence
			&& max_ang_velocity < max_ang_velocity_convergence
			&& counterRegenerate_before == counterRegenerate
			){
			mechanical_equilibrium =  true;
		}
	} else {
		/* Check the equilibrium condition
		 * (1) The change of the measured stress is enough small.
		 * (2) The difference between the stresses at the two walls is enough small.
		 * (3) The maximum velocity of the particles are enough small.
		 * (4) The maximum angular velocity of the particles are enough small.
		 *  ---------> Bond breakup requires (1)-(4).
		 * (5) No restructuring.
		 *  ---------> The equilibrium requires (1)-(5).
		 */

		if (stress_z_change < stress_change_convergence
			&& diff_stress_z < diff_stress_convergence
			&& max_velocity < max_velocity_convergence
			&& max_ang_velocity < max_ang_velocity_convergence
			&& counterRegenerate_before == counterRegenerate
			){
			mechanical_equilibrium =  true;
		}
	}
	counterRegenerate_before = counterRegenerate;
	return mechanical_equilibrium;
}

void System::preProcesses(){
	initialprocess = true;
	setSimulationViscosity();
	preparationOutput();
	setBondGenerationDistance(2.0);
	setWall();
	initGrid();
	//outputParameter(sy, fout_yap);
	initDEM();
	makeInitialBond(2.0001);
	initialprocess = false;
	checkActiveElements();
	calcVolumeFraction();
	outputConfiguration('n');
	while(checkPercolation() == false){
		cerr << "not percolation" << endl;
		shiftForPercolation();
		outputConfiguration('n');
		makeNeighbor();
		generateBond();
		wl[1]->addNewContact(particle_active);
		wl[0]->addNewContact(particle_active);

	}
	calcVolumeFraction();
	dt = dt_max;
	setFirstTarget();
	checkState();
	outputData();
	outputConfiguration('e');
}

void System::timeEvolution(){
	double t_next = time + interval_convergence_check*dt_max;
	int calc_count = 0;
	counter_relax_for_restructuring = 0;
	while ( time < t_next ){
		if (reachStrainTarget() ||
			counter_relax_for_restructuring < relax_for_restructuring){
			prog_strain = false;
		} else {
			prog_strain = true;
		}
		/* Main time evolution */
		if (simulation == 'c') {
			TimeDevStrainControlCompactionEuler();
		} else {
			TimeDevStrainControlShearEuler();
		}
		/* Breakup of bonds */
		if (counter_relax_for_restructuring++ >= relax_for_restructuring){
			checkBondFailure();
			if (!regeneration_bond.empty()){
				//outputRestructuring();
				regeneration_onebyone();
				counter_relax_for_restructuring = 0;
			}
			if (!rupture_bond.empty()){
				//outputRestructuring();
				rupture();
				counter_relax_for_restructuring = 0;
			}
		}
		/* Bond generation */
		if ( calc_count++ % interval_makeNeighbor == 0 ){
			makeNeighbor();
		}
		generateBond();
		wl[0]->addNewContact(particle_active);
		wl[1]->addNewContact(particle_active);
	}
	if (simulation == 'c'){
		calcVolumeFraction();
	}
}

void System::middleProcedures(){
	simuAdjustment();
	calcStress();

	checkState();
	checkStressDiff();
	makeNeighbor();
	optimalTimeStep();
	if ( checkOutputStrain(1.0) ){
		outputConfiguration('n');
	}
	output_log();
}

void System::strainControlSimulation(){
	preProcesses();
	while(true){
		timeEvolution();
		middleProcedures();
		outlog();
		if ( reachStrainTarget() ){
			if (mechanicalEquilibrium()){
				cerr <<"Equilibrium!"<< endl;
				calcStress();
				/*  Output */
				outputData();
				outputConfiguration('e');
				monitorDeformation('e');
				/*  Prepare next step */
				setTarget();
				if (checkEndCondition()){
					break;
				}
			}
		}
	}
	return;
}

void System::bendingSimulation(){
	initialprocess = true;
	setSimulationViscosity();
	preparationOutput();
	setBondGenerationDistance(2.0);
	
	initGrid();
	//outputParameter(sy, fout_yap);
	initDEM();
	makeInitialBond(2.0002);
	initialprocess = false;
	checkActiveElements();
	
	cerr << "particle_active " << particle_active.size() << endl;
	cerr << "bond_active " << bond_active.size() << endl;
	dt = dt_max;
	checkState();
	int cnt  = 0 ;
	double diff_x_previous = 0;
	double time_interval = 300;
	double time_max = time_interval;
	for (int i = 0; i < 2 ; i++){
		if (i==1){
			f_ex.set(0,0,0);
		}
		
		while(true){
			TimeDevBendingTest();
			checkBondFailure();
			if (!regeneration_bond.empty()){
				//outputRestructuring();
				regeneration_onebyone();
				counter_relax_for_restructuring = 0;
			}
			if (!rupture_bond.empty()){
				//outputRestructuring();
				rupture();
				counter_relax_for_restructuring = 0;
			}
			if (cnt++ % 100 ==0){
				
				fout_conf << "# " << f_ex.x << endl;
				foreach( vector<Particle *>, particle_active, it_particle){
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
				double diff_x = particle_active[5]->p.x - particle_active[0]->p.x ;
				
				fout_data << time << ' ' << diff_x << ' ' << f_ex.x << ' ' << counterRegenerate << endl;
				diff_x_previous = diff_x;
			}
			if ( time > time_max ){
				break;
			}
		}
		time_max +=time_interval;
	}
	return;
}


void System::outlog(){
	if (simulation=='s'){
		cerr << strain_x << " --> " << strain_target  << endl;
	} else if (simulation=='c'){
		cerr << volume_fraction << " --> " << vf_target  << endl;
	}
}

bool System::reachStrainTarget(){
	if (simulation=='s'){
		// shear
		if ( strain_x >= strain_target){
			return true;
		} else {
			return false;
		}
	} else {
		// compaction
		if ( volume_fraction >= vf_target){
			return true;
		} else {
			return false;
		}
	}
}

void System::setTarget(){
	if (simulation=='s'){
		strain_target += step_strain_x;
	} else {
		vf_target *= volumefraction_increment;
	}
}

void System::setFirstTarget(){
	if (simulation == 's'){
		strain_target = step_strain_x;
	} else {
		vf_target = volume_fraction / initial_compaction;
	}
}

bool System::checkEndCondition(){
	bool end_simulation = false;
	if (simulation =='s' ){
		if (strain_x > max_strain_x){
			end_simulation = true;
		}
	} else if (simulation =='c' ){
		if (volume_fraction > max_volume_fraction){
			end_simulation = true;
		}
	}
	return end_simulation;
}

void System::checkStressDiff(){
	if (simulation == 's'){
		if (diff_stress_x == 0){
			stress_x_change = 0;
			stress_z_change = 0;
		} else {
			stress_x_change = abs(stress_x - stress_x_before)/stress_x;
			stress_z_change = abs(stress_z - stress_z_before)/stress_z;
		}
		strain_x = abs( wl[0]->x - wl[1]->x )/lz;
		strain_z = (lz_init - lz)/lz_init;
		
	} else {
		stress_z_change = abs(stress_z - stress_z_before)/stress_z;
		stress_x_change = abs(stress_x - stress_x_before)/stress_x;
		strain_z = (lz_init - lz)/lz_init;
	}
	stress_z_before = stress_z;
	stress_x_before = stress_x;
}

bool System::checkOutputStrain(double diff){
	if (simulation == 's'){
		if (strain_x > strain_x_last_output + diff){
			lz_last_output = lz;
			return true;
		}
	} else {
		if ( lz < lz_last_output - diff ){
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

void System::initDEM(){
	ct->set(n_particle);
	rup_normal = rup_shear = rup_bend = rup_torsion = 0;
	rupture_bond.clear();
	regeneration_bond.clear();
	time = 0.;
	counterBreak = 0;
	counterRegenerate = 0;
	stress_z = 0.;
	stress_x = 0.;
	strain_x = 0.;
	strain_z = 0.;
	stress_z_before = 0;
	stress_x_before = 0;
	counter_relax_for_restructuring = 0;
	lz_last_output = lz_init;
	counterRegenerate_before = 0;
}

void System::preparationOutput(){
	char fn_snapshot[128];
	char fn_data[128];
	char fn_log[128];
	char fn_conf[128];
	char fn_deform[128];
	
	char name_string[128];
	if(simulation == 'c'|| simulation == 's'){
		sprintf(name_string, "%s_%s_%s",
				parameters, init_cluster, version.c_str());
	} else if(simulation == 'b'){
		sprintf(name_string, "%s_%1.2f", parameters, f_ex.x );
	}

	string fn_common = (string)name_string ;
	cerr << version << endl;
	cerr << fn_common.c_str() << endl;

	sprintf(fn_data, "data_%s.dat", fn_common.c_str());
	sprintf(fn_log, "log_%s.dat", fn_common.c_str());
	sprintf(fn_conf, "conf_%s.dat", fn_common.c_str());
	sprintf(fn_deform, "deform_%s.dat", fn_common.c_str());
	
	if ( simulation == 'b'){
		sprintf(fn_snapshot, "ss_%s.yap", fn_common.c_str());
		fout_yap.open(fn_snapshot);
	}
	fout_data.open(fn_data);
	fout_log.open(fn_log);
	fout_conf.open(fn_conf);
	fout_deform.open(fn_deform);
}

void System::setInitClusterFile(string init_cluster_file_)
{
	sprintf(init_cluster_file, "%s", init_cluster_file_.c_str());
	string s_init_cluster_file = init_cluster_file;
	unsigned long i_backslash = s_init_cluster_file.find_last_of( "/") + 1;
	unsigned long i_extention = s_init_cluster_file.find( ".dat" );
	sprintf(init_cluster, "%s",
			(s_init_cluster_file.substr( i_backslash, i_extention-i_backslash)).c_str());
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

void System::generateBond(){
	int bond_number_before = n_bond;
	for (int i=0; i < n_particle; i++){
		particle[i]->generateBond();
	}
	if (bond_number_before < n_bond){
		for (int i = bond_number_before; i < n_bond; i++){
			bond_active.push_back( bond[i] );
			bond[i]->number_activebond  = bond_active.size() - 1;
		}
	}
	return;
}

void System::generateBondAll(){
	for (int i=0; i < n_particle; i++){
		particle[i]->generateBond();
	}
}

void System::setSimulationViscosity()
{
	eta = eta_factor*bond0.c_norm;
	eta_rot = (4.0/3.0)*eta;
}

void System::optimalTimeStep(){
	double dt_new1 = max_move_step / max_velocity;
	double dt_new2 = max_move_step / max_ang_velocity;
	double dt_new;
	if (dt_new1 < dt_new2){
		dt_new = dt_new1;
	} else {
		dt_new = dt_new2;
	}
	if ( dt_new > dt_max ){
		dt_new = dt_max;
	} else {
		cerr << "small time step"  << dt_new << endl;
	}
	dt = dt_new;
}

void System::calcStress(){
	double stress_z_bot;// = wl[0]->stressSensor_z();
	double stress_z_top;// = wl[1]->stressSensor_z();
	double stress_x_bot;// = wl[0]->stressSensor_x();
	double stress_x_top;// = wl[1]->stressSensor_x();
	wl[0]->stressSensor(stress_x_bot, stress_z_bot);
	wl[1]->stressSensor(stress_x_top, stress_z_top);
	stress_z = 0.5*(abs(stress_z_top-stress_z_bot));
	diff_stress_z = abs( stress_z_bot + stress_z_top )/stress_z;
	stress_x = 0.5*(stress_x_top - stress_x_bot);
	diff_stress_x = abs(stress_x_bot + stress_x_top )/stress_x;
//	if (stress_x < 1e-8){
//		diff_stress_x = 0;
//	}else{
//		diff_stress_x = abs(stress_x_bot + stress_x_top )/stress_x;
//	}
	
}


void prepareBond(BondParameter & _bond){
	double c = 1;
	_bond.n_max = _bond.fnc/_bond.kn;
	_bond.s_max = _bond.fsc/_bond.ks;
	_bond.b_max = _bond.mbc/_bond.kb;
	_bond.t_max = _bond.mtc/_bond.kt;
	_bond.c_norm = c*2*sqrt( _bond.kn);
	_bond.c_slid = c*2*sqrt( _bond.ks);
	_bond.c_bend = c*2*sqrt( (2./5.)*_bond.kb);
	_bond.c_tort = c*2*sqrt( (2./5.)*_bond.kt);
	
	cerr << "k:" << _bond.kn << ' ' << _bond.ks << ' ' << _bond.kb << ' ' << _bond.kt << endl;
	cerr << "critical_strain:" << _bond.n_max << ' ' << _bond.s_max << ' ';
	cerr << _bond.b_max << ' ' << _bond.t_max << endl;
	cerr << "c:" << _bond.c_norm << ' ' << _bond.c_slid << ' ' << _bond.c_bend << ' ' << _bond.c_tort << endl;
}

void readBond(const string &codeword, const string &value, BondParameter &bondparameter)
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
	keylist["reinforce_factor:"]=10; const int _reinforce_factor = 10;
	cerr << codeword << ' ' << value << endl;
	switch(keylist[codeword]){
		case _fnc: bondparameter.fnc = atof(value.c_str()); break;
		case _fsc: bondparameter.fsc = atof(value.c_str()); break;
		case _mbc: bondparameter.mbc = atof(value.c_str()); break;
		case _mtc: bondparameter.mtc = atof(value.c_str()); break;
		case _kn: bondparameter.kn = atof(value.c_str()); break;
		case _ks: bondparameter.ks = atof(value.c_str()); break;
		case _kb: bondparameter.kb = atof(value.c_str()); break;
		case _kt: bondparameter.kt = atof(value.c_str()); break;
		case _kn3: bondparameter.kn3 = atof(value.c_str()); break;
		case _reinforce_factor: bondparameter.reinforce_factor = atof(value.c_str()); break;
		default:
			cerr << "Bond: The codeword " << codeword << " is'nt associated with an parameter" << endl;
			exit(1);
	}
	
}


void setBondParameter( BondParameter &bondparameter, string &bond_file){
	string path_bond = bond_file;
	ifstream fin;
	if(fin.is_open())
		fin.close();
	fin.open(path_bond.c_str());
	string codeword;
	string value;
	while (!fin.eof()){
		fin >> codeword;
		if ( codeword == "#") {
			char buf[1024]; fin.get(buf, 1024);
		} else if (codeword == "!"){
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

void System::readBondParameter(){
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
	if (fin_parameter.is_open()){
		cerr << parameters_file << " is opend." << endl;
	}
	string codeword;
	string value;
	while (!fin_parameter.eof()){
		fin_parameter >> codeword;
		if ( codeword == "#") {
			char buf[1024]; fin_parameter.get(buf, 1024);
		} else if (codeword == "!"){
			break;
		} else {
			fin_parameter >> value ;
			if (fin_parameter.eof()) break;
			//set_value(codeword, value, unit);
			readParameter(codeword, value);
		}
	}
	fin_parameter.close();
	if (!fin_parameter.is_open()){
		cerr << parameters_file << " is closed." << endl;
	}
}

void System::readParameter(const string &codeword, const string &value)
{
	map<string,int> keylist;
	keylist["bond0_file:"]=2; const int _bond0_file = 2;
	keylist["bond1_file:"]=3; const int _bond1_file = 3;
	keylist["initial_compaction:"]=11; const int _initial_compaction=11;
	keylist["volumefraction_increment:"]=12; const int _volumefraction_increment=12;
	keylist["max_volume_fraction:"]=13; const int _max_volume_fraction=13;
	
	keylist["step_strain_x:"]=14; const int _step_strain_x=14;
	keylist["max_strain_x:"]=15; const int _max_strain_x=15;
	keylist["wall_velocity:"]=20;const int _wall_velocity=20;
	keylist["dt_max:"]=21; const int _dt_max=21;
	keylist["max_move_step:"]=22; const int _max_move_step=22;
	keylist["eta_factor:"]=23;const int _eta_factor=23;
	keylist["relax_for_restructuring:"]=30; const int _relax_for_restructuring=30;
	keylist["interval_convergence_check:"]=31; const int _interval_convergence_check=31;
	keylist["interval_makeNeighbor:"]=32; const int _interval_makeNeighbor=32;
	keylist["max_velocity_convergence:"]=33; const int _max_velocity_convergence=33;
	keylist["max_ang_velocity_convergence:"]=34; const int _max_ang_velocity_convergence=34;
	keylist["diff_stress_convergence:"]=35;const int _diff_stress_convergence=35;
	keylist["stress_change_convergence:"]=36;const int _stress_change_convergence=36;
	
	keylist["stress_minimum:"]=37; const int _stress_minimum=37;
	cerr << codeword << ' ' << value << endl;
	switch(keylist[codeword]){
		case _bond0_file: bond0_file = value; break;
		case _bond1_file: bond1_file = value; break;
		case _initial_compaction: initial_compaction = atof(value.c_str()) ; break;
		case _interval_convergence_check: interval_convergence_check = atoi(value.c_str()); break;
		case _interval_makeNeighbor: interval_makeNeighbor = atoi(value.c_str()); break;
		case _max_velocity_convergence: max_velocity_convergence = atof(value.c_str()); break;
		case _wall_velocity: wall_velocity = atof(value.c_str()); break;
		case _dt_max: dt_max = atof(value.c_str()); break;
		case _eta_factor: eta_factor = atof(value.c_str()); break;
		case _max_ang_velocity_convergence: max_ang_velocity_convergence = atof(value.c_str()); break;
		case _max_volume_fraction: max_volume_fraction = atof(value.c_str()); break;
		case _step_strain_x: step_strain_x = atof(value.c_str()); break;
		case _max_strain_x: max_strain_x = atof(value.c_str()); break;
		case _max_move_step: max_move_step = atof(value.c_str()); break;

			
		case _volumefraction_increment: volumefraction_increment = atof(value.c_str());; break;
		case _diff_stress_convergence: diff_stress_convergence = atof(value.c_str());; break;
		case _stress_change_convergence: stress_change_convergence = atof(value.c_str());; break;
		case _relax_for_restructuring: relax_for_restructuring = atoi(value.c_str()); break;
		case _stress_minimum: stress_minimum = atof(value.c_str()); break;
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

	fin.open( path.c_str());
	if (! fin.is_open() ){
		cerr << "no initial file" << endl;
		exit(1);
	}
	
	lx = atof(width_str);
	lz_init = atof(height_str);

#ifdef TWODIMENSION
	ly = 0;
	ly0 = 0;
#else
	ly = lx;
	ly0 = ly/2;
#endif
	lx0 = lx/2;
	lz0 = lz_init/2;
	
	double x, y, z;
	int i_cluster;
	if (simulation == 'c' || simulation == 's'){
		do{
#ifdef TWODIMENSION
			fin >> x >> z >> i_cluster ;
#else
			fin >> x >> y >> z;
#endif
			cerr << x << ' '<< z << endl;
#ifdef TWODIMENSION
			vec3d new_p(x, y, z);
#else
			vec3d new_p(x, 0, z);
#endif
			init_aggregate.push_back(new_p);
			init_aggregate_cluster.push_back(i_cluster);
			
		}while (!fin.eof());
	}
	init_aggregate.pop_back();
	int count_particle_number = 0;
	cerr << "N=" << init_aggregate.size() << endl;
	for (int i = 0 ; i < init_aggregate.size(); i++){
		vec3d new_p = init_aggregate[i];
		particle.push_back(new Particle(count_particle_number, new_p,
										init_aggregate_cluster[i],
										*this) );
		++ count_particle_number;
	}
	n_particle = count_particle_number;
	if (n_particle == 0){
		cerr << "failed to import the initial configuration.\n";
		exit(1);
	}
}

void System::setLinearChain(int number_of_particles){
	int m = number_of_particles;
	int count_particle_number = 0;
	for (int i = 0 ; i < m; i++){
		vec3d new_p(0.5*lx, 0, 2*i - number_of_particles + 1 + 0.5*lz);
		particle.push_back(new Particle(count_particle_number, new_p,
										0, *this) );
		++ count_particle_number;
	}
	n_particle = count_particle_number;
	if (n_particle == 0){
		cerr << "failed to import the initial configuration.\n";
		exit(1);
	}

	
	
}


void System::shiftCenterOfMass(vector<vec3d> &p)
{
	vec3d cm(0,0,0);
	foreach( vector<vec3d>, p, p_iter){
		cm += (*p_iter);
	}
	cm *= 1.0/p.size();
	foreach( vector<vec3d>, p, p_iter){
		*p_iter -= cm;
		p_iter->x += lx/2;
#ifdef TWODIMENSION
		p_iter->y += ly/2;
#endif
		p_iter->z += lz/2;
	}
}

/* Checking active bonds at the first step.
 *
 *
 */
void System::checkActiveElements(){
	for (int i=0; i < n_bond; i++){
		if (bond[i]->status){
			bond_active.push_back( bond[i] );
			bond[i]->number_activebond  = bond_active.size() - 1;
		}
	}
	for (int i=0; i<n_particle; i++){
		if ( particle[i]->wall == false){
			particle_active.push_back( particle[i]);
		}
	}
}

void System::calcLocalStrains(){
	/* The minimum and maximum distances between two particles.
	 */
	r_min = 2;
	r_max = 2;
	bending_angle_max = 0;
	torsional_angle_max = 0;
	sliding_disp_max = 0;
	for (int i=0; i < n_bond; i++){
		if (bond[i]->status){
			if (bond[i]->r < r_min){
				r_min = bond[i]->r;
			}
			if (bond[i]->r > r_max){
				r_max = bond[i]->r;
			}
			if (bond[i]->valSlidingDisplacement() > sliding_disp_max ){
				sliding_disp_max = bond[i]->valSlidingDisplacement();
			}
			if (bond[i]->valBendingAngle() > bending_angle_max){
				bending_angle_max = bond[i]->valBendingAngle();
			}
#ifndef TWODIMENSION
			if (bond[i]->valTorsionalAngle() > torsional_angle_max){
				torsional_angle_max = bond[i]->valTorsionalAngle();
			}
#endif
		}
	}
	return;
}

void System::initGrid(){
	if (n_particle == 0)
	{
		cerr << "before setting system box, ";
		cerr << "it is required to import configuration of particles.";
		cerr << endl;
		exit(1);
	}
	double cell_h = 2.2;
	grid->init(n_particle, lx, ly, lz_init, cell_h);
}

//void System::setWallStress(double value_x, double value_y,  double value_z){
//	if (value_x == -1){
//		wl[0]->x_movable = false;
//		wl[1]->x_movable = false;
//		stress_x = 0;
//	} else {
//		wl[0]->x_movable = true;
//		wl[1]->x_movable = true;
//		stress_x = value_x;
//	}
//	if (value_y == -1){
//		wl[0]->y_movable = false;
//		wl[1]->y_movable = false;
//		stress_y = 0;
//	} else {
//		wl[0]->y_movable = true;
//		wl[1]->y_movable = true;
//		stress_y = value_y;
//	}
//	
//	if (value_z == -1){
//		wl[0]->z_movable = false;
//		wl[1]->z_movable = false;
//		stress_z = 0;
//	} else {
//		wl[0]->z_movable = true;
//		wl[1]->z_movable = true;
//		stress_z = value_z;
//	}
//
//#ifndef TWODIMENSION
//	double fx_wall = stress_x*lx*ly;
//	double fy_wall = stress_y*lx*ly;
//	double fz_wall = stress_z*lx*ly;
//#else
//	double fx_wall = stress_x*lx*2;
//	double fy_wall = stress_y*lx*2;
//	double fz_wall = stress_z*lx*2;
//#endif
//	
//	force_wall.set(-fx_wall, -fy_wall, fz_wall);
//	cerr << "stress = " << stress_z << ',' << stress_x << " Pa ";
//	cerr << " Fwall=" << fz_wall << ',' << fx_wall << " N" << endl;
//	
//}

void System::renew_Lz(){
	lz = wl[1]->z - wl[0]->z;
}

void System::reset_strain_x(){
	wl[1]->x = lx0;
	wl[0]->x = lx0;
	strain_x = 0;
}

void System::outputData(){
	static bool firsttime = true;
	static int rup_normal_previous = 0;
	static int rup_shear_previous = 0;
	static int rup_bend_previous = 0;
	static int rup_torsion_previous = 0;
	static double strain_previous = 0;
	double del_strain = 0;
	
	if ( firsttime ){
		firsttime = false;
		fout_data << "#1 volume_fraction\n";
		fout_data << "#2 stress_z\n";
		fout_data << "#3 stress_x\n";
		fout_data << "#4 strain_z\n";
		fout_data << "#5 strain_x\n";
		fout_data << "######  rupture: total number \n";
		fout_data << "#6  rup_normal(total number)\n";
		fout_data << "#7  rup_shear(total number)\n";
		fout_data << "#8  rup_bend(total number)\n";
		fout_data << "#9 rup_torsion(total number)\n";
		fout_data << "######  rupture: rate\n";
		fout_data << "#10 rate_rup_normal\n";
		fout_data << "#11 rate_rup_shear\n";
		fout_data << "#12 rate_rup_bend\n";
		fout_data << "#13 rate_rup_torsion\n";
		fout_data << "###### \n";
		fout_data << "#14 average_contact_number\n";
		fout_data << "#15 del_strain\n";
	}
	
	if (simulation == 'c'){
		del_strain = (strain_z - strain_previous);
		strain_previous = strain_z;
	} else if (simulation == 's'){
		del_strain = (strain_x - strain_previous);
		strain_previous = strain_x;
	}
	double rate_rup_normal, rate_rup_shear, rate_rup_bend, rate_rup_torsion;
	if (del_strain == 0){
		rate_rup_normal = 0;
		rate_rup_shear = 0;
		rate_rup_bend = 0;
		rate_rup_torsion = 0;
	} else {
		rate_rup_normal = (rup_normal - rup_normal_previous)/del_strain;
		rate_rup_shear = (rup_shear - rup_shear_previous)/del_strain;
		rate_rup_bend = (rup_bend - rup_bend_previous)/del_strain;
		rate_rup_torsion = (rup_torsion - rup_torsion_previous)/del_strain ;
	}
	fout_data << volume_fraction << ' ' ; //1
	fout_data << stress_z << ' ' ;  // 2
	fout_data << stress_x << ' ' ; // 3
	fout_data << strain_z << ' ' ; // 4
	fout_data << strain_x << ' ' ; // 5
	///////////////////////////////////////
	fout_data << rup_normal << ' '; //6
	fout_data << rup_shear << ' ';  //7
	fout_data << rup_bend << ' ';   //8
	fout_data << rup_torsion << ' ';//9
	///////////////////////////////////////
	fout_data << rate_rup_normal << ' '; //10
	fout_data << rate_rup_shear << ' ';  //11
	fout_data << rate_rup_bend << ' ';   //12
	fout_data << rate_rup_torsion << ' ';//13
	///////////////////////////////////////
	fout_data << average_contact_number << ' '; //14
	fout_data << del_strain << ' '; // 15
	fout_data << endl;
	rup_normal_previous = rup_normal;
	rup_shear_previous = rup_shear;
	rup_bend_previous = rup_bend;
	rup_torsion_previous = rup_torsion;
}

void System::output_log()
{
	static bool firsttime = true;
	if ( firsttime ){
		firsttime = false;
		fout_log << "#1 time\n";
		fout_log << "#2 lz\n";
		fout_log << "#3 volume_fraction\n";
		fout_log << "#4 stress_z\n";
		fout_log << "#5 stress_x\n";
		fout_log << "#6 strain_z\n";
		fout_log << "#7 strain_x\n";
		fout_log << "#8 average_contact_number\n";
		fout_log << "#9 ave_force\n";
		fout_log << "#10 max_force\n";
		fout_log << "#11 r_min\n";
		fout_log << "#12 r_max\n";
		fout_log << "#13 sliding_disp_max\n";
		fout_log << "#14 bending_angle_max\n";
		fout_log << "#15 torsional_angle_max\n";
		fout_log << "#16 dt\n";
		fout_log << "#17 counterRegenerate\n";
		fout_log << "#18 counterBreak\n";
		fout_log << "#19 bond number\n";
		fout_log << "#20 max_velocity\n";
		fout_log << "#21 max_ang_velocity\n";
		fout_log << "#22 ave_bondforce\n";
		fout_log << "#23 diff_stress_z\n";
		fout_log << "#24 diff_stress_x\n";
		fout_log << "#25 stress_z_change\n";
		fout_log << "#26 stress_x_change\n";
		fout_log << "#27 kinetic_energy\n";
	}
	
	fout_log << time << ' '; //1
	fout_log << lz << ' '; //2
	fout_log << volume_fraction << ' ';// 3
	fout_log << stress_z << ' ' ; // 4
	fout_log << stress_x << ' ' ; // 5
	fout_log << strain_z << ' ' ; // 6
	fout_log << strain_x << ' ' ; // 7
	fout_log << average_contact_number << ' '; //8
	fout_log << ave_force << ' '; //9
	fout_log << max_force << ' '; //10
	fout_log << r_min << ' ';// 11
	fout_log << r_max << ' ';// 12
	fout_log << sliding_disp_max << ' '; // 13
	fout_log << bending_angle_max << ' '; // 14
	fout_log << torsional_angle_max << ' '; // 15
	fout_log << max_move_step / max_velocity << ' '; // 16
	fout_log << counterRegenerate <<' '; //17
	fout_log << counterBreak << ' ';//18
	fout_log << n_bond - counterBreak << ' '; // 19
	fout_log << max_velocity << ' ';//20
	fout_log << max_ang_velocity << ' '; // 21
	fout_log << ave_bondforce << ' '; // 22
	fout_log << diff_stress_z << ' ' ; //23
	fout_log << diff_stress_x << ' ' ; //24
	fout_log << stress_z_change << ' '; // 25
	fout_log << stress_x_change << ' '; // 26
	fout_log << kinetic_energy  << endl; // 27
}

void System::outputConfiguration(char equilibrium){
	lz0 = (wl[1]->z + wl[0]->z) / 2;
	fout_conf << "# " << equilibrium;
	fout_conf << ' ' << volume_fraction;
	fout_conf << ' ' << stress_z;
	fout_conf << ' ' << stress_x;
	fout_conf << ' ' << strain_z;
	fout_conf << ' ' << strain_x;
	if (simulation == 'c' || simulation == 's'){
		fout_conf << ' ' << wl[0]->z - lz0;
		fout_conf << ' ' << wl[1]->z - lz0;
	} else {
		fout_conf << ' ' << 0 ;
		fout_conf << ' ' << 0 ;
	}
	fout_conf << ' ' << rup_normal;
	fout_conf << ' ' << rup_shear;
	fout_conf << ' ' << rup_bend;
	fout_conf << ' ' << rup_torsion;
	fout_conf << ' ' << average_contact_number;
	fout_conf << endl;
	fout_conf << "P " << particle.size() << endl;
	ForAllParticle{
		fout_conf << (*p_iter)->p.x - lx0  << ' ';
		fout_conf << (*p_iter)->p.y - ly0  << ' ';
		fout_conf << (*p_iter)->p.z - lz0  << ' ';
		fout_conf << (*p_iter)->orientation.q[0] << ' ';
		fout_conf << (*p_iter)->orientation.q[1] << ' ';
		fout_conf << (*p_iter)->orientation.q[2] << ' ';
		fout_conf << (*p_iter)->orientation.q[3] << ' ';
		fout_conf << (*p_iter)->init_cluster << ' ' ;
		fout_conf << (*p_iter)->wall_connected << ' ';
		if ((*p_iter)->wall){
			fout_conf << 1 << ' ';
		} else {
			fout_conf << 0 << ' ';
		}
		fout_conf << (*p_iter)->valCn_size() << endl;
	}
	int number_of_live_bonds = n_bond - counterBreak ;
	fout_conf << "B " << number_of_live_bonds << endl;
	

	ForAllBond{
		if( (*bond_iter)->status ){
			int i_p0, i_p1;
			int initial_bond;
			
			if ((*bond_iter)->initial_bond){
				initial_bond = 1;
			} else {
				initial_bond = 0;
			}
			(*bond_iter)->whichparticle(i_p0, i_p1);
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

void System::monitorDeformation(char equilibrium){
	fout_deform << "# " << equilibrium;
	fout_deform << ' ' << stress_z ;
	fout_deform << ' ' << stress_x ;
	fout_deform << ' ' << volume_fraction;
	fout_deform << ' ' << wl[0]->z - lz0;
	fout_deform << ' ' << wl[1]->z - lz0;
	fout_deform << ' ' << rup_normal;
	fout_deform << ' ' << rup_shear;
	fout_deform << ' ' << rup_bend;
	fout_deform << ' ' << rup_torsion;
	fout_deform << ' ' << average_contact_number;
	fout_deform << endl;
	int count_bond = 0;
	for(int i=0; i < n_bond; i++){
		if ( bond[i]->initial_bond ){
			count_bond ++;
		}
	}
	fout_deform << "B " << count_bond   <<endl;
	for(int i=0; i < n_bond; i++){
		if ( bond[i]->initial_bond ){
			bond[i]->monitor_state(fout_deform);
		}
	}
}
void System::calcVolumeFraction(){
#ifndef TWODIMENSION
	double volume = lx*ly*lz;
	static const double volume1 = (4.0/3.0)*M_PI;
	volume_fraction = n_particle*volume1 / volume;
#else 
	volume_fraction = n_particle* M_PI / (lx * lz);
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
void System::checkState(){
	
	unsigned long n_particle_active = particle_active.size();
	int sum_num_of_contacts_for_active_particle = 0;
	for (int i=0; i < n_particle_active; i++){
		sum_num_of_contacts_for_active_particle += particle_active[i]->valCn_size();
	}
	average_contact_number = (1.*sum_num_of_contacts_for_active_particle) / n_particle_active;
	
	unsigned long n_bond_active = bond_active.size();
	for (int i=0; i < n_bond_active; i++){
		bond_active[i]->addContactForce();
	}
	double bforce_max = 0;
	double bforce_sum = 0;
	for (int i=0; i < n_bond_active; i++){
		double bforce = bond_active[i]->forceNorm();
		bforce_sum += bforce;
		if ( bforce > bforce_max ){
			bforce_max = bforce ;
		}
	}
	ave_bondforce = bforce_sum / n_bond_active;
	double sum_force = 0.;
	max_force = 0.;
	max_velocity = 0.;
	max_ang_velocity = 0.;
	for (int i=0; i < n_particle_active; i++){
		double force = particle_active[i]->valForce();
		double velocity = particle_active[i]->valVelocity();
		double ang_velocity = particle_active[i]->valOmega();
		sum_force += force;
		if ( force > max_force )
			max_force = force ;
		if ( velocity > max_velocity)
			max_velocity = velocity;
		if ( ang_velocity > max_ang_velocity )
			max_ang_velocity = ang_velocity;
	}
	ave_force = sum_force / n_particle_active;
	for (int i=0; i < n_particle_active; i++){
		particle_active[i]->resetForce();
	}
	calcLocalStrains();
	kinetic_energy = 0;
	for (int i=0; i < n_particle; i++){
		kinetic_energy += particle[i]->kineticEnergy();
	}
}

bool System::checkPercolation(){
	percolation = false;
	for (int i=0; i < n_particle ;i++){
		particle[i]->wall_connected = -1;
	}
	wl[0]->wall_group.clear();
	wl[0]->markWallConnected();
	
	wl[1]->wall_group.clear();
	wl[1]->markWallConnected();

	return percolation;
}

void System::shiftForPercolation(){
	/* This version can deal only 2D
	 * This is only for the cases that a percolated cluster is separated by
	 * top and bottom walls.
	 *
	 * Particles connecting to top wall will be moved.
	 * The fist collision is considered for all other particles.
	 *
	 */
	double lambda_min = 1000;
	vec3d shift;
	if (simulation == 's'){
		shift.set(1,0,0);
	}else if (simulation == 'c'){
		shift.set(0,0,-1);
	}
//	outputConfiguration('e');
	for (int i = 0; i < wl[1]->wall_group.size(); i++){
		int i_uppder = wl[1]->wall_group[i];
		vec3d p_up = particle[i_uppder]->p;
		for (int j = 0; j < n_particle ; j++){
			if (particle[j]->wall_connected != 2){
				vec3d p_lo = particle[j]->p;
				vec3d dp = p_up - p_lo;
				double dot_d_dp = dot(shift, dp );
				double discriminant = dot_d_dp*dot_d_dp - (dp.sq_norm() - 4);
				if (discriminant > 0 ){
					double lambda1 = - dot_d_dp + sqrt(discriminant );
					double lambda2 = - dot_d_dp - sqrt(discriminant );
					if (lambda1>0 && lambda_min > lambda1){
						lambda_min = lambda1;
					}
					if (lambda2>0 && lambda_min > lambda2){
						lambda_min = lambda2;
					}
				}
			}
		}
	}
	double z_highest = 0;
	double z_lowest = 10;
	int cnt = 0;
	for (int i = 0; i < particle.size(); i++){
		if (particle[i]->wall_connected != 2  ){
			cnt ++;
			if (particle[i]->p.z > z_highest){
				cerr << z_highest << endl;
				z_highest = particle[i]->p.z;
			}
		}else {
			if (particle[i]->p.z < z_lowest){
				z_lowest = particle[i]->p.z;
			}
		}
	}

	double lambda_wl_top = wl[1]->z - 1 - z_highest;
	if ( lambda_min > lambda_wl_top){
		lambda_min = lambda_wl_top  ;
	}

	double lambda_wl_bot = z_lowest - (wl[0]->z + 1);
	if ( lambda_min > lambda_wl_bot){
		lambda_min = lambda_wl_bot  ;
	}
	if ( lambda_min != 1000){
		lambda_min += 1e-5; // This avoid to fail new bond generation.
		
		wl[1]->moveWallGroup( lambda_min*shift );
		
		cerr << " lambda_min : " << lambda_min << endl;
	}
	
}

/*
 *	Adjust
 *	- center of mass is set to the center
 *  - Bond vector
 */
void System::simuAdjustment(){
	ForAllParticle{
		(*p_iter)->setNorm_u();
	}
	if (simulation == 'c'){
		lz0 = (wl[1]->z + wl[0]->z) / 2;
	}
}

void System::makeNeighbor(){
	if (simulation == 'c'||simulation == 's'){
		grid->remake_with_walls((wl[0]->z) + 2, (wl[1]->z) - 2, particle);
		ForAllParticle{
			(*p_iter)->makeNeighbor();
		}
		wl[0]->makeNeighbor();
		wl[1]->makeNeighbor();
	} else {
		grid->remake(particle);
		ForAllParticle{
			(*p_iter)->makeNeighbor();
		}
	}
}

void System::setWall()
{
	
	double z_min = lz_init/2;
	double z_max = lz_init/2;
	ForAllParticle{
		if ( (*p_iter)->p.z > z_max){
			z_max = (*p_iter)->p.z;
		}
		if ( (*p_iter)->p.z < z_min){
			z_min = (*p_iter)->p.z;
		}
	}
	n_bot = n_particle;
	n_top = n_particle+1;
	//wl.push_back(new Wall(n_bot, z_min - 1.0, *this)); // wl[0] --> bot
	wl.push_back(new Wall(n_bot, z_min-1.0, *this)); // wl[0] --> bot
	//  cerr << "wall created at z = " << z_min - 1.0 << endl;
	//wl.push_back(new Wall(n_top, z_max + 1.0, *this)); // wl[1] --> top
	wl.push_back(new Wall(n_top, z_max+1.0, *this)); // wl[1] --> top
	//cerr << "wall created at z = " << z_max + 1.0 << endl;
	renew_Lz();
	
	cerr << z_min << endl;
	cerr << z_max << endl;
	cerr << "Lz = " << lz << endl;

}

void System::regeneration(){
	// To be checked
	foreach(vector<int>, regeneration_bond, b){
		bond[*b]->regeneration();
		counterRegenerate ++;
	}
	regeneration_bond.clear();
}

void System::regeneration_onebyone(){
	int most_stressed_bond = regeneration_bond[0];
	double D_max = bond[most_stressed_bond]->D_function;
	unsigned long n_regeneration_bond = regeneration_bond.size();
	if (n_regeneration_bond > 1){
		for (int i=1; i < n_regeneration_bond; i++){
			if ( D_max < bond[regeneration_bond[i]]->D_function ){
				D_max = bond[regeneration_bond[i]]->D_function;
				most_stressed_bond = regeneration_bond[i];
			}
		}
	}
	bond[most_stressed_bond]->regeneration();
	counterRegenerate ++;
	regeneration_bond.clear();
}

void System::rupture(){
	int most_stressed_bond = rupture_bond[0] ;
	double D_max = bond[most_stressed_bond]->D_function;
	unsigned long n_rupture_bond = rupture_bond.size();
	if ( n_rupture_bond > 1){
		for (int i=1; i< n_rupture_bond; i++){
			if ( D_max < bond[ rupture_bond[i] ]->D_function ){
				D_max =  bond[ rupture_bond[i] ]->D_function;
				most_stressed_bond = rupture_bond[i];
			}
		}
	}
	bond[ most_stressed_bond ]->rupture();
	counterBreak ++;

	/* I shold rewrite STL algorithm.
	 *
	 */
	unsigned long i_bondactive = bond[ most_stressed_bond ]->number_activebond;
	bond_active[ i_bondactive ] = bond_active.back();
	bond_active[ i_bondactive ]->number_activebond = i_bondactive;
	bond_active.pop_back();
	rupture_bond.clear();
}

void System::makeInitialBond(double generation_distance){
	double tmp = sq_dist_generate;
	sq_dist_generate = sq(generation_distance);
	makeNeighbor();
	generateBondAll();
	/* z 方向の周期境界条件で生成した初期配置。
	 * 粒子の中心座標は z
	 */
	if (simulation == 'c'||simulation == 's'){
		vector<int> bot_list;
		vector<int> top_list;
		for(int i=0; i < n_particle; i++){
			if ( particle[i]->p.z  <= wl[0]->z + 1.01){
				bot_list.push_back(i);
			}
			if ( particle[i]->p.z  >=  wl[1]->z - 1.01){
				top_list.push_back(i);
			}
		}

		for (int i = 0; i < bot_list.size(); i++){
			wl[0]->initWallParticle(bot_list[i]);
		}
			
		for (int i = 0; i < top_list.size(); i++){
			wl[1]->initWallParticle(top_list[i]);
		}
	}
	sq_dist_generate = tmp;
}

void System::checkBondFailure(){
	foreach( vector<Bond *>, bond_active, it_bond){
		(*it_bond)->cheackBondStress();	
	}
}

void System::outputRestructuring(){

	for (int i = 0; i < n_particle ; i++){
		particle[i]->outputBond();
	}
	
	for (int i = 0; i < regeneration_bond.size(); i++){
		bond[regeneration_bond[i]]->outputRuptureMark();
	}
	for (int i = 0; i < rupture_bond.size(); i++){
		bond[rupture_bond[i]]->outputRuptureMark();
	} 
	cout << endl;
}














