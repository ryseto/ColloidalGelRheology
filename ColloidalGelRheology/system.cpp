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
	sprintf(version ,"0");
	dt = 0.;
	stress_z = 0.;
	stress_x = 0.;
	strain_x = 0.;
    strain_z = 0.;
	max_force = 0.;
	d_strain_x = 0;
	d_strain_z = 0;
    ct = new ConTable;
    grid = new Grid;
    //    modified_area_fraction = -1;
}

System::~System(){
	//fout_yap.close();
    fout_conf.close();
	fout_data.close();
	fout_log.close();
    fout_deform.close();
    delete ct;
    delete grid;
}

void System::TimeDevStrainControlCompactionEuler(vector<Particle *> &particle_active, 
                                                 vector<Bond *> &bond_active,
                                                 bool wall_move){
    if(wall_move){
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
}

void System::TimeDevStrainControlShearEuler(vector<Particle *> &particle_active, 
                                            vector<Bond *> &bond_active,
                                            bool wall_move){
    if(wall_move){
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
}


void System::setParameterFile(char *parameters_file_)
{
	sprintf(parameters_file, "%s", parameters_file_);
	string s_parameters_file = parameters_file;
	int i_backslash = s_parameters_file.find_last_of( "/") + 1;
	int i_extention = s_parameters_file.find( ".txt" );	
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
				parameters, init_cluster, version);
        //	} else if(simulation == 'f'){
        //		sprintf(name_string, "CF%s_%s_H%d_%s",
        //				parameters, init_cluster, (int)lz_init, version);		
	} else if(simulation == 't'){
        sprintf(name_string, "doublet");
    }
    
	string fn_common = (string)name_string ;
	sprintf(fn_snapshot, "ss_%s.yap", fn_common.c_str()); 
	sprintf(fn_data, "data_%s.dat", fn_common.c_str()); 
	sprintf(fn_log, "log_%s.dat", fn_common.c_str());	
	sprintf(fn_conf, "conf_%s.dat", fn_common.c_str());
    sprintf(fn_deform, "deform_%s.dat", fn_common.c_str());
    
    if ( simulation == 't')
        fout_yap.open(fn_snapshot);    
	fout_data.open(fn_data);	
	fout_log.open(fn_log);
	fout_conf.open(fn_conf);
    fout_deform.open(fn_deform);
}

void System::setInitClusterFile(char *init_cluster_file_)
{
	sprintf(init_cluster_file, "%s", init_cluster_file_);
	string s_init_cluster_file = init_cluster_file;
	int i_backslash = s_init_cluster_file.find_last_of( "/") + 1;
	int i_extention = s_init_cluster_file.find( ".dat" );
	sprintf(init_cluster, "%s", 
			(s_init_cluster_file.substr( i_backslash, i_extention-i_backslash)).c_str());
    //    int i_L = s_init_cluster_file.find("_L")+2;
    //  int j_L = s_init_cluster_file.find( "_", i_L );
    //    sprintf(string_L, "%s", 
    //			(s_init_cluster_file.substr( i_L, j_L - i_L)).c_str());
    
	//string input = s_init_cluster_file.substr(i_backslash, i_extention-i_backslash);
}

/*
 * to distinguish the output files :
 * file name is FileName_***.dat
 * "***" is given by command argument. 
 * Default: version = 0
 */ 
void System::setVersion(char *version_)
{
	sprintf(version, "%s", version_);
}



void System::setBondGenerationDistance(double distance)
{
	dist_generate = distance;
	sq_dist_generate = sq(distance);
}

void System::generateBond(vector<Bond *> &bond_active){
	int bond_number_before = n_bond;
	for (int i=0; i < n_particle; i++){
		particle[i]->generateBond();
	}
    if (bond_number_before == n_bond){
        return;
    }
	for (int i = bond_number_before; i < n_bond; i++){
		bond_active.push_back( bond[i] );
        bond[i]->number_activebond  = bond_active.size() - 1;
	}
}

void System::generateBond(){
	for (int i=0; i < n_particle; i++){
		particle[i]->generateBond();
	}
}

void System::prepareCalculationParameters()
{
    double k_max = 0;
    if (k_max < bond0.kn )
        k_max = bond0.kn;
    if (k_max < bond0.ks )
        k_max = bond0.ks;
    if (k_max < bond0.kb )
        k_max = bond0.kb;
    if (k_max < bond0.kt )
        k_max = bond0.kt;
    if (k_max < bond1.kn )
        k_max = bond1.kn;
    if (k_max < bond1.ks )
        k_max = bond1.ks;
    if (k_max < bond1.kb )
        k_max = bond1.kb;
    if (k_max < bond1.kt )
        k_max = bond1.kt;     
    
    dt_max = (1./50.)*(1./sqrt(k_max));
    cerr << "dt_max = " << dt_max  << endl;
    dt = dt_max;
    eta = 0.05*bond0.c_bend;
    //eta = bond0.c_bend;
    eta_rot = (4.0/3.0)*eta;
    if ( simulation == 'c' ||  simulation == 's'){
		lx0 = lx/2;
		ly0 = ly/2;
        lz0 = lz_init/2;
	}
}

void System::optimalTimeStep(){
    //double dt_new = sqrt( max_move_step /  max_force );
    double dt_new = max_move_step /  max_velocity;
    if ( dt_new > dt_max ){
        dt_new = dt_max;
    } else {
        cerr << "small time step"  << dt_new << endl;
    }
    
    
    dt = dt_new;
}

void System::calcStress(){
    double stress_z_bot = wl[0]->stressSensor_z();
    double stress_z_top = wl[1]->stressSensor_z();
    stress_z = 0.5*(abs(stress_z_bot)+abs(stress_z_top));
    diff_stress_z = abs( stress_z_bot + stress_z_top )/stress_z;
}

void System::calcShearStress(){
    double stress_x_bot = wl[0]->stressSensor_x();
    double stress_x_top = wl[1]->stressSensor_x();
    stress_x = 0.5*(abs(stress_x_bot) + abs(stress_x_top));
    if (stress_x < 1e-3){
        diff_stress_x = 0;
    }else{
        diff_stress_x = abs(stress_x_bot + stress_x_top )/stress_x;
    }
    cerr << stress_x << ' ' << diff_stress_x << endl;
    //    cerr << stress_x << endl;
    double stress_z_bot = wl[0]->stressSensor_z();
    double stress_z_top = wl[1]->stressSensor_z();
    stress_z = 0.5*(abs(stress_z_bot)+abs(stress_z_top));
    diff_stress_z = abs( stress_z_bot + stress_z_top )/stress_z;
}


void prepareBond(BondParameter & _bond){
    double a0 = 1.0;
    if (_bond.fnc != 0){
        _bond.kn = _bond.fnc / _bond.n_max;
        _bond.ks = _bond.fsc / _bond.s_max;
        _bond.kb = _bond.mbc / (_bond.b_max * a0 * a0);
        _bond.kt = _bond.mtc / (_bond.t_max * a0 * a0);
        _bond.kn3 = _bond.overlap_force_factor * _bond.kn / sq(_bond.n_max);
        
        _bond.c_norm = 2*sqrt( _bond.kn);
        _bond.c_slid = 2*sqrt( _bond.ks);
        _bond.c_bend = 2*sqrt( (2./5.)*_bond.kb);
        _bond.c_tort = 2*sqrt( (2./5.)*_bond.kt);
        //_bond.c_bend = 2*sqrt( _bond.kb);
        //_bond.c_tort = 2*sqrt( _bond.kt);
        cerr << "k:" << _bond.kn << ' ' << _bond.ks << ' ' << _bond.kb << ' ' << _bond.kt << endl;
        cerr << "c:" << _bond.c_norm << ' ' << _bond.c_slid << ' ' << _bond.c_bend << ' ' << _bond.c_tort << endl;
    } else {
        _bond.fnc = 0.;
        _bond.kn = 1.0 / _bond.n_max;
        _bond.ks = 0;
        _bond.kb = 0;
        _bond.kt = 0;
        _bond.kn3 = _bond.overlap_force_factor * _bond.kn / sq(_bond.n_max);
        _bond.c_norm = 2*sqrt( _bond.kn);
        _bond.c_slid = 0;
        _bond.c_bend = 0;
        _bond.c_tort = 0;
    }
}

void setBondParameter( BondParameter &bondparameter, string &bond_file, string &dir){
    string path_bond = dir + "/" + bond_file;
    ifstream fin;
    if(fin.is_open())
        fin.close();
    fin.open(path_bond.c_str());
    string codeword;
    fin >> codeword >> bondparameter.fnc;
    fin >> codeword >> bondparameter.fsc;
    fin >> codeword >> bondparameter.mbc;
    fin >> codeword >> bondparameter.mtc;
    fin >> codeword >> bondparameter.n_max;
    fin >> codeword >> bondparameter.s_max;
    fin >> codeword >> bondparameter.b_max;
    fin >> codeword >> bondparameter.t_max;
    fin >> codeword >> bondparameter.overlap_force_factor;
    fin >> codeword >> bondparameter.robust_bond_compression;
    fin.close();
    prepareBond(bondparameter);
    return;
}

void System::readBondParameter(){
    setBondParameter(bond0, bond0_file, rootdir);
    setBondParameter(bond1, bond1_file, rootdir);
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
    
    keylist["rootdir:"]=1; const int _rootdir = 1;
    keylist["bond0_file:"]=2; const int _bond0_file = 2;
    keylist["bond1_file:"]=3; const int _bond1_file = 3;
    keylist["lx:"]=9; const int _lx = 9;
    keylist["lz_init:"]=10; const int _lz_init = 10;
    keylist["initial_compaction:"]=11; const int _initial_compaction=11;
    keylist["volumefraction_increment:"]=12; const int _volumefraction_increment=12;
	keylist["max_volume_fraction:"]=13; const int _max_volume_fraction=13;
    keylist["wall_velocity:"]=20;const int _wall_velocity=20;
    //    keylist["dt_max:"]=21; const int _dt_max=21;
    keylist["max_move_step:"]=22; const int _max_move_step=22;
    keylist["relax_for_restructuring:"]=30; const int _relax_for_restructuring=30;
    //    keylist["tolerance_convergence:"]=31; const int _tolerance_convergence=31;
	keylist["interval_convergence_check:"]=31; const int _interval_convergence_check=31;
    keylist["max_velocity_convergence:"]=32; const int _max_velocity_convergence=32;
    keylist["max_ang_velocity_convergence:"]=33; const int _max_ang_velocity_convergence=33;
    keylist["diff_stress_convergence:"]=34;const int _diff_stress_convergence=34;    
    keylist["stress_change_convergence:"]=35;const int _stress_change_convergence=35;
	cerr << codeword << ' ' << value << endl;
	switch(keylist[codeword]){
        case _rootdir: rootdir = value; break;
        case _bond0_file: bond0_file = value; break;
        case _bond1_file: bond1_file = value; break;
        case _lx: lx = atof(value.c_str()); break;
        case _lz_init: lz_init = atof( value.c_str()) ; break;
            //		case _dt_max: dt_max = atof( value.c_str()) ; break;
        case _initial_compaction: initial_compaction = atof(value.c_str()) ; break;
            //case _tolerance_convergence: tolerance_convergence = atof(value.c_str()); break;
		case _interval_convergence_check: interval_convergence_check = atoi(value.c_str()); break;
        case _max_velocity_convergence: max_velocity_convergence = atof(value.c_str()); break;
        case _wall_velocity: wall_velocity = atof(value.c_str()); break;
        case _max_ang_velocity_convergence: max_ang_velocity_convergence = atof(value.c_str());; break;
		case _max_volume_fraction: max_volume_fraction = atof(value.c_str());; break;
		case _max_move_step: max_move_step = atof(value.c_str()); break;
        case _volumefraction_increment: volumefraction_increment = atof(value.c_str());; break;
        case _diff_stress_convergence: diff_stress_convergence = atof(value.c_str());; break;
        case _stress_change_convergence: stress_change_convergence = atof(value.c_str());; break;
        case _relax_for_restructuring: relax_for_restructuring = atoi(value.c_str()); break;
		default:
			cerr << "The codeword " << codeword << " is'nt associated with an parameter" << endl;
			exit(1);
	}
}

void System::importPositions()
{	
    string path = rootdir + "/" + init_cluster_file;
	vector<vec3d> init_aggregate;
    vector<int> init_aggregate_cluster;
	ifstream fin;
    fin.open( path.c_str());
    if (! fin.is_open() ){
        cerr << "no initial file" << endl;
        exit(1);
    }
    cerr << "lx = " << lx  << endl;
#ifdef TWODIMENSION
    ly = 2.;
#else 
    ly = lx;
#endif
    
	double x, y, z;
    int i_cluster;
	if (simulation == 'c' || simulation == 's'){
		do{			
#ifdef TWODIMENSION	
			fin >> x >> z >> i_cluster ;
			y = 1.0;
#else
			fin >> x >> y >> z;
#endif
            vec3d new_p(x, y, z);
            init_aggregate.push_back(new_p);
            init_aggregate_cluster.push_back(i_cluster);
		}while (!fin.eof());
	}
	
    //	init_aggregate.pop_back();
	int count_particle_number = 0;
    
    //	while ( !init_aggregate.empty() ){
    for (int i = 0 ; i < init_aggregate.size()-1 ; i++){
		//vec3d new_p = init_aggregate.back();
        vec3d new_p = init_aggregate[i];
        //init_aggregate.pop_back();
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

void System::setRod(int m){
	
	for (int i = 0 ; i < m ; i++){
		vec3d new_p(0, 0, 2*i - m);
		particle.push_back(new Particle(i, new_p, 0, *this) );
	}
	n_particle = m;
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
		p_iter->y += ly/2;
		p_iter->z += lz/2;
	}
    
}

/* Checking active bonds at the first step.
 *
 *
 */
void System::check_active(vector<Particle *> &particle_active, vector<Bond *> &bond_active){
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
    r_min = 2.;
    r_max = 2.;
    bending_angle_max = 0.;
    torsional_angle_max = 0.;
    sliding_disp_max = 0.;
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
	grid->init(n_particle, lx, ly, lz_init, cell_h, simulation);
}

void System::setWallStress(double value_x, double value_y,  double value_z){
	if (value_x == -1){
		wl[0]->x_movable = false;
		wl[1]->x_movable = false;
		stress_x = 0;
	} else {
		wl[0]->x_movable = true;
		wl[1]->x_movable = true;
		stress_x = value_x;
	}
	if (value_y == -1){
		wl[0]->y_movable = false;
		wl[1]->y_movable = false;
		stress_y = 0;
	} else {
		wl[0]->y_movable = true;
		wl[1]->y_movable = true;
		stress_y = value_y;
	}
	
	if (value_z == -1){
		wl[0]->z_movable = false;
		wl[1]->z_movable = false;
		stress_z = 0;
	} else {
		wl[0]->z_movable = true;
		wl[1]->z_movable = true;
		stress_z = value_z;
	}
    
	double fx_wall = stress_x*lx*ly;
	double fy_wall = stress_y*lx*ly;
	double fz_wall = stress_z*lx*ly;
	force_wall.set(-fx_wall, -fy_wall, fz_wall);
	cerr << "stress = " << stress_z << ',' << stress_x << " Pa ";
	cerr << " Fwall=" << fz_wall << ',' << fx_wall << " N" << endl;
	
}

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
    double del_strain;
    
	if ( firsttime ){
		firsttime = false;
        fout_data << "#1 volume_fraction\n";
		fout_data << "#2 stress_z\n";
        fout_data << "#3 stress_x\n";
        fout_data << "#4 strain_z\n";
        fout_data << "#5 strain_x\n";
        fout_data << "#6 area_fraction\n";
		fout_data << "######  rupture: total number \n";
		fout_data << "#7  rup_normal(total number)\n";
		fout_data << "#8  rup_shear(total number)\n";
		fout_data << "#9  rup_bend(total number)\n";
		fout_data << "#10 rup_torsion(total number)\n";
		fout_data << "######  rupture: rate\n";
		fout_data << "#11 rate_rup_normal\n";
		fout_data << "#12 rate_rup_shear\n";
		fout_data << "#13 rate_rup_bend\n";
		fout_data << "#14 rate_rup_torsion\n";
		fout_data << "###### \n";
		fout_data << "#15 average_contact_number\n";
		fout_data << "#16 del_strain\n";
	}
	//double previous_volume_fraction = volume_fraction_equilibrium.back() ;
    //double compressive_strain = (volume_fraction - previous_volume_fraction)/previous_volume_fraction;
    
    if (simulation == 'c'){
        del_strain = (strain_z - strain_previous);
        strain_previous = strain_z;
    } else if (simulation == 's'){
        del_strain = (strain_x - strain_previous);
        strain_previous = strain_x;
    }
    
	double rate_rup_normal = (rup_normal - rup_normal_previous)/del_strain;
	double rate_rup_shear = (rup_shear - rup_shear_previous)/del_strain;
	double rate_rup_bend = (rup_bend - rup_bend_previous)/del_strain;
	double rate_rup_torsion = (rup_torsion - rup_torsion_previous)/del_strain ;
	fout_data << volume_fraction << ' ' ; //1
    fout_data << stress_z << ' ' ;  // 2
	fout_data << stress_x << ' ' ; // 3
    fout_data << strain_z << ' ' ; // 4
    fout_data << strain_x << ' ' ; // 5
    fout_data << area_fraction << ' '; // 6
    ///////////////////////////////////////
	fout_data << rup_normal << ' '; //7
	fout_data << rup_shear << ' ';  //8
	fout_data << rup_bend << ' ';   //9
	fout_data << rup_torsion << ' ';//10
	///////////////////////////////////////
	fout_data << rate_rup_normal << ' '; //11
	fout_data << rate_rup_shear << ' ';  //12
	fout_data << rate_rup_bend << ' ';   //13
	fout_data << rate_rup_torsion << ' ';//14
	///////////////////////////////////////
	fout_data << average_contact_number << ' '; //15
	fout_data << del_strain << ' '; // 16
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
        fout_log << "#27 d_strain_z\n";
        fout_log << "#28 d_strain_x\n"; 
        fout_log << "#29 calc_count\n";
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
    fout_log << d_strain_z << ' '; // 27
    fout_log << d_strain_x << ' '; // 28
    fout_log << calc_count << ' '; // 29
	fout_log << endl; 
}

void System::outputConfiguration(char equilibrium){
    lz0 = (wl[1]->z + wl[0]->z) / 2;
    
	fout_conf << "# " << equilibrium;
	fout_conf << ' ' << volume_fraction;
	fout_conf << ' ' << area_fraction;
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
        fout_conf << (*p_iter)->wall << endl;
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
            fout_conf << (*bond_iter)->val_M_tors() << endl;
            //fout_bconf << (*bond_iter)->orientation.q[3] << endl;
        }
    }
}

void System::monitorDeformation(char equilibrium){
    fout_deform << "# " << equilibrium;
    fout_deform << ' ' << stress_z ;
    fout_deform << ' ' << stress_x ;
    fout_deform << ' ' << volume_fraction;
    fout_deform << ' ' << area_fraction;
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

//double partialSphereVolume(double lack_z){
//    double volume = M_PI*(4.0/3 + lack_z*lack_z*lack_z/3 - lack_z*lack_z);
//    return volume;    
//}
//
//double partialDiskArea(double lack_z){
//    double area = 0.5*(M_PI + 2.*(1. - lack_z)*sqrt(lack_z*(2.-lack_z))  \
//                         + 2.*asin(1.-lack_z));
//    return area;    
//}


//void System::calcModifiedVolumeFraction(){
//    /* Due to the implimentation of non-slip walls,
//     * the lower dencity aound wall should be removed.
//     * However, I stop to modify in the main simulation.
//     */ 
//    int n_inside = 0;
//    double partial_volume = 0.;
//#ifdef TWODIMENSION
//    double partial_area = 0.;
//#endif
//    double margin_cutoff = 3.;
//    double botlevel = wl[0]->z + margin_cutoff;
//    double toplevel = wl[1]->z - margin_cutoff;
//    double lack_z;
//    for (int i = 0; i < n_particle ; i++){
//        double z = particle[i]->p.z;
//        if ( z > botlevel + 1 && z < toplevel - 1){
//            n_inside ++;
//        } else if ( z > botlevel - 1 && z < botlevel +1){
//            lack_z = botlevel + 1 - z;
//            partial_volume += partialSphereVolume(lack_z);
//#ifdef TWODIMENSION
//            partial_area += partialDiskArea(lack_z);
//#endif
//        } else if ( z < toplevel + 1 && z > toplevel -1){
//            lack_z = z - (toplevel - 1);
//            partial_volume += partialSphereVolume(lack_z);
//#ifdef TWODIMENSION
//            partial_area += partialDiskArea(lack_z);
//#endif
//        }   
//    }
//    double solid_volume = n_inside*M_PI*4.0/3 + partial_volume;
//    modified_volume_fraction = solid_volume / (lx*ly*(toplevel-botlevel));
//#ifdef TWODIMENSION
//    double solid_area = n_inside*M_PI + partial_area;
//    modified_area_fraction = solid_area / (lx * (toplevel-botlevel));
//#endif
//}

void System::calcVolumeFraction(){
    double volume = lx*ly*lz;
    static const double volume1 = (4.0/3.0)*M_PI;
    volume_fraction = n_particle*volume1 / volume;
#ifdef TWODIMENSION
    area_fraction = n_particle* M_PI / (lx * lz);
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
void System::checkState(vector <Particle *> &particle_active, 
                        vector <Bond *> &bond_active){
    int n_particle_active = particle_active.size();
    int sum_num_of_contacts_for_active_particle = 0;
    for (int i=0; i < n_particle_active; i++){
        sum_num_of_contacts_for_active_particle += particle_active[i]->valCn_size();
    }
    average_contact_number = (1.*sum_num_of_contacts_for_active_particle) / n_particle_active;
    //
    //////////////////////////////////////////////////////////////////////////
    // Calc forces
    int n_bond_active = bond_active.size();
	for (int i=0; i < n_bond_active; i++){
        bond_active[i]->addContactForce();
	}
	double bforce_max = 0;
	double bforce_sum = 0;
	for (int i=0; i < n_bond_active; i++){
        double bforce = bond_active[i]->strenghOfForce();
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
    //
    //////////////////////////////////////////////////////////////////////////
    //
    calcLocalStrains();
    //
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
}

bool System::checkPercolation(){
    percolation = false;    
    for (int i=0; i < n_particle ;i++){
        particle[i]->wall_connected = -1;
    }
    wl[0]->wall_group.clear();
    wl[1]->wall_group.clear();
    wl[0]->markWallConnected();
    wl[1]->markWallConnected();
    return percolation;    
}

void System::shiftForPercolation(){
    //This version can deal only 2D
    double del_x_min = 1000;
    cerr << wl[0]->wall_group.size() << ' ' << wl[1]->wall_group.size() << endl; 
    
    for (int i = 0; i < wl[1]->wall_group.size(); i++){
        int i_uppder = wl[1]->wall_group[i]; 
        vec3d p_up = particle[i_uppder]->p;
        for (int j = 0; j < wl[0]->wall_group.size(); j++){
            int i_lower = wl[0]->wall_group[j];
            vec3d p_lo = particle[i_lower]->p;
            if (p_lo.x < p_up.x){
                p_lo.x += lx; 
            }            
            // p1.x - p0.x > 0
            double yy_zz = sq(p_lo.y - p_up.y) + sq(p_lo.z - p_up.z);
            if (4 - yy_zz > 0){
                double del_x = p_lo.x - p_up.x - sqrt( 4 - yy_zz);
                //                cerr << del_x << endl;
                if ( del_x < del_x_min){
                    del_x_min = del_x;
                }
            }            
        }
        //cerr << i << " " << del_x_min << endl;
    }
    
    for (int i = 0; i < wl[1]->wall_group.size(); i++){
        int k = wl[1]->wall_group[i];
        particle[k]->x_shift(del_x_min);
    }
    cerr << " del_x_min : " << del_x_min << endl;
}



/*
 *	Adjust 
 *	- center of mass is set to the center
 *  - Bond vector 
 */ 
void System::simuAdjustment(){
	//shiftCenterOfMass(particle);
	ForAllParticle{
		(*p_iter)->setNorm_u();
	}
    if (simulation == 'c'||simulation == 's'){
        lz0 = (wl[1]->z + wl[0]->z) / 2;
    }
}

void System::makeNeighbor(){
    if (simulation == 'c'||simulation == 's'){
        double close_to_bot_wall = (wl[0]->z) + 2.;
        double close_to_top_wall = (wl[1]->z) - 2.;
        
        grid->remake_with_walls(close_to_bot_wall, close_to_top_wall, particle);
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
    double z_min = 100, z_max = 0;
    ForAllParticle{
        if ( (*p_iter)->z() > z_max){z_max = (*p_iter)->z();}
        if ( (*p_iter)->z() < z_min){z_min = (*p_iter)->z();}
    }
    //  cerr << "z_max " << z_max << endl;
    //    cerr << "z_min " << z_min << endl;
    n_bot = n_particle;
    n_top = n_particle+1;
    wl.push_back(new Wall(n_bot, z_min - 1.0, *this)); // wl[0] --> bot
    //wl.push_back(new Wall(n_bot, 0, *this)); // wl[0] --> bot
    //  cerr << "wall created at z = " << z_min - 1.0 << endl;
    wl.push_back(new Wall(n_top, z_max + 1.0, *this)); // wl[1] --> top
    //wl.push_back(new Wall(n_top, lz_init, *this)); // wl[1] --> top
    //    cerr << "wall created at z = " << z_max + 1.0 << endl;
    renew_Lz();
    cerr << "Lz_init = " << lz << endl;
}

void System::regeneration(){
    // To be checked
	foreach(vector<int>, regeneration_bond, b){
		bond[*b]->regeneration();
		counterRegenerate ++;
	}
}

void System::regeneration_onebyone(){
    double D_max = 0.;
    int most_stressed_bond = -1;
    int n_regeneration_bond = regeneration_bond.size();
    for (int i=0; i < n_regeneration_bond; i++){
        if ( D_max < bond[ regeneration_bond[ i ] ]->D_function ){
            D_max =  bond[ regeneration_bond[ i ] ]->D_function;
            most_stressed_bond = regeneration_bond[ i ];
        }
    }
    bond[ most_stressed_bond  ]->regeneration();
    counterRegenerate ++;
}

void System::rupture(vector<Bond *> &bond_active){
    // To be check
	foreach(vector<int>, rupture_bond, b){
		bond[*b]->rupture();	
		counterBreak ++;
	}
	/*	
	 * This is a dangerous trick.
	 * I should change this part more suitable way.
	 */
    int n_bond_active = bond_active.size() ;
	for (int i = 0 ; i < n_bond_active; i++){
		if ( bond_active[i]->status == 0 ){
			bond_active[i] = bond_active.back();
            bond_active[i]->number_activebond = i;
			bond_active.pop_back();
            n_bond_active --;
			i --;
		}
	}    
}

void System::rupture_onebyone(vector<Bond *> &bond_active){
    double D_max = 0;
    int most_stressed_bond = rupture_bond[0] ;
    int n_rupture_bond =  rupture_bond.size();
    for (int i=0; i< n_rupture_bond; i++){
        if ( D_max < bond[ rupture_bond[ i ] ]->D_function ){
            D_max =  bond[ rupture_bond[ i ] ]->D_function;
            most_stressed_bond = rupture_bond[ i ];
        }
    }
    bond[ most_stressed_bond ]->rupture();
    counterBreak ++;
    int i_bondactive = bond[ most_stressed_bond ]->number_activebond;
    bond_active[ i_bondactive ] = bond_active.back();
    bond_active[ i_bondactive ]->number_activebond = i_bondactive;
    bond_active.pop_back();
}

void System::makeInitialBond(double generation_distance){
	double tmp = sq_dist_generate;
	sq_dist_generate = sq(generation_distance);
	makeNeighbor();
	generateBond();	
    
	if (simulation == 'c'||simulation == 's'){
        vector<int> bot_list;
        vector<int> top_list;
        for(int i=0; i < n_particle; i++){
            if ( particle[i]->p.z  <= 3 ){
                bot_list.push_back(i);
            }
            if ( particle[i]->p.z  >= lz_init - 3 ){
                top_list.push_back(i);
            }
        }
        for (int i = 0; i < bot_list.size(); i++){
            for (int j = 0; j < top_list.size(); j++){
                vec3d p_bot = particle[bot_list[i]]->p;
                vec3d p_top = particle[top_list[j]]->p;
                p_bot.z += lz_init - 2;
                vec3d del = p_bot - p_top;
                if ( del.sq_norm() <= 4.0001 ){
                    wl[0]->initWallParticle(bot_list[i]);
                    wl[1]->initWallParticle(top_list[j]);
                }
            }
        }
	}
	sq_dist_generate = tmp;
}

void System::checkBondFailure(vector<Bond *> &bond_active){
	foreach( vector<Bond *>, bond_active, it_bond){
        (*it_bond)->failureCondition();	
	}
}
