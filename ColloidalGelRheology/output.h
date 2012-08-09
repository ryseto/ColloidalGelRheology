//
//  output.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_output_h
#define ColloidalGelRheology_output_h
#include <iostream>
#include <fstream>
#include <iomanip>
#include "common.h"
#include "bond.h"
#include "particle.h"
#include "system.h"
using namespace std;


void System::outputYaplot()
{
	static bool firsttime = true;
	if ( firsttime ){firsttime = false;}else{fout_yap << endl;}
	//	fout_yap << "y 9" << endl;
	//	fout_yap << "@ 2" << endl;
	ForAllParticle{
		if( ! (*p_iter)->wall )
			(*p_iter)->output(fout_yap);
	}
	//	fout_yap << "t " << sy.lx/2 + 10 << ' ' << 0 << ' ' << -1 << ' ';
	//	fout_yap << "N= " << sy.n_particle << "  Nb= " << sy.average_contact_number << endl;
	
	if (simulation == 'c'){
		fout_yap << "y 11" << endl;
		fout_yap << "@ 0" << endl;
		wl[0]->output(fout_yap);
		wl[1]->output(fout_yap);
	}
	if (simulation == 'f'){
		fout_yap << "y 11" << endl;
		fout_yap << "@ 0" << endl;
		wl[0]->output(fout_yap);
	}
	
	if (simulation != 'r'){
		fout_yap << "y 12" << endl;
		fout_yap << "@ 3" << endl;
		fout_yap << "t " << lx/2 + 5 << ' ' << 0 << ' ' << 10 << ' ';
		fout_yap << "N= " << n_particle  << " Nb = " << average_contact_number << endl;
		fout_yap << "t " << lx/2 + 5 << ' ' << 0 << ' ' << 5 << ' ';
		fout_yap << "rup: " << rup_normal << ' ' << rup_shear << ' ';
		fout_yap << rup_bend << ' ' << rup_torsion << endl;
		fout_yap << "t " << lx/2 + 5 << ' ' << 0 << ' ' << 0 << ' ';
		fout_yap << "Pz= " << stress_z ;
		fout_yap << "  Px= " << stress_x << endl;
		fout_yap << "t " << lx/2 + 5 << ' ' << 0 << ' ' << -5 << ' ';
		fout_yap << "phi= " << volume_fraction << endl;
	}
	/*
	 fout_yap << "y 10" << endl;
	 fout_yap << "@ 3" << endl;
	 ForAllParticle{
	 (*p_iter)->outputVelocity(fout_yap);
	 }
	 fout_yap << "y 11" << endl;
	 fout_yap << "@ 4" << endl;
	 ForAllParticle{
	 (*p_iter)->outputAngleVelo(fout_yap);
	 }
	 */
}



#endif
