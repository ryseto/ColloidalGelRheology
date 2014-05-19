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
	if (firsttime) {
		firsttime = false;
	} else {
		fout_yap << endl;
	}
	//	fout_yap << "y 9" << endl;
	//	fout_yap << "@ 2" << endl;
	ForAllParticle {
		(*p_iter)->output(fout_yap);
	}
	//	fout_yap << "t " << sy.lx/2 + 10 << ' ' << 0 << ' ' << -1 << ' ';
	//	fout_yap << "N= " << sy.n_particle << "  Nb= " << sy.average_contact_number << endl;
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
