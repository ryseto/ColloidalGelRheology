//
//  contable.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_contable_h
#define ColloidalGelRheology_contable_h
#include <cmath>
#include <cstdlib>
#define DELETE(x) if(x){delete [] x; x = NULL;}

class ConTable {
	bool allocate;
	bool **tbl;
	int n;
public:
	ConTable():allocate(false) {}
	~ConTable();
	bool connect(int i, int j);
	void set(int particleNumber);
	void reset();
	void on_connect(int i, int j);
	void off_connect(int i, int j);
};
#endif
