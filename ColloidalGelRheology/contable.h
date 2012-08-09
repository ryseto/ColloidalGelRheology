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
using namespace std;

class ConTable{
	bool allocate;
	bool **tbl;
	int n;
public:
	ConTable():allocate(false) {}
	~ConTable();
	inline bool connect(int i, int j){
		return tbl[i][j];
	}
	
	void set(int particleNumber);
	void reset();
	void on_connect(int i, int j);
	void off_connect(int i, int j);
};
#endif
