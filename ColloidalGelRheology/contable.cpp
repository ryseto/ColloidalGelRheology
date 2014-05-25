//
//  contable.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

//#include <iostream>
#include "contable.h"
#include <iostream>
#include <iomanip>
using namespace std;

ConTable::~ConTable()
{
	if (allocate) {
		for (int i=0; i<n; i++) {
			DELETE(tbl[i]);
		}
		DELETE(tbl);
	}
}

void ConTable::set(int particleNumber)
{
	allocate = true;
	n = particleNumber;
	tbl = new bool* [n];
	for (int i=0; i<n; i++) {
		tbl[i] = new bool[n];
	}
	for (int i=0; i<n; i++) {
		for (int j=0; j<i; j++) {
			tbl[i][j] = false;
			tbl[j][i] = false;
		}
		tbl[i][i] = true;
	}
}

bool ConTable::connect(int i, int j)
{
	return tbl[i][j];
}

void ConTable::on_connect(int i, int j)
{
	tbl[i][j] = true;
	tbl[j][i] = true;
}

void ConTable::off_connect(int i, int j)
{
	tbl[i][j] = false;
	tbl[j][i] = false;
}

void ConTable::reset()
{
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			tbl[i][j] = false;
		}
		tbl[i][i] = true;
	}
}
