//
//  vec2d.h
//  ColloidalGelRheology
//
//  Created by Ryohei Seto on 2012/08/17.
//
//

#ifndef __ColloidalGelRheology__vec2d__
#define __ColloidalGelRheology__vec2d__

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
const double L_SIZE = 5;
extern double L;

class vec2d {
	double xvalue;
	double yvalue;
public:
	vec2d (double x, double y):xvalue(x),yvalue(y){}
	vec2d (void):xvalue(0), yvalue(0){}
	~vec2d(void){}
	inline void set(double x, double y){xvalue=x;yvalue=y;}
	inline double x()const {return xvalue;}
	inline double y()const {return yvalue;}
	inline double x(double L);
	inline void shift(double x, double y){
		xvalue += x;
		yvalue += y;
	}
	friend vec2d operator+(vec2d, vec2d);
	friend vec2d operator-(vec2d, vec2d);
	friend vec2d operator/(vec2d, double);
	friend vec2d operator-(vec2d);
	friend vec2d operator*(vec2d, double);
	friend vec2d operator*(double, vec2d);
	friend void operator+=(vec2d &, vec2d );
	friend bool operator==(vec2d, vec2d);
	friend bool operator!=(vec2d, vec2d);
	friend double dot(vec2d, vec2d);
	friend double dist(vec2d a1, vec2d a2);
	friend double sq_dist(vec2d a1, vec2d a2);
	friend double sq_dist_pd(vec2d a1, vec2d a2);
	
	
	inline friend vec2d n_vec(vec2d, vec2d);
	inline friend vec2d x_plus(vec2d v, double b){
		return vec2d( v.xvalue + b, v.yvalue);}
	friend double angle(vec2d, vec2d);
	inline vec2d pd(double xo, double L) {
		vec2d p_periodic;
		if (xo < L/2){
			if (xvalue > xo + L/2 )
				p_periodic.set(xvalue-L, yvalue);
			else
				p_periodic.set(xvalue, yvalue);
		} else {
			if (xvalue < xo - L/2 )
				p_periodic.set(xvalue+L, yvalue);
			else
				p_periodic.set(xvalue, yvalue);
		}
		return p_periodic;
	}
	inline vec2d plus(double dx, double dy){
		return vec2d(xvalue + dx, yvalue + dy);
	}
	void periodic_range(double L);
	void periodic_range_xy(double Lx, double Ly);
	
	double abs();
	double sq_norm(){return xvalue*xvalue + yvalue*yvalue;}
	void cout();
	void cerr(){
		std::cerr << xvalue << ' ' << yvalue << std::endl;
		
	}
	void out_data();
};

#endif /* defined(__ColloidalGelRheology__vec2d__) */
