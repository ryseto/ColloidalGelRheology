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

const double pd_threshold = 20;

extern double Lx; // box size
extern double Lz_pd;

class vec2d {
	double xvalue;
	double zvalue;
public:
	vec2d (double x, double z):xvalue(x),zvalue(z){}
	vec2d (void):xvalue(0), zvalue(0){}
	~vec2d(void){}
	inline void set(double x, double z){xvalue=x; zvalue=z;}
	inline double x()const {return xvalue;}
	inline double z()const {return zvalue;}
	inline void shift(double dx, double dz)
	{
		xvalue += dx;
		zvalue += dz;
		
		if (xvalue < 0) {
			xvalue += Lx;
		} else if ( xvalue >= Lx) {
			xvalue -= Lx;
		}
		if (zvalue < 0 ) {
			zvalue += Lz_pd;
		} else if ( zvalue >= Lz_pd) {
			zvalue -= Lz_pd;
		}
	}
	friend vec2d operator+(const vec2d&, const vec2d&);
	friend vec2d operator-(const vec2d&, const vec2d&);
	friend vec2d operator/(const vec2d&, double);
	friend vec2d operator-(const vec2d&);
	friend vec2d operator*(const vec2d&, double);
	friend vec2d operator*(double, const vec2d&);
	friend void operator+=(vec2d &, const vec2d &);
	friend bool operator==(const vec2d &, const vec2d &);
	friend bool operator!=(const vec2d &, const vec2d &);
	friend double dot(const vec2d &a1, const vec2d &a2)
	{
		return a1.xvalue*a2.xvalue + a1.zvalue*a2.zvalue;
	}
	
	friend double dist(const vec2d &a1,
					   const vec2d &a2)
	{
		return sqrt(sq_dist(a1, a2));
	}
	friend double sq_dist(const vec2d &p1, const vec2d &p2)
	{
		return (p1.xvalue-p2.xvalue)*(p1.xvalue-p2.xvalue) +
		(p1.zvalue-p2.zvalue)*(p1.zvalue-p2.zvalue);
	}
	friend double sq_dist_pd(const vec2d &p1, const vec2d &p2)
	{
		double dx = p1.xvalue - p2.xvalue;
		double dz = p1.zvalue - p2.zvalue;
		if (dx > pd_threshold) {
			dx -= Lx;
		} else if (dx < -pd_threshold) {
			dx += Lx;
		}
		
		if (dz > pd_threshold) {
			dz -= Lz_pd;
		} else if (dz < - pd_threshold){
			dz += Lz_pd;
		}
		return dx*dx + dz*dz;
	}

	friend vec2d diff_pd(const vec2d &p1, const vec2d &p2){
		double dx = p1.xvalue - p2.xvalue;
		double dz = p1.zvalue - p2.zvalue;
		if (dx > pd_threshold) {
			dx -= Lx;
		} else if (dx < -pd_threshold) {
			dx += Lx;
		}
		if (dz > pd_threshold) {
			dz -= Lz_pd;
		} else if (dz < -pd_threshold) {
			dz += Lz_pd;
		}
		return vec2d(dx, dz);
	}
	
	inline friend vec2d n_vec(vec2d, vec2d);
	inline friend vec2d x_plus(vec2d v, double b)
	{
		return vec2d( v.xvalue + b, v.zvalue);
	}
	friend double angle(vec2d, vec2d);
//	inline vec2d pd(double xo, double L) {
//		vec2d p_periodic;
//		if (xo < L/2){
//			if (xvalue > xo + Lx/2 )
//				p_periodic.set(xvalue - Lx, zvalue);
//			else
//				p_periodic.set(xvalue, zvalue);
//		} else {
//			if (xvalue < xo - L/2 )
//				p_periodic.set(xvalue + Lx, zvalue);
//			else
//				p_periodic.set(xvalue, zvalue);
//		}
//		return p_periodic;
//	}
	inline vec2d plus(double dx, double dz)
	{
		return vec2d(xvalue + dx, zvalue + dz);
	}
	void periodic_range(double Lx_);
	void periodic_range_xz(double Lx_, double Lz_pd_);
	
	double abs();
	double sq_norm(){return xvalue*xvalue + zvalue*zvalue;}
	void cout();
	void cerr()
	{
		std::cerr << xvalue << ' ' << zvalue << std::endl;
	}
	void out_data();
};

#endif /* defined(__ColloidalGelRheology__vec2d__) */
