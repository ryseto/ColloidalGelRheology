//
//  vec2d.cpp
//  ColloidalGelRheology
//
//  Created by Ryohei Seto on 2012/08/17.
//
//

#include "vec2d.h"
#include "vec2d.h"

inline double sq(double x){
	return x*x;
}
//
//double vec2d::x(double L){
//	if (xvalue < 0)
//		return xvalue + L;
//	else if (xvalue > L)
//		return xvalue - L;
//	else
//		return xvalue;
//}

vec2d operator+(const vec2d &a1, const vec2d &a2){
	return vec2d( a1.xvalue + a2.xvalue, a1.zvalue + a2.zvalue);
}

vec2d operator-(const vec2d &a1, const vec2d &a2){
	return vec2d( a1.xvalue - a2.xvalue, a1.zvalue - a2.zvalue);
}

vec2d operator-(const vec2d &a){
	return vec2d( -a.xvalue, -a.zvalue);
}

vec2d operator/(const vec2d &a, double b){
	return vec2d( a.xvalue/b, a.zvalue/b);
}

vec2d operator*(const vec2d &v, double a){
	return vec2d( a*v.xvalue, a*v.zvalue);
}

vec2d operator*(double a, const vec2d &v){
	return vec2d( a*v.xvalue, a*v.zvalue);                                                                                                                                                                                                                                             
}

//double dot(vec2d &a1, vec2d &a2){
//}

double angle(const vec2d &a1, const vec2d &a2){
	double xx, zz, phi;
	xx = a2.x() - a1.x();
	zz = a2.z() - a1.z();
	if( xx > 0){
		if ( zz >= 0) phi = atan( zz / xx );
		else  phi = 2*M_PI - atan(- zz / xx );
	}else if ( xx < 0 ) {
		if ( zz >= 0) phi = M_PI - atan( - zz / xx );
		else  phi = M_PI + atan(  zz / xx );
	} else {
		if ( zz >= 0) phi= M_PI/2;
		else  phi = 3*M_PI/2;
	}
	return phi;
}

void vec2d::periodic_range(double Lx_) {
	if (xvalue < 0 ){
		xvalue += Lx_;
	} else if (xvalue > Lx_){
		xvalue -= Lx_;
	}
}

void vec2d::periodic_range_xz(double Lx_, double Lz_pd_) {
	if (xvalue < 0 ){
		xvalue += Lx_;
	} else if (xvalue >= Lx_){
		xvalue -= Lx_;
	}
	if (zvalue < 0 ){
		zvalue += Lz_pd_;
	} else if (zvalue >= Lz_pd_){
		zvalue -= Lz_pd_;
	}
}

bool operator==(const vec2d &a1, const vec2d &a2){
	if (a1.xvalue == a2.xvalue && a1.zvalue == a2.zvalue)
		return true;
	return false;
}

bool operator!=(const vec2d &a1, const vec2d &a2){
	if (a1.xvalue != a2.xvalue || a1.zvalue != a2.zvalue )
		return true;
	return false;
}

void operator+=(vec2d &a1, const vec2d &a2){
	a1.xvalue += a2.xvalue;
	a1.zvalue += a2.zvalue;
}

vec2d n_vec(vec2d v1, vec2d v2){
	vec2d v = v1 - v2;
	return v / v.abs();
}

void vec2d::cout(){
	//	std::cout << '(' << xvalue << ',' << yvalue << ')' << std::endl;
	//	std::cout << "r 0.5 " << endl;
	std::cout << "c " << xvalue << ' ' << 0 << ' ' << zvalue << std::endl;
}

void vec2d::out_data(){
	std::cout << xvalue << ' ' << zvalue << std::endl;
}

double vec2d::abs(){
	return sqrt(xvalue*xvalue + zvalue*zvalue);
}
