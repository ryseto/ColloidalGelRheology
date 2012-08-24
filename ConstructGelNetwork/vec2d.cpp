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

double vec2d::x(double L){
	if (xvalue < 0)
		return xvalue + L;
	else if (xvalue > L)
		return xvalue - L;
	else
		return xvalue;
}

vec2d operator+(vec2d a1, vec2d a2){
	return vec2d( a1.xvalue + a2.xvalue, a1.yvalue + a2.yvalue);
}

vec2d operator-(vec2d a1, vec2d a2){
	return vec2d( a1.xvalue - a2.xvalue, a1.yvalue - a2.yvalue);
}

vec2d operator-(vec2d a){
	return vec2d( -a.xvalue, -a.yvalue);
}

vec2d operator/(vec2d a, double b){
	return vec2d( a.xvalue/b, a.yvalue/b);
}

vec2d operator*(vec2d v, double a){
	return vec2d( a*v.xvalue, a*v.yvalue);
}

vec2d operator*(double a, vec2d v){
	return vec2d( a*v.xvalue, a*v.yvalue);
}

double dot(vec2d a1, vec2d a2){
	return a1.xvalue*a2.xvalue + a1.yvalue*a2.yvalue;
}

double angle(vec2d a1, vec2d a2){
	// ベクトル a1 → a2 の角
	double xx, yy, phi;
	xx = a2.x() - a1.x();
	yy = a2.y() - a1.y();
	if( xx > 0){
		if ( yy >= 0) phi = atan( yy / xx );
		else  phi = 2*M_PI - atan(- yy / xx );
	}else if ( xx < 0 ) {
		if ( yy >= 0) phi = M_PI - atan( - yy / xx );
		else  phi = M_PI + atan(  yy / xx );
	} else {
		if ( yy >= 0) phi= M_PI/2;
		else  phi = 3*M_PI/2;
	}
	return phi;
}

void vec2d::periodic_range(double L) {
	if (xvalue < 0 ){
		xvalue += L;
	} else if (xvalue > L){
		xvalue -= L;
	}
}

void vec2d::periodic_range_xy(double Lx, double Ly) {
	if (xvalue < 0 ){
		xvalue += Lx;
	} else if (xvalue >= Lx){
		xvalue -= Lx;
	}
	if (yvalue < 0 ){
		yvalue += Ly;
	} else if (yvalue >= Ly){
		yvalue -= Ly;
	}
}

bool operator==(vec2d a1, vec2d a2){
	if (a1.xvalue == a2.xvalue && a1.yvalue == a2.yvalue)
		return true;
	return false;
}

bool operator!=(vec2d a1, vec2d a2){
	if (a1.xvalue != a2.xvalue || a1.yvalue != a2.yvalue )
		return true;
	return false;
}

void operator+=(vec2d &a1, vec2d a2){
	a1.xvalue += a2.xvalue;
	a1.yvalue += a2.yvalue;
}

vec2d n_vec(vec2d v1, vec2d v2){
	vec2d v = v1 - v2;
	return v / v.abs();
}

void vec2d::cout(){
	//	std::cout << '(' << xvalue << ',' << yvalue << ')' << std::endl;
	//	std::cout << "r 0.5 " << endl;
	std::cout << "c " << xvalue << ' '
	<< 0 << ' ' << yvalue << std::endl;
}

void vec2d::out_data(){
	std::cout << xvalue << ' ' << yvalue << std::endl;
}

double vec2d::abs(){
	return sqrt(xvalue*xvalue + yvalue*yvalue);
}

double dist(vec2d a1, vec2d a2){
	return sqrt(sq(a1.xvalue-a2.xvalue) + sq(a1.yvalue-a2.yvalue));
}

double sq_dist(vec2d a1, vec2d a2){
	return sq(a1.xvalue-a2.xvalue) + sq(a1.yvalue-a2.yvalue);
}

double sq_dist_pd(vec2d a1, vec2d a2){
	double xx = a1.xvalue - a2.xvalue;
	if ( xx > L_SIZE )
		xx -= L;
	else if ( xx < -L_SIZE )
		xx += L;
	return sq(xx) + sq(a1.yvalue-a2.yvalue);
}

