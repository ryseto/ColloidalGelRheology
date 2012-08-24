//
//  vec3d.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_vec3d_h
#define ColloidalGelRheology_vec3d_h

#include <iostream>
#include <iomanip>
#include <cmath>
class vec3d {
public:
	/* variables */
	double x;
	double y;
	double z;
	/* constructor/destructor */

	inline vec3d (void): x(0), y(0), z(0){}
#ifndef TWODIMENSION
	inline vec3d (const double &_x,
				  const double &_y,
				  const double &_z): x(_x), y(_y), z(_z) {}
#else
	inline vec3d (const double &_x,
				  const double &_y,
				  const double &_z): x(_x), y(0), z(_z) {}
#endif
	inline ~vec3d(void){}
	
	/* operators */
	inline vec3d& operator = (const vec3d& v){
#ifndef TWODIMENSION
		x = v.x, y = v.y, z = v.z;
#else
		x = v.x, z = v.z;
#endif
		return *this;
	}
	
	inline friend bool operator==(const vec3d &v1, const vec3d &v2){
#ifndef TWODIMENSION
		if (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z)
			return true;
#else
		if (v1.x == v2.x && v1.z == v2.z)
			return true;
#endif
		return false;
	}
	inline friend bool operator!=(const vec3d &v1, const vec3d &v2){
#ifndef TWODIMENSION
		if (v1.x != v2.x || v1.y != v2.y || v1.z != v2.z )
			return true;
#else
		if (v1.x != v2.x || v1.z != v2.z )
			return true;
#endif
		return false;
	}
	
	inline friend vec3d operator + (const vec3d &a1, const vec3d &a2){
#ifndef TWODIMENSION
		return vec3d( a1.x + a2.x, a1.y + a2.y, a1.z + a2.z);
#else
		return vec3d( a1.x + a2.x, 0, a1.z + a2.z);
#endif
		
	}
	inline friend vec3d operator + (const vec3d &v){
		return v;
	}
	
	/* subtraction */
	inline friend vec3d operator - (const vec3d &a1, const vec3d &a2){
#ifndef TWODIMENSION
		return vec3d( a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
#else
		return vec3d( a1.x - a2.x, 0 , a1.z - a2.z);
#endif
	}
	inline friend vec3d operator - (const vec3d &v){
#ifndef TWODIMENSION
		return vec3d( -v.x, -v.y, -v.z);
#else
		return vec3d( -v.x, 0 , -v.z);
#endif
		
	}
	
	/* multiplication */
	inline friend vec3d operator * (const double &d, const vec3d &v){
#ifndef TWODIMENSION
		return vec3d( d*v.x, d*v.y, d*v.z);
#else
		return vec3d( d*v.x, 0, d*v.z);
#endif
	}
	inline friend vec3d operator * (const vec3d &v, const double &d){
		return d*v;
	}
	inline friend vec3d operator * (const int &i, const vec3d &v){
		double d_tmp = (double)i;
		return d_tmp*v;
	}
	inline friend vec3d operator * (const vec3d &v, const int &i){
		double d_tmp = (double)i;
		return d_tmp*v;
	}
	
	/* scalar product */
	inline friend double dot(const vec3d &a1, const vec3d &a2){
#ifndef TWODIMENSION
		return a1.x*a2.x + a1.y*a2.y + a1.z*a2.z;
#else
		return a1.x*a2.x + a1.z*a2.z;
#endif
	}
	inline friend double dot_2d(const vec3d &a1, const vec3d &a2){
		return a1.x*a2.x + a1.z*a2.z;
	}
	
	
	/* vector product */
	inline friend vec3d cross(const vec3d &v1, const vec3d &v2){
#ifndef TWODIMENSION
		return vec3d(v1.y*v2.z - v1.z*v2.y,
					 v1.z*v2.x - v1.x*v2.z,
					 v1.x*v2.y - v1.y*v2.x);
#else
		return vec3d(0, v1.z*v2.x - v1.x*v2.z, 0);
#endif
		
	}
	inline friend double cross_2d(const vec3d &v1, const vec3d &v2){
		return v1.z*v2.x - v1.x*v2.z;
	}
	
	/* division */
	inline friend vec3d operator / (const vec3d &v, const double &d){
		double d_tmp = 1.0/d;
		return d_tmp*v;
	}
	inline friend vec3d operator / (const vec3d &v, const int &i){
		double d_tmp = 1.0/i;
		return d_tmp*v;
	}
	
	// assign operator
	inline vec3d& operator +=(const vec3d &v){
#ifndef TWODIMENSION
		x += v.x, y += v.y, z += v.z;
#else
		x += v.x, z += v.z;
#endif
		return *this;
	}
	inline vec3d& operator -=(const vec3d &v){
#ifndef TWODIMENSION
		x -= v.x, y -= v.y, z -= v.z;
#else
		x -= v.x, z -= v.z;
#endif
		return *this;
	}
	inline vec3d& operator *=(const double & d){
#ifndef TWODIMENSION
		x *= d, y *= d, z *= d;
#else
		x *= d, z *= d;
#endif
		return *this;
	}
	inline vec3d& operator *=(const int & i){
		double d_tmp = (double)i;
#ifndef TWODIMENSION
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
#else
		x *= d_tmp, z *= d_tmp;
#endif
		return *this;
	}
	inline vec3d& operator /=(const double & d){
		double d_tmp = 1.0/d;
#ifndef TWODIMENSION
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
#else
		x *= d_tmp, z *= d_tmp;
#endif
		return 	*this;
	}
	inline vec3d& operator /=(const int & i){
		double d_tmp = 1.0/i;
#ifndef TWODIMENSION
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
#else
		x *= d_tmp, z *= d_tmp;
#endif
		return *this;
	}
	
	/* utility */

	inline void set(const double &_x, const double &_y, const double &_z){x=_x, y=_y, z=_z;}
	inline void reset(){ x = y = z = 0.;}

//	inline bool eq_zero(){
//		if ( x==0 && y==0 && z==0 )
//			return true;
//		return false;
//	}
	
#ifndef TWODIMENSION
	//3D
	inline vec3d unitvector(){
		return (*this) / norm();
	}
#else
	inline vec3d unitvector(){
		return (*this) / norm_2d();
	}
#endif
	
	inline void sign_reverse(){ (*this) = -(*this);}
	inline double sq_norm(){
#ifndef TWODIMENSION
		return x*x + y*y + z*z;
#else
		return x*x + z*z;
#endif
	}
	inline double sq_norm_xy(){return x*x + y*y;}
	inline double sq_norm_xz(){return x*x + z*z;}
	
	inline double norm(){
#ifndef TWODIMENSION
		return sqrt(x*x + y*y + z*z);
#else
		return sqrt(x*x + z*z);
#endif
	}
	
	//#ifdef TWODIMENSION
	inline double norm_2d(){return sqrt(x*x+z*z);}
	//#endif
	
	inline vec3d division_2d(const double &d){
		double d_tmp = 1.0/d;
		return vec3d( x*d_tmp, 0, z*d_tmp);
	}
	
	
	inline friend double dist(const vec3d &a1, const vec3d &a2){
		return (a1 - a2).norm();}
	inline friend double sq_dist(const vec3d &a1, const vec3d &a2){
		return (a1 - a2).sq_norm() ;
	}
	inline friend double sq_dist_2d(const vec3d &a1, const vec3d &a2){
		double dx = a1.x-a2.x;
		double dz = a1.z-a2.z;
		return dx*dx + dz*dz;
	}
	
#ifndef TWODIMENSION
	inline void rotateInfinitesimal(const vec3d &dphi){
		/* dphi must be small vector. */
		//3D
		(*this) += cross(dphi, *this);
		//(*this) += cross_2d(dphi, *this);
	}
#else
	inline void rotateInfinitesimal(const double &dphi_y){
		/* dphi must be small vector. */
		(*this) += vec3d( dphi_y * z, 0 , - dphi_y * x);
	}
#endif
	inline void cerr(){ std::cerr << x << ' ' << y << ' ' << z << std::endl;}
};


#endif




