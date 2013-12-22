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
#ifdef TWODIMENSION
	// 2D
	inline vec3d (const double &_x,
				  const double &_y,
				  const double &_z): x(_x), y(0), z(_z) {}
	inline vec3d (const double &_x,
				  const double &_z): x(_x), z(_z) {}
#else
	// 3D
	inline vec3d (const double &_x,
				  const double &_y,
				  const double &_z): x(_x), y(_y), z(_z) {}
#endif
	inline ~vec3d(void){}
	
	/* operators */
	inline vec3d& operator = (const vec3d& v){
#ifdef TWODIMENSION
		// 2D
		x = v.x, z = v.z;
#else
		// 3D
		x = v.x, y = v.y, z = v.z;
#endif
		return *this;
	}
	
	inline friend bool operator==(const vec3d &v1, const vec3d &v2){
#ifdef TWODIMENSION
		// 2D
		if (v1.x == v2.x && v1.z == v2.z)
			return true;
#else
		// 3D
		if (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z)
			return true;
#endif
		return false;
	}
	inline friend bool operator!=(const vec3d &v1, const vec3d &v2){
#ifdef TWODIMENSION
		// 2D
		if (v1.x != v2.x || v1.z != v2.z )
				return true;
#else
		// 3D
		if (v1.x != v2.x || v1.y != v2.y || v1.z != v2.z )
			return true;
#endif
		return false;
	}
	
	inline friend vec3d operator + (const vec3d &a1, const vec3d &a2){
#ifdef TWODIMENSION
		// 2D
		return vec3d(a1.x+a2.x, a1.z+a2.z);
#else
		// 3D
		return vec3d(a1.x+a2.x, a1.y+a2.y, a1.z+a2.z);
#endif
	}

	inline friend vec3d operator + (const vec3d &v){
		return v;
	}
	
	/* subtraction */
	inline friend vec3d operator - (const vec3d &a1, const vec3d &a2){
#ifdef TWODIMENSION
		// 2D
		return vec3d(a1.x-a2.x, a1.z-a2.z);
#else
		// 3D
		return vec3d(a1.x-a2.x, a1.y-a2.y, a1.z-a2.z);
#endif
	}
	inline friend vec3d operator - (const vec3d &v){
#ifdef TWODIMENSION
		// 2D
		return vec3d(-v.x, -v.z);
#else
		// 3D
		return vec3d(-v.x, -v.y, -v.z);
#endif
		
	}
	
	/* multiplication */
	inline friend vec3d operator * (const double &d, const vec3d &v){
#ifdef TWODIMENSION
		// 2D
		return vec3d(d*v.x, d*v.z);
	
#else
		// 3D
		return vec3d(d*v.x, d*v.y, d*v.z);
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
#ifdef TWODIMENSION
		// 2D
		return a1.x*a2.x+a1.z*a2.z;
#else
		// 3D
		return a1.x*a2.x+a1.y*a2.y+a1.z*a2.z;
#endif
	}
	inline friend double dot_2d(const vec3d &a1, const vec3d &a2){
		return a1.x*a2.x + a1.z*a2.z;
	}
	
	
	/* vector product */
	inline friend vec3d cross(const vec3d &v1, const vec3d &v2){
#ifdef TWODIMENSION
		// 2D
		return vec3d(0, v1.z*v2.x - v1.x*v2.z, 0);
#else
		// 3D
		return vec3d(v1.y*v2.z-v1.z*v2.y,
					 v1.z*v2.x-v1.x*v2.z,
					 v1.x*v2.y-v1.y*v2.x);
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
#ifdef TWODIMENSION
		// 2D
		x += v.x, z += v.z;
#else
		// 3D
		x += v.x, y += v.y, z += v.z;
#endif
		return *this;
	}
	inline vec3d& operator -=(const vec3d &v){
#ifdef TWODIMENSION
		// 2D
		x -= v.x, z -= v.z;
#else
		// 3D
		x -= v.x, y -= v.y, z -= v.z;
#endif
		return *this;
	}
	inline vec3d& operator *=(const double & d){
#ifdef TWODIMENSION
		// 2D
		x *= d, z *= d;
#else
		// 3D
		x *= d, y *= d, z *= d;
#endif
		return *this;
	}
	inline vec3d& operator *=(const int & i){
		double d_tmp = (double)i;
#ifdef TWODIMENSION
		// 2D
		x *= d_tmp, z *= d_tmp;
#else
		// 3D
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
#endif
		return *this;
	}
	inline vec3d& operator /=(const double & d){
		double d_tmp = 1.0/d;
#ifdef TWODIMENSION
		// 2D
		x *= d_tmp, z *= d_tmp;
#else
		// 3D
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
#endif
		return 	*this;
	}
	inline vec3d& operator /=(const int & i){
		double d_tmp = 1.0/i;
#ifdef TWODIMENSION
		// 2D
		x *= d_tmp, z *= d_tmp;
#else
		// 3D
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
#endif
		return *this;
	}
	
	/* utility */

	inline void set(const double &_x, const double &_y, const double &_z){x=_x, y=_y, z=_z;}
	inline void reset(){ x = y = z = 0.;}

#ifdef TWODIMENSION
	// 2D
	inline vec3d unitvector(){
		return (*this) / norm_2d();
	}
#else
	//3D
	inline vec3d unitvector(){
		return (*this) / norm();
	}
#endif
	
	inline void sign_reverse(){ (*this) = -(*this);}
	inline double sq_norm(){
#ifdef TWODIMENSION
		// 2D
		return x*x + z*z;
#else
		// 3D
		return x*x + y*y + z*z;
#endif
	}
	inline double sq_norm_xy(){return x*x + y*y;}
	inline double sq_norm_xz(){return x*x + z*z;}
	
	inline double norm(){
#ifdef TWODIMENSION
		// 2D
		return sqrt(x*x + z*z);
#else
		// 3D
		return sqrt(x*x + y*y + z*z);
#endif
	}
	
	inline double norm_2d(){return sqrt(x*x+z*z);}
	
	inline vec3d division_2d(const double &d){
		double d_tmp = 1.0/d;
		return vec3d( x*d_tmp, 0, z*d_tmp);
	}
	
	
	inline friend double dist(const vec3d &a1, const vec3d &a2){
		return (a1 - a2).norm();
	}
	inline friend double sq_dist(const vec3d &a1, const vec3d &a2){
		return (a1 - a2).sq_norm() ;
	}
	inline friend double sq_dist_2d(const vec3d &a1, const vec3d &a2){
		double dx = a1.x-a2.x;
		double dz = a1.z-a2.z;
		return dx*dx + dz*dz;
	}
	
#ifdef TWODIMENSION
	// 2D
	inline void rotateInfinitesimal(const double &dphi_y){
		/* dphi must be small vector. */
		(*this) += vec3d( dphi_y * z, 0 , - dphi_y * x);
	}
#else
	// 3D
	inline void rotateInfinitesimal(const vec3d &dphi){
		/* dphi must be small vector. */
		(*this) += cross(dphi, *this);
	}
#endif
	
	inline void cerr(){ std::cerr << x << ' ' << y << ' ' << z << std::endl;}
};


#endif




