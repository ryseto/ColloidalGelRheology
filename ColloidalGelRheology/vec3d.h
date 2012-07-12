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
	inline vec3d (void){}
	inline vec3d (const double &_x, 
                  const double &_y, 
                  const double &_z): x(_x), y(_y), z(_z) {}
	inline ~vec3d(void){}
    
	/* operators */
	inline vec3d& operator = (const vec3d& v){
		x = v.x, y = v.y, z = v.z;
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
		return vec3d( a1.x + a2.x, a1.y + a2.y, a1.z + a2.z);
        
	}
	inline friend vec3d operator + (const vec3d &v){
		return v;
	}
	
	/* subtraction */
	inline friend vec3d operator - (const vec3d &a1, const vec3d &a2){
		return vec3d( a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
	}
	inline friend vec3d operator - (const vec3d &v){
		return vec3d( -v.x, -v.y, -v.z);
	}
	
	/* multiplication */
	inline friend vec3d operator * (const double &d, const vec3d &v){
		return vec3d( d*v.x, d*v.y, d*v.z);
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
		return a1.x*a2.x + a1.y*a2.y + a1.z*a2.z;
	}
	inline friend double dot_2d(const vec3d &a1, const vec3d &a2){
		return a1.x*a2.x + a1.z*a2.z;
	}
	
	
	/* vector product */
	inline friend vec3d cross(const vec3d &v1, const vec3d &v2){
		return vec3d(v1.y*v2.z - v1.z*v2.y, 
					 v1.z*v2.x - v1.x*v2.z, 
					 v1.x*v2.y - v1.y*v2.x);
		
	}	
	inline friend vec3d cross_2d(const vec3d &v1, const vec3d &v2){
		return vec3d(0, v1.z*v2.x - v1.x*v2.z, 0 );
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
		x += v.x, y += v.y, z += v.z;
		return *this;
	}
	inline vec3d& operator -=(const vec3d &v){
		x -= v.x, y -= v.y, z -= v.z;
		return *this;
	}
	inline vec3d& operator *=(const double & d){
        
		x *= d, y *= d, z *= d;
		return *this;
	}
	inline vec3d& operator *=(const int & i){
		double d_tmp = (double)i;
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
		return *this;
	}
	inline vec3d& operator /=(const double & d){
		double d_tmp = 1.0/d;
        
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
		return 	*this;
	}
	inline vec3d& operator /=(const int & i){
		double d_tmp = 1.0/i;
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
		return *this;
	}
    
	/* utility */
	inline void set(const double &_x, const double &_y, const double &_z){x=_x, y=_y, z=_z;}
	inline void reset(){ x = y = z = 0.;}
    
    inline bool eq_zero(){
        if ( x==0 && y==0 && z==0 )
            return true;
        return false;
    }
    
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
		return x*x + y*y + z*z; 
	}
	inline double sq_norm_xy(){return x*x + y*y;}
	inline double sq_norm_xz(){return x*x + z*z;}
	
    inline double norm(){
        return sqrt(x*x + y*y + z*z);
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
		double dy = a1.z-a2.z;
		return dx*dx+ dy*dy;
	}
	inline void rotateInfinitesimal(const vec3d &dphi){
		/* dphi must be small vector. */
#ifndef TWODIMENSION
		//3D
		(*this) += cross(dphi, *this);
#else
		//2D 
		(*this) += vec3d(dphi.y*z,0,- dphi.y*x);
		//(*this) += cross_2d(dphi, *this);
#endif
	}
	inline void cerr(){ std::cerr << x << ' ' << y << ' ' << z << std::endl;}
};


#endif
