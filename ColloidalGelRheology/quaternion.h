//
//  quaternion.h
//  ColloidalGelRheology
//
//  Created by Ryohei SETO on 6/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ColloidalGelRheology_quaternion_h
#define ColloidalGelRheology_quaternion_h

#include "vec3d.h"

class quaternion {
    
public:
	/* variables */
	double q[4];
    
    quaternion (void){
        q[0] = 1.; 
        q[1] = q[2] = q[3] = 0.;
    }
    
	quaternion (const double &_q0, 
                const double &_q1, 
                const double &_q2, 
                const double &_q3){
        q[0] = _q0, q[1] = _q1, q[2] = _q2 , q[3] = _q3;}
    
    quaternion (const double &_q0, 
                const vec3d &_u){
        q[0] = _q0, q[1] = _u.x, q[2] = _u.y , q[3] = _u.z;}
    
    
	~quaternion (void){}
    
    void set(const double &_q0, 
             const double &_q1, 
             const double &_q2, 
             const double &_q3){
        q[0] = _q0, q[1] = _q1, q[2] = _q2 , q[3] = _q3;}
    void set(const double &_q0, const vec3d &_v){
        q[0] = _q0, q[1] = _v.x, q[2] = _v.y , q[3] = _v.z;
    }
    
    inline quaternion& operator = (const quaternion& qq){
        for (int i = 0; i< 4; i++)
            q[i] = qq.q[i];
		return *this;
	}
    
    inline friend quaternion operator * (const quaternion &qa, const quaternion &qb){
        double pq[4];
        pq[0] = qa.q[0]*qb.q[0] - qa.q[1]*qb.q[1] - qa.q[2]*qb.q[2] - qa.q[3]*qb.q[3];
        pq[1] = qa.q[1]*qb.q[0] + qa.q[0]*qb.q[1] - qa.q[3]*qb.q[2] + qa.q[2]*qb.q[3];
        pq[2] = qa.q[2]*qb.q[0] + qa.q[3]*qb.q[1] + qa.q[0]*qb.q[2] - qa.q[1]*qb.q[3];
        pq[3] = qa.q[3]*qb.q[0] - qa.q[2]*qb.q[1] + qa.q[1]*qb.q[2] + qa.q[0]*qb.q[3];
		return quaternion(pq[0], pq[1], pq[2], pq[3]);
	}
    
    inline vec3d ori_backward( const vec3d & v){
        double c = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
        vec3d qvec(q[1],q[2],q[3]);
        return (c*v + 2*q[0]*cross(v,qvec) + 2*dot(v, qvec)*qvec);
    }
    
    inline vec3d ori_forward( const vec3d & v){
        double c = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
        vec3d qvec(q[1],q[2],q[3]);
        return (c*v - 2*q[0]*cross(v,qvec) + 2*dot(v, qvec)*qvec);
    }    
    
    
    inline quaternion cong (){  
        return quaternion(q[0], - q[1], - q[2], -q[3]);
    }
    inline vec3d vec(){
        return vec3d(q[1],q[2],q[3]);
    }
    
    inline vec3d direction_pcoodinate(const vec3d &v){
        // TO BE ADDED
        // The direction of bond is recoreded.
        // with quatanion.
        
        return vec3d(1,0,0) ;
    }
    
    inline void infinitesimalRotation(vec3d & dw){
        q[0] += 0.5*(            - dw.x * q[1] - dw.y * q[2] - dw.z * q[3]);
        q[1] += 0.5*(dw.x * q[0]               - dw.z * q[2] + dw.y * q[3]);
        q[2] += 0.5*(dw.y * q[0] + dw.z * q[1]             - dw.x * q[3]);
        q[3] += 0.5*(dw.z * q[0] - dw.y * q[1] + dw.x * q[2]            );
    }
    
    inline void normalize(){
        double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
        for(int i=0; i < 4; i++)
            q[i] = q[i] / norm;
        
    }
    
    inline void rotmatrix(double **rot){
        double q12_2 = 2*q[1]*q[2];
        double q03_2 = 2*q[0]*q[3];
        double q13_2 = 2*q[1]*q[3];
        double q02_2 = 2*q[0]*q[2];
        double q23_2 = 2*q[2]*q[3];
        double q01_2 = 2*q[0]*q[1];
        double q0qq = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] - q[3]*q[3]; 
        rot[0][0] = 2*q[1]*q[1] + q0qq;
        rot[0][1] = q12_2   - q03_2; // 2*(q1*q2 - q0*q3);
        rot[0][2] = q13_2   - q02_2; // 2*(q1*q3 + q0*q2);
        
        rot[1][0] = q12_2   + q03_2; // 2*(q1*q2 + q0*q3);
        rot[1][1] = 2*q[2]*q[2] + q0qq;
        rot[1][2] = q23_2   - q01_2; // 2*(q2*q3 - q0*q1);
        
        rot[2][0] = q13_2   - q02_2; //2*(q1*q3 - q0*q2);
        rot[2][1] = q23_2   + q01_2; //2*(q2*q3 + q0*q1);
        rot[2][2] = 2*q[3]*q[3] + q0qq;
    }
    
};


#endif
