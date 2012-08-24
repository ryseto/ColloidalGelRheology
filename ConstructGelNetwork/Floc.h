//
//  Floc.h
//  ColloidalGelRheology
//
//  Created by Ryohei Seto on 2012/08/20.
//
//

#ifndef ColloidalGelRheology_Floc_h
#define ColloidalGelRheology_Floc_h
#include "vec2d.h"
#include <vector>

class Floc{
private:
	unsigned long floc_size;
	vector<vec2d> positions;
	
	double Lx;
	double Ly;

	
protected:

public:
	Floc(){
		p_center.set(0, 0);
	}
	~Floc(){
	}
	
	vec2d p_center;
	int parent;
	void setSystemSize(double Lx_, double Ly_){
		Lx = Lx_;
		Ly = Ly_;
	}

	void setPosition(vector<vec2d> &pos);
	void setCenterPosition(vec2d p_center_){
		p_center = p_center_;
	}
	vector<int> connect;
	void move(vec2d translation){
		//		for (int i=0; i < floc_size; i++){
		//			positions[i] += translation;
		//		}
		p_center += translation;
		if (p_center.x() < 0){
			p_center.shift( Lx, 0);
		}else if (p_center.x() > Lx){
			p_center.shift(-Lx, 0);
		}
		if (p_center.y() < 0){
			p_center.shift(0, Ly);
		} else if (p_center.y() > Ly){
			p_center.shift(0, -Ly);
		}
	}
	
	unsigned long size(){
		return floc_size;
	}
	vec2d pos(int i){
		vec2d p = p_center + positions[i];

		if (p.x() < 0){
			p.shift(Lx, 0);
		}else if (p.x() > Lx){
			p.shift(-Lx, 0);
		}
		if (p.y() < 0){
			p.shift(0, Ly);
		}else if (p.y() > Ly){
			p.shift(0, -Ly);
		}
		return p;
	}
};
#endif
