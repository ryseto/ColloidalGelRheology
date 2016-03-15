//
//  tools.h
//  ColloidalGelRheology
//
//  Created by Ryohei Seto on 2012/08/17.
//
//

#ifndef ColloidalGelRheology_tools_h
#define ColloidalGelRheology_tools_h
int ipow(int p, int q){
	int x = 1;
	for (int i=0; i < q; i++){
		x *= p;
	}
	return x;
}
#endif
