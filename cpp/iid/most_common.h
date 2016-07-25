#include "../utils.h"

double most_common(byte data[]){

	map<byte, int> p;
	map_init(p);

	for(int i = 0; i < SIZE; i++){
		p[data[i]]++;
	}

	int mode = 0;
	map<byte, int>::iterator itr;
	for(itr = p.begin(); itr != p.end(); ++itr){
		if(itr->second > mode){
			mode = itr->second;
		}
	}

	double pmax = mode / (double)SIZE;
	double ubound = pmax + 2.576 * sqrt(pmax * (1.0 - pmax) / (double)SIZE);
	double pu = min(1.0, ubound);

	return -log2(pu);
}