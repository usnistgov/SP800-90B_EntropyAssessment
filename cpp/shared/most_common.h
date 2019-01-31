#pragma once
#include "../shared/utils.h"

// Section 6.3.1 - Most Common Value Estimate
double most_common(byte* data, const long len, const int alph_size){

	long counts[alph_size];
	long i, mode;
	double pmax, ubound;

	for(i = 0; i < alph_size; i++) counts[i] = 0;
	for (i = 0; i < len; i++) counts[data[i]]++;

	mode = 0;
	for(i = 0; i < alph_size; i++){
		if(counts[i] > mode) mode = counts[i];
	}

	pmax = mode/(double)len;
	ubound = pmax + ZALPHA*sqrt(pmax*(1.0-pmax)/(len-1.0));

	return -log2(min(1.0, ubound));
}
