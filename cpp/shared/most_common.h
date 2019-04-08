#pragma once
#include "../shared/utils.h"

// Section 6.3.1 - Most Common Value Estimate
double most_common(byte* data, const long len, const int alph_size, const bool verbose){

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

	ubound = min(1.0,pmax + ZALPHA*sqrt(pmax*(1.0-pmax)/(len-1.0)));
	if(verbose) printf("MCV Estimate: mode = %ld, p-hat = %.17g, p_u = %.17g\n", mode, pmax, ubound);

	return -log2(ubound);
}
