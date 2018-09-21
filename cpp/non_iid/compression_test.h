#pragma once
#include "../shared/utils.h"

double G(double z, long v, int d, long num_blocks){
	long u, t;
	double ret, inner_sum;

	// precompute inner sum
	inner_sum = 0.0;
	for(u = 2; u <= d; u++) inner_sum += log2(u) * z*z*pow(1.0-z, u-1);
	
	// compute full sum
	ret = 0.0;
	for(t = d+1; t <= num_blocks; t++){
		ret += log2(t) * z*pow(1.0-z, t-1) + inner_sum;
		inner_sum += log2(t) * z*z*pow(1.0-z, t-1);
	}
	
	return ret/v;
}

double com_exp(double p, double q, unsigned int alph_size, long v, int d, long num_blocks){
        return G(p, v, d, num_blocks) + alph_size * G(q, v, d, num_blocks);
}

// Section 6.3.4 - Compression Estimate
// data is assumed to be binary (e.g., bit string)
double compression_test(byte* data, long len){
	int j, d, b = 6;
	long i, num_blocks, v;
	unsigned int block, alph_size = 1 << b; 
	unsigned int dict[alph_size];
	double X, sigma, p, p_lo, p_hi, eps, exp;

	d = 1000;
	num_blocks = len/b;
	X = 0;
	sigma = 0;

	if(num_blocks <= d){
		printf("\t*** Warning: not enough samples to run compression test (need more than %d) ***\n", d);
		return -1.0;
	}

	// create dictionary
	for(i = 0; i < alph_size; i++) dict[i] = 0;
	for(i = 0; i < d; i++){
		block = 0;
		for(j = 0; j < b; j++) block |= (data[i*b + j] & 0x1) << (b-j-1);
		dict[block] = i+1;
	}

	// test data against dictionary
	v = num_blocks - d;
	for(i = d; i < num_blocks; i++){
		block = 0;
		for(j = 0; j < b; j++) block |= (data[i*b + j] & 0x1) << (b-j-1);
		X += log2(i+1-dict[block]);
		sigma += log2(i+1-dict[block])*log2(i+1-dict[block]);
		dict[block] = i+1;
	}

	// compute mean and stdev
	X /= v;
	sigma = 0.5907 * sqrt(sigma/(v-1.0) - X*X);

        // binary search for p
	X -= 2.576 * sigma/sqrt(v);
        eps = 1.0 / (1 << 20); // 2^-20
        p_lo = 1.0 / alph_size;
        p_hi = 1.0 - eps; // avoid division by zero
        do{
                p = (p_lo + p_hi) / 2.0;
                exp = com_exp(p, (1.0-p)/alph_size, alph_size, v, d, num_blocks);

                if(X < exp) p_lo = p;
                else p_hi = p;

                if((com_exp(p_lo, (1.0-p_lo)/alph_size, alph_size, v, d, num_blocks) < X) || (com_exp(p_hi, (1.0-p_hi)/alph_size, alph_size, v, d, num_blocks) > X)){
                        // binary search failed, settle for one bit of entropy
			printf("\t *** WARNING: binary search for compression test failed ***\n");
                        p = 1.0/alph_size; 
                        break;
                }
        }while(fabs(p_hi - p_lo) > eps);

        return -log2(p)/b;
}
