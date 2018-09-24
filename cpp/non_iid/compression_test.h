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
        return G(p, v, d, num_blocks) + (alph_size-1) * G(q, v, d, num_blocks);
}

// Section 6.3.4 - Compression Estimate
// data is assumed to be binary (e.g., bit string)
double compression_test(byte* data, long len){
	int j, d, b = 6;
	long i, num_blocks, v;
	unsigned int block, alph_size = 1 << b; 
	unsigned int dict[alph_size];
	double X, sigma, p, p_lo, p_hi, eps, exp;
	double ldomain, hdomain, lbound, hbound, lvalue, hvalue, pVal, lastP;

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
	X -= ZALPHA * sigma/sqrt(v);

	ldomain = 1.0 / alph_size;
	hdomain = 1.0;

        lbound = ldomain;
        hbound = hdomain;

        lvalue = DBL_INFINITY;
        hvalue = -DBL_INFINITY;

        //Note that the bounds are in [0,1], so overflows aren't an issue
        //But underflows are.
        p = (lbound + hbound) / 2.0;
	pVal = com_exp(p, (1.0-p)/(alph_size-1), alph_size, v, d, num_blocks);

        //We don't need the initial pVal invariant, as our initial bounds are infinite.
        //We don't need the initial bounds, as they are set to the domain bounds
        for(j=0; j<ITERMAX; j++) {
                //Have we reached "equality"?
                if(relEpsilonEqual(pVal, X, ABSEPSILON, RELEPSILON, 4)) break;

                //Now update based on the found pVal
                if(X < pVal) {
                        lbound = p;
                        lvalue = pVal;
                } else {
                        hbound = p;
                        hvalue = pVal;
                }

                //We now verify that ldomain <= lbound < p < hbound <= hdomain
                //and that target in [ lvalue, hvalue ]
                if(lbound >= hbound) {
                        p = fmin(fmax(lbound, hbound),hdomain);
                        break;
                }

                //invariant. If this isn't true, then we can't evaluate here.
                if(!(INCLOSEDINTERVAL(lbound, ldomain, hdomain) && INCLOSEDINTERVAL(hbound,  ldomain, hdomain))) {
                        p = hdomain;
                        break;
                }

                //invariant. If this isn't true, then seeking the value within this interval doesn't make sense.
                if(!INCLOSEDINTERVAL(X, lvalue, hvalue)) {
                        p = hdomain;
                        break;
                }

                //Update p
                lastP = p;
                p = (lbound + hbound) / 2.0;

                //invariant. If this isn't true, then further calculation isn't really meaningful.
                if(!INOPENINTERVAL(p,  lbound, hbound)) {
                        p = hbound;
                        break;
                }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
                //Look for a cycle
                if(lastP == p) {
                        p = hbound;
                        break;
                }
#pragma GCC diagnostic pop

		pVal = com_exp(p, (1.0-p)/(alph_size-1), alph_size, v, d, num_blocks);

                //invariant: If this isn't true, then this isn't loosly monotonic
                if(!INCLOSEDINTERVAL(pVal, lvalue, hvalue)) {
                        p = hbound;
                        break;
                }
        }//for loop

        return -log2(p)/b;
}
