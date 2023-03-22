#pragma once
#include "../shared/utils.h"

inline void kahan_add(double &sum, double &comp, double in){
	double y, t; 

	y = in - comp;
	t = sum + y;
	comp = (t - sum) - y;
	sum = t;
}

//There is some cleverness associated with this calculation of G; in particular,
//one doesn't need to calculate all the terms independently (they are inter-related!)
//See UL's implementation comments here: https://bit.ly/UL90BCOM 
//Look in the section "Compression Estimate G Function Calculation"
double G(double z, int d, long num_blocks){
	double Ai=0.0, Ai_comp=0.0;
	double firstSum=0.0, firstSum_comp=0.0;
	long v = num_blocks - d;
	double Ad1;

	long double Bi;
	long double Bterm;
	long double ai;
	long double aiScaled;
	bool underflowTruncate;

	assert(d>0);
	assert(num_blocks>d);

	//i=2
	Bterm = (1.0L-(long double)z);
	//Note: B_1 isn't needed, as a_1 = 0
	//B_2
	Bi = Bterm;

	//Calculate A_{d+1}
	for(int i=2; i<=d; i++) {
		//calculate the a_i term
		kahan_add(Ai, Ai_comp, log2l((long double)i)*Bi);

		//Calculate B_{i+1}
		Bi *= Bterm;
	}

	//Store A_{d+1}
	Ad1 = Ai;

	underflowTruncate = false;
	//Now calculate A_{num_blocks} and the sum of sums term (firstsum)
	for(long i=d+1; i<=num_blocks-1; i++) {
		//calculate the a_i term
		ai = log2l((long double)i)*Bi;

		//Calculate A_{i+1}
		kahan_add(Ai, Ai_comp, (double)ai);
		//Sum in A_{i+1} into the firstSum

		//Calculate the tail of the sum of sums term (firstsum)
		aiScaled = (long double)(num_blocks-i) * ai;
		if((double)aiScaled > 0.0) {
			kahan_add(firstSum, firstSum_comp, (double)aiScaled);
		} else {
			underflowTruncate = true;
			break;
		}

		//Calculate B_{i+1}
		Bi *= Bterm;
	}

	//Ai now contains A_{num_blocks} and firstsum contains the tail
	//finalize the calculation of firstsum
	kahan_add(firstSum, firstSum_comp, ((double)(num_blocks-d))*Ad1);

	//Calculate A_{num_blocks+1}
	if(!underflowTruncate) {
		ai = log2l((long double)num_blocks)*Bi;
		kahan_add(Ai, Ai_comp, (double)ai);
	}

	return 1/(double)v * z*(z*firstSum + (Ai - Ad1));
}

double com_exp(double p, unsigned int alph_size, int d, long num_blocks){
	double q = (1.0-p)/((double)alph_size-1.0);
        return G(p, d, num_blocks) + ((double)alph_size-1.0) * G(q, d, num_blocks);
}

// Section 6.3.4 - Compression Estimate
// data is assumed to be binary (e.g., bit string)
double compression_test(uint8_t* data, long len, const int verbose, const char *label){
	int j, d, b = 6;
	long i, num_blocks, v;
	unsigned int block, alph_size = 1 << b; 
	unsigned int dict[alph_size];
	double X=0.0, X_comp=0.0;
	double sigma=0.0, sigma_comp=0.0;
	double p, entEst;
	double ldomain, hdomain, lbound, hbound, lvalue, hvalue, pVal, lastP;

	d = 1000;
	num_blocks = len/b;

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
		kahan_add(X, X_comp, log2(i+1-dict[block]));
		kahan_add(sigma, sigma_comp, log2(i+1-dict[block])*log2(i+1-dict[block]));
		dict[block] = i+1;
	}

	// compute mean and stdev
	X /= v;
	sigma = 0.5907 * sqrt(sigma/(v-1.0) - X*X);

	if(verbose == 2) {
		printf("%s Compression Estimate: X-bar = %.17g, ", label, X);
		printf("sigma-hat = %.17g, ", sigma);
	} else if(verbose == 3) {
		printf("%s Compression Estimate: X-bar = %.17g\n", label, X);
		printf("%s Compression Estimate: sigma-hat = %.17g\n", label, sigma);
	}

        // binary search for p
	X -= ZALPHA * sigma/sqrt(v);

	if(verbose == 3) printf("%s Compression Estimate: X-bar' = %.17g\n", label, X);

	if(com_exp(1.0/(double)alph_size, alph_size, d, num_blocks) > X) {
		ldomain = 1.0 / (double)alph_size;
		hdomain = 1.0;

		lbound = ldomain;
		hbound = hdomain;

		lvalue = DBL_INFINITY;
		hvalue = -DBL_INFINITY;

		//Note that the bounds are in [0,1], so overflows aren't an issue
		//But underflows are.
		p = (lbound + hbound) / 2.0;
		pVal = com_exp(p, alph_size, d, num_blocks);

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
				//This is a search failure. We need to return "full entropy"  (as directed in step #8).
				p = ldomain;
				break;
			}

			//invariant. If this isn't true, then seeking the value within this interval doesn't make sense.
			if(!INCLOSEDINTERVAL(X, lvalue, hvalue)) {
				//This is a search failure. We need to return "full entropy"  (as directed in step #8).
				p = ldomain;
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

			pVal = com_exp(p, alph_size, d, num_blocks);

			//invariant: If this isn't true, then this isn't loosely monotonic
			if(!INCLOSEDINTERVAL(pVal, lvalue, hvalue)) {
				p = hbound;
				break;
			}
		}//for loop
	} else {
		p = -1.0;
	}

	if(p > 1.0 / (double)alph_size) {
        	entEst = -log2(p)/b;
	
		if(verbose == 3) printf("%s Compression Estimate: Found p.\n", label);
	} else {
		p = 1.0 / (double)alph_size;
		entEst = 1.0;
		if(verbose == 3) printf("%s Compression Estimate: Could Not Find p. Proceeding with the lower bound for p.\n", label);
	}

	if(verbose == 2) printf("p = %.17g\n", p);
	else if(verbose == 3) {
		printf("%s Compression Estimate: p = %.17g\n", label, p);
		printf("%s Compression Estimate: min entropy = %.17g\n", label, entEst);
	}

        return entEst;
}
