#pragma once
#include "../shared/utils.h"

// Computed using efficient implementation in Appendix G.1.1
double F(double q){
   return q*(2.0*q*q+2.0*q+1.0);
}

double col_exp(double p){
	double q = 1.0 - p;

	return (p/(q*q))*(1.0 + 0.5*(1.0/p - 1.0/q))*F(q) - (p/q)*0.5*(1.0/p - 1.0/q);
}

// Section 6.3.2 - Collision Estimate
// data is assumed to be binary (e.g., bit string)
double collision_test(uint8_t* data, long len, const int verbose, const char *label){
	long v, i;
	int t_v;
	double X, s, p;
	double entEst;

	i = 0;
	v = 0;
	s = 0.0;

	// compute wait times until collisions
	while(i < len-1){
		if(data[i] == data[i+1]) t_v = 2; // 00 or 11
		else if(i < len-2) t_v = 3; // 101, 011, 100, or 101
		else break;
		
		v++;
		s += t_v*t_v;
		i += t_v;
	}

	// X is mean of t_v's, s is sample stdev, where
	// s^2 = (sum(t_v^2) - sum(t_v)^2/v) / (v-1)
	X = i / (double)v;
	if(verbose == 2) printf("%s Collision Estimate: X-bar = %.17g, ", label, X);
	s = sqrt((s - (i*X)) / (v-1));
	if(verbose == 2) printf("sigma-hat = %.17g, ", s);

	if(verbose == 3) {
		printf("%s Collision Estimate: v = %ld\n", label, v);
		printf("%s Collision Estimate: Sum t_i = %ld\n", label, i);
		printf("%s Collision Estimate: X-bar = %.17g\n", label, X);
		printf("%s Collision Estimate: sigma-hat = %.17g\n", label, s);
	}

	// Directly calculate p
	X -= ZALPHA * s/sqrt(v);
	//2 is the smallest meaninful value here.
	if(X < 2.0) X = 2.0;

	if(verbose == 3)
		printf("%s Collision Estimate: X-bar' = %.17g\n", label, X);

	//Uyen Dinh observed that (with the simpler F function described in UL comments) we can simplify the entire expression much further than in 90B.
	//The whole mess in 90B step 7 reduces to X'-bar = -2p^2 + 2p + 2, which we can solve using the quadratic formula.
	//We only care about the root greater than 0.5, so we only care about the "+" branch.
	//If the meanbound > 2.5, then the roots become complex, so this isn't well defined (and is processed as per the error handling specified in 90B).
	if(X < 2.5) {
		p = 0.5 + sqrt(1.25 - 0.5 * X);
		entEst = -log2(p);
		if(verbose == 3) printf("%s Collision Estimate: Found p.\n", label);
	} else {
		if(verbose == 3) printf("%s Collision Estimate: Could Not Find p. Proceeding with the lower bound for p.\n", label);
		p = 0.5;
		entEst = 1.0;
	}

	if(verbose == 2) printf("p = %.17g\n", p);
	else if(verbose == 3) {
		printf("%s Collision Estimate: p = %.17g\n", label, p);
		printf("%s Collision Estimate: min entropy = %.17g\n", label, entEst);
	}

	return entEst;
}
