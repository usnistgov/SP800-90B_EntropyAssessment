#pragma once
#include "../shared/utils.h"

// Computed using efficient implementation in Appendix G.1.1
double F(double q){
	int i;	
	double z, ret;
	
	i = 1000;
	z = 1.0/q;
	ret = z + (i-2.0)/(i+1.0);
	
	while(i > 0){
		i--;
		ret = z + (i-2.0) / (1.0 + ((i+1.0)/ret));
	}

	return  1.0/ret;
}

double col_exp(double p, double q){
	return (p/(q*q))*(1.0 + 0.5*(1.0/p - 1.0/q))*F(q) - (p/q)*0.5*(1.0/p - 1.0/q);
}

// Section 6.3.2 - Collision Estimate
// data is assumed to be binary (e.g., bit string)
double collision_test(byte* data, long len){
	long v, i;
	int t_v;
	double X, s, p, p_hi, p_lo, eps, exp;
	
	i = 0;
	v = 0;
	X = 0.0;
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
	s = sqrt((s - (i*X)) / (v-1));

	// binary search for p
	X -= 2.576 * s/sqrt(v);
	eps = 1.0 / (1 << 20); // 2^-20
	p_lo = 0.5;
	p_hi = 1.0 - eps; // avoid division by zero
	do{
		p = (p_lo + p_hi) / 2.0;
		exp = col_exp(p, 1.0-p);

		if(X < exp) p_lo = p;
		else p_hi = p;

		if((col_exp(p_lo, 1.0-p_lo) < X) || (col_exp(p_hi, 1.0-p_hi) > X)){
			// binary search failed, settle for one bit of entropy
			printf("\t *** WARNING: binary search for collision test failed ***\n");
			p = 0.5; 
			break;
		}
	}while(fabs(p_hi - p_lo) > eps);

	return -log2(p);
}
