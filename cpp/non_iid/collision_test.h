#pragma once
#include "../shared/utils.h"

// Computed using efficient implementation in Appendix G.1.1
double F(double q){
   return q*(2.0*q*q+2.0*q+1.0);
}

double col_exp(double p, double q){
	return (p/(q*q))*(1.0 + 0.5*(1.0/p - 1.0/q))*F(q) - (p/q)*0.5*(1.0/p - 1.0/q);
}

// Section 6.3.2 - Collision Estimate
// data is assumed to be binary (e.g., bit string)
double collision_test(byte* data, long len){
	long v, i, j;
	int t_v;
	double X, s, p, lastP, pVal;
        double lvalue, hvalue;
        double hbound, lbound;
        double hdomain, ldomain;

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
	X -= ZALPHA * s/sqrt(v);

	ldomain = 0.5;
        hdomain = 1.0;

        lbound = ldomain;
        hbound = hdomain;

        lvalue = DBL_INFINITY;
        hvalue = -DBL_INFINITY;

        //Note that the bounds are in [0,1], so overflows aren't an issue
        //But underflows are.
        p = (lbound + hbound) / 2.0;
	pVal = col_exp(p, 1.0-p);

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
			//This is a search failure (as directed in step #8)
                        p = 0.5;
                        break;
                }

                //invariant. If this isn't true, then seeking the value within this interval doesn't make sense.
                if(!INCLOSEDINTERVAL(X, lvalue, hvalue)) {
			//This is a search failure (as directed in step #8)
                        p = 0.5;
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

		pVal = col_exp(p, 1.0-p);

                //invariant: If this isn't true, then this isn't loosly monotonic
                if(!INCLOSEDINTERVAL(pVal, lvalue, hvalue)) {
                        p = hbound;
                        break;
                }
        }//for loop

	return -log2(p);
}
