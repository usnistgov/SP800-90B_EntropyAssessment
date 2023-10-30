#pragma once

#include "utils.h"
#include <climits>
#include <divsufsort.h>
#include <divsufsort64.h>

#define SAINDEX_MAX INT32_MAX
#define SAINDEX64_MAX INT64_MAX

//Using the Kasai (et al.) O(n) time "13n space" algorithm.
//"Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications", by Kasai, Lee, Arimura, Arikawa, and Park
//https://doi.org/10.1007/3-540-48194-X_17
//http://web.cs.iastate.edu/~cs548/references/linear_lcp.pdf
//The default implementation uses 4 byte indexes
static void sa2lcp32(const uint8_t text[], long int n, const vector<saidx_t> &sa, vector<saidx_t> &lcp) {
	saidx_t h;
	vector<saidx_t> rank(n+1,-1);

	assert(n>1);

	lcp[0] = -1;
	lcp[1] = 0;

	// compute rank = sa^{-1}
	for(saidx_t i=0; i<=(saidx_t)n; i++) {
		rank[sa[i]] = i;
	}

	// traverse suffixes in rank order
	h=0;

	for(saidx_t i=0; i<(saidx_t)n; i++) {
		saidx_t k = rank[i]; // rank of s[i ... n-1]
		if(k>1) {
			saidx_t j = sa[k-1]; // predecessor of s[i ... n-1]
			while((i+h<(saidx_t)n) && (j+h<(saidx_t)n) && (text[i+h]==text[j+h])) {
				h++;
			}

			lcp[k] = h;
		}
		if(h>0) {
			h--;
		}
	}
}

//Using the Kasai (et al.) O(n) time "25n space" algorithm (with 64-bit indicies)
static void sa2lcp64(const uint8_t text[], long int n, const vector<saidx64_t> &sa, vector<saidx64_t> &lcp) {
	saidx64_t h;
	vector<saidx64_t> rank(n+1,-1);

	assert(n>1);

	lcp[0] = -1;
	lcp[1] = 0;

	// compute rank = sa^{-1}
	for(saidx64_t i=0; i<=(saidx64_t)n; i++) {
		rank[sa[i]] = i;
	}

	// traverse suffixes in rank order
	h=0;

	for(saidx64_t i=0; i<(saidx64_t)n; i++) {
		saidx64_t k = rank[i]; // rank of s[i ... n-1]
		if(k>1) {
			saidx64_t j = sa[k-1]; // predecessor of s[i ... n-1]
			while((i+h<(saidx64_t)n) && (j+h<(saidx64_t)n) && (text[i+h]==text[j+h])) {
				h++;
			}

			lcp[k] = h;
		}
		if(h>0) {
			h--;
		}
	}
}


void calcSALCP32(const uint8_t text[], long int n, vector<saidx_t> &sa, vector<saidx_t> &lcp) {
	int32_t res;

	assert(n < SAINDEX_MAX);
	assert(n > 0); 
	assert(sa.size() == (size_t)(n+1));
	assert(sa.size() == (size_t)(n+1));

	sa[0] = (saidx_t)n;

	res=divsufsort((const sauchar_t *)text, (saidx_t *)(sa.data()+1), (saidx_t)n);
	assert(res==0);
   	sa2lcp32(text, n, sa, lcp);
}

void calcSALCP64(const uint8_t text[], long int n, vector<saidx64_t> &sa, vector<saidx64_t> &lcp) {
	int32_t res;

	assert(n < SAINDEX64_MAX);
	assert(n > 0);
	assert(sa.size() == (size_t)(n+1));
	assert(sa.size() == (size_t)(n+1));

	sa[0] = (saidx64_t)n;

	res=divsufsort64((const sauchar_t *)text, (saidx64_t *)(sa.data()+1), (saidx64_t)n);
	assert(res==0);
   	sa2lcp64(text, n, sa, lcp);
}
/* Based on the algorithm outlined by Aaron Kaufer
 * This is described here:
 * http://www.untruth.org/~josh/sp80090b/Kaufer%20Further%20Improvements%20for%20SP%20800-90B%20Tuple%20Counts.pdf
 */
void SAalgs32(const uint8_t text[], long int n, int k, double &t_tuple_res, double &lrs_res, const int verbose, const char *label) {
	vector <saidx_t> sa(n+1, -1); //each value is at most n-1
	vector <saidx_t> L(n+2, -1); //each value is at most n-1

   	long int u; //The length of a string: 1 <= u <= v+1 <= n
   	long int v; //The length of the LRS. 1 <= v <= n-1
	long int c; //contains a count from A
	long int j; //0 <= j <= v+1 <= n
	saidx_t t; //Takes values from LCP array. 0 <= t < n

	long double Pmax;
	long double pu;

	assert(n>0);
	assert(k>0);
	assert(n < SAINDEX_MAX);
	assert((UINT64_MAX / (uint64_t)n) >= ((uint64_t)n+1U)); // (mult assert)

	calcSALCP32(text, n, sa, L);

	//to conform with Kaufer's conventions
	L.erase(L.begin());
	L[n] = 0;
	assert(L[0] == 0);

	//Find the length of the LRS, v

	v=0;
	for(long int i=0; i<n; i++) {
		if(L[i]>v) v = L[i];
	}

	assert((v>0) && (v < n));
	//v is now set correctly

	vector <saidx_t> Q(v+1, 1); //Contains an accumulation of positive counts 1 <= Q[i] <= n
	vector <saidx_t> A(v+2, 0); //Contains an accumulation of positive counts 0 <= A[i] <= n
	//I is set from L
	//Note that I is indexed by at most j+1.
	// j takes the value 0 to v+1  (so I[v+2] should work)
	//(I stores indices of A, and there are only v+2 of these)
	vector <saidx_t> I(v+3, 0); //each value is most 0 <= I[i] <= v+2 <= n+1

	j = 0;
	for(long int i = 1; i <= n; i++) {
		c = 0;
		//Note L[0] is already verified to be 0
		assert(L[i] >= 0);

		if(L[i] < L[i-1]) {
			t = L[i-1];
			assert(j>0);
			j--;
			assert(j<=v);

			while(t > L[i]) {
				assert((t>0) && (t <= v));
				if((j > 0) && (I[j] == t)) {
					/* update count for non-zero entry of A */
					A[I[j]] += A[I[j+1]];
					A[I[j+1]] = 0;
					j--;
				}

				if(Q[t] >= A[I[j+1]]+1) {
					/*
					 * Q[t] is at least as large as current count,
					 * and since Q[t] <= Q[t-1] <= ... <= Q[1],
					 * there is no need to check zero entries of A
					 * until next non-zero entry
					 */
					if(j > 0) {
						/* skip to next non-zero entry of A */
						t = I[j];
					} else {
						/*
						 * no more non-zero entries of A,
						 * so skip to L[i] (terminate while loop)
						 */
						t = L[i];
					}
				} else {
					/* update Q[t] with new maximum count */
					Q[t--] = A[I[j+1]]+1;
				}
			}

			c = A[I[j+1]]; /* store carry over count */
			A[I[j+1]] = 0;
		}

		if(L[i] > 0) {
			if((j < 1) || (I[j] < L[i])) {
				/* insert index of next non-zero entry of A */
				assert(j<v);
				I[++j] = L[i];
			}
			A[I[j]] += c+1; /* update count for t = I[j] = L[i] */
		}
	}

	//Calculate u
	for(u=1; (u<=v) && (Q[u] >= 35); u++);

	assert(u > 0);
	assert(((u == v+1) || ((u <= v) && (Q[u] < 35))));
	assert(((u == 1) || (Q[u-1] >= 35)));
	//u is now correctly set.

	//at this point, Q is completely calculated.
	/*Calculate the various Pmax[i] values. We need not save the actual values, only the largest*/
	Pmax = -1.0;
	for(long int i=1; i<u; i++) {
		long double curP = ((long double)(Q[i]))/((long double)(n-i+1));
		long double curPMax = powl(curP, 1.0L/(long double)i);

		if(curPMax > Pmax) {
			Pmax=curPMax;
		}
	}

	//finalize the t-tuple estimate
	if(Pmax > 0.0L) {
		//We encountered a valid t, so we can run the test
		pu = Pmax + ZALPHA_L*sqrtl(Pmax*(1.0L - Pmax)/((long double)(n - 1)));
		if(pu > 1.0L) {
			pu = 1.0L;
		}

		t_tuple_res = (double)-log2l(pu);

		if(verbose == 2) printf("%s t-Tuple Estimate: t = %ld, p-hat_max = %.22Lg, p_u = %.22Lg\n", label, u-1, Pmax, pu);
		else if(verbose == 3) {
			printf("%s t-Tuple Estimate: t = %ld\n", label, u-1);
			printf("%s t-Tuple Estimate: p-hat_max = %.22Lg\n", label, Pmax);
			printf("%s t-Tuple Estimate: p_u = %.22Lg\n", label, pu);
			printf("%s t-Tuple Estimate: min entropy = %.17g\n", label, t_tuple_res);
		}

	} else {
		if(verbose > 1) printf("t-Tuple Estimate: No strings are repeated 35 times. t-Tuple estimate failed.\n");
		t_tuple_res = -1.0;
	}

	//calculate the LRS estimate
	if(v>=u) {
		vector <uint64_t> S(v+1, 0);
		memset(A.data(), 0, sizeof(saidx_t)*((size_t)v+2));

		for(long int i = 1; i <= n; i++) {
			if((L[i-1] >= u) && (L[i] < L[i-1])) {
				saidx_t b = L[i];

				//A[u] stores the number of u-length tuples. We need to eventually clear down to A[u]=A[b+1].
				if(b < u) b = u-1;

				for(t = L[i-1]; t > b; t--) {
					uint64_t priorS;
					uint64_t choices;
					A[t] += A[t+1];
					A[t+1] = 0;

					assert(A[t] >= 0);
					// update sum
					// Note that (c choose 2) is just (c)(c-1)/2.
					// The numerator of this expression is necessarily even
					// Dividing an even quantity by 2 is the same as right shifting by 1.
					// Check for overflows when adding to S[t] (unsigned 64 bit integers)
					// Note, A[t] <= n, so the assert marked "(mult assert)" tells us that the multiplication won't rollover.

					priorS = S[t];
					choices = ((((uint64_t)(A[t]+1) * (uint64_t)(A[t]))))>>1;
					S[t] = priorS + choices;
					assert(S[t] >= priorS);
				}

				if(b >= u) A[b] += A[b+1]; /* carry over count for t = L[i] */
				A[b+1] = 0;
			}

			if(L[i] >= u) A[L[i]]++; /* update count for t = L[i] */
		}

		//We now have a complete set of numerators in S
		Pmax = 0.0;
		for(long int i=u; i<=v; i++) {
			// Note, the assert marked "(mult assert)" tells us that the multiplication won't rollover.
			uint64_t choices = (((uint64_t)n-(uint64_t)i)*((uint64_t)n-(uint64_t)i+1U))>>1;
			long double curP = ((long double)S[i]) / (long double)choices;
			long double curPMax = pow(curP, 1.0/((long double)i));

			if(Pmax < curPMax) {
				Pmax = curPMax;
			}
		}

		pu = Pmax + ZALPHA_L*sqrtl(Pmax*(1.0L - Pmax)/((long double)(n - 1)));
		if(pu > 1.0L) {
			pu = 1.0L;
		}

		lrs_res = (double)-log2l(pu);

		if(verbose == 2) printf("%s LRS Estimate: u = %ld, v = %ld, p-hat = %.17Lg, p_u = %.17Lg\n", label, u, v, Pmax, pu);
		else if(verbose == 3) {
			printf("%s LRS Estimate: u = %ld\n", label, u);
			printf("%s LRS Estimate: v = %ld\n", label, v);
			printf("%s LRS Estimate: p-hat = %.22Lg\n", label, Pmax);
			printf("%s LRS Estimate: p_u = %.22Lg\n", label, pu);
			printf("%s LRS Estimate: min entropy = %.17g\n", label, lrs_res);
		}

	} else {
		printf("LRS Estimate: v<u. Can't Run LRS Test.\n");
		lrs_res = -1.0;
		return;
	}

	return;
}

void SAalgs64(const uint8_t text[], long int n, int k, double &t_tuple_res, double &lrs_res, const int verbose, const char *label)
{
	vector <saidx64_t> sa(n+1, -1); //each value is at most n-1
	vector <saidx64_t> L(n+2, -1); //each value is at most n-1

   	long int u; //The length of a string: 1 <= u <= v+1 <= n
   	long int v; //The length of the LRS. 1 <= v <= n-1
	long int c; //contains a count from A
	long int j; //0 <= j <= v+1 <= n
	saidx64_t t; //Takes values from LCP array. 0 <= t < n

	long double Pmax;
	long double pu;

	assert(n>0);
	assert(k>0);
	assert(n <= SAINDEX64_MAX - 1);
	assert((UINT128_MAX / (uint128_t)n) >= ((uint128_t)n+1U)); // (mult assert)

	calcSALCP64(text, n, sa, L);

	//to conform with Kaufer's conventions
	L.erase(L.begin());
	L[n] = 0;
	assert(L[0] == 0);

	//Find the length of the LRS, v

	v=0;
	for(long int i=0; i<n; i++) {
		if(L[i]>v) v = L[i];
	}

	assert((v>0) && (v < n));
	//v is now set correctly

	vector <saidx64_t> Q(v+1, 1); //Contains an accumulation of positive counts 1 <= Q[i] <= n
	vector <saidx64_t> A(v+2, 0); //Contains an accumulation of positive counts 0 <= A[i] <= n
	//I is set from L
	//Note that I is indexed by at most j+1.
	// j takes the value 0 to v+1  (so I[v+2] should work)
	//(I stores indices of A, and there are only v+2 of these)
	vector <saidx64_t> I(v+3, 0); //each value is most 0 <= I[i] <= v+2 <= n+1

	j = 0;
	for(long int i = 1; i <= n; i++) {
		c = 0;
		//Note L[0] is already verified to be 0
		assert(L[i] >= 0);

		if(L[i] < L[i-1]) {
			t = L[i-1];
			assert(j>0);
			j--;
			assert(j<=v);

			while(t > L[i]) {
				assert((t>0) && (t <= v));
				if((j > 0) && (I[j] == t)) {
					/* update count for non-zero entry of A */
					A[I[j]] += A[I[j+1]];
					A[I[j+1]] = 0;
					j--;
				}

				if(Q[t] >= A[I[j+1]]+1) {
					/*
					 * Q[t] is at least as large as current count,
					 * and since Q[t] <= Q[t-1] <= ... <= Q[1],
					 * there is no need to check zero entries of A
					 * until next non-zero entry
					 */
					if(j > 0) {
						/* skip to next non-zero entry of A */
						t = I[j];
					} else {
						/*
						 * no more non-zero entries of A,
						 * so skip to L[i] (terminate while loop)
						 */
						t = L[i];
					}
				} else {
					/* update Q[t] with new maximum count */
					Q[t--] = A[I[j+1]]+1;
				}
			}

			c = A[I[j+1]]; /* store carry over count */
			A[I[j+1]] = 0;
		}

		if(L[i] > 0) {
			if((j < 1) || (I[j] < L[i])) {
				/* insert index of next non-zero entry of A */
				assert(j<v);
				I[++j] = L[i];
			}
			A[I[j]] += c+1; /* update count for t = I[j] = L[i] */
		}
	}

	//Calculate u
	for(u=1; (u<=v) && (Q[u] >= 35); u++);

	assert(u > 0);
	assert(((u == v+1) || ((u <= v) && (Q[u] < 35))));
	assert(((u == 1) || (Q[u-1] >= 35)));
	//u is now correctly set.

	//at this point, Q is completely calculated.
	/*Calculate the various Pmax[i] values. We need not save the actual values, only the largest*/
	Pmax = -1.0;
	for(long int i=1; i<u; i++) {
		long double curP = ((long double)(Q[i]))/((long double)(n-i+1));
		long double curPMax = powl(curP, 1.0L/(long double)i);

		if(curPMax > Pmax) {
			Pmax=curPMax;
		}
	}

	//finalize the t-tuple estimate
	if(Pmax > 0.0) {
		//We encountered a valid t, so we can run the test
		pu = Pmax + ZALPHA_L*sqrtl(Pmax*(1.0L - Pmax)/((long double)(n - 1)));
		if(pu > 1.0) {
			pu = 1.0;
		}

		t_tuple_res = (double)-log2l(pu);

		if(verbose == 2) printf("%s t-Tuple Estimate: t = %ld, p-hat_max = %.22Lg, p_u = %.22Lg\n", label, u-1, Pmax, pu);
		else if(verbose == 3) {
			printf("%s t-Tuple Estimate: t = %ld\n", label, u-1);
			printf("%s t-Tuple Estimate: p-hat_max = %.22Lg\n", label, Pmax);
			printf("%s t-Tuple Estimate: p_u = %.22Lg\n", label, pu);
			printf("%s t-Tuple Estimate: min entropy = %.17g\n", label, t_tuple_res);
		}

	} else {
		if(verbose > 1) printf("t-Tuple Estimate: No strings are repeated 35 times. t-Tuple estimate failed.\n");
		t_tuple_res = -1.0;
	}

	//calculate the LRS estimate
	if(v>=u) {
		vector <uint128_t> S(v+1, 0);
		memset(A.data(), 0, sizeof(saidx64_t)*((size_t)v+2));

		for(long int i = 1; i <= n; i++) {
			if((L[i-1] >= u) && (L[i] < L[i-1])) {
				saidx64_t b = L[i];

				//A[u] stores the number of u-length tuples. We need to eventually clear down to A[u]=A[b+1].
				if(b < u) b = u-1;

				for(t = L[i-1]; t > b; t--) {
					uint128_t priorS;
					uint128_t choices;
					A[t] += A[t+1];
					A[t+1] = 0;

					assert(A[t] >= 0);
					// update sum
					// Note that (c choose 2) is just (c)(c-1)/2.
					// The numerator of this expression is necessarily even
					// Dividing an even quantity by 2 is the same as right shifting by 1.
					// Check for overflows when adding to S[t] (unsigned 64 bit integers)
					// Note, A[t] <= n, so the assert marked "(mult assert)" tells us that the multiplication won't rollover.

					priorS = S[t];
					choices = ((((uint128_t)(A[t]+1) * (uint128_t)(A[t]))))>>1;
					S[t] = priorS + choices;
					assert(S[t] >= priorS);
				}

				if(b >= u) A[b] += A[b+1]; /* carry over count for t = L[i] */
				A[b+1] = 0;
			}

			if(L[i] >= u) A[L[i]]++; /* update count for t = L[i] */
		}

		//We now have a complete set of numerators in S
		Pmax = 0.0;
		for(long int i=u; i<=v; i++) {
			// Note, the assert marked "(mult assert)" tells us that the multiplication won't rollover.
			uint128_t choices = (((uint128_t)n-(uint128_t)i)*((uint128_t)n-(uint128_t)i+1U))>>1;
			long double curP = ((long double)S[i]) / (long double)choices;
			long double curPMax = powl(curP, 1.0L/((long double)i));


			if(Pmax < curPMax) {
				Pmax = curPMax;
			}
		}

		pu = Pmax + ZALPHA_L*sqrtl(Pmax*(1.0L - Pmax)/((long double)(n - 1)));
		if(pu > 1.0) {
			pu = 1.0;
		}

		lrs_res = -log2(pu);

		if(verbose == 2) printf("%s LRS Estimate: u = %ld, v = %ld, p-hat = %.22Lg, p_u = %.22Lg\n", label, u, v, Pmax, pu);
		else if(verbose == 3) {
			printf("%s LRS Estimate: u = %ld\n", label, u);
			printf("%s LRS Estimate: v = %ld\n", label, v);
			printf("%s LRS Estimate: p-hat = %.22Lg\n", label, Pmax);
			printf("%s LRS Estimate: p_u = %.22Lg\n", label, pu);
			printf("%s LRS Estimate: min entropy = %.17g\n", label, lrs_res);
		}

	} else {
		printf("LRS Estimate: v<u. Can't Run LRS Test.\n");
		lrs_res = -1.0;
		return;
	}

	return;
}

void SAalgs(const uint8_t text[], long int n, int k, double &t_tuple_res, double &lrs_res, const int verbose, const char *label) {
	if(n<SAINDEX_MAX) {
		SAalgs32(text, n, k, t_tuple_res, lrs_res, verbose, label);
	} else {
		SAalgs64(text, n, k, t_tuple_res, lrs_res, verbose, label);
	}
}

long int len_LRS32(const uint8_t text[], const int sample_size){
	vector <saidx_t> sa(sample_size+1, -1);
	vector <saidx_t> lcp(sample_size+1, -1);
	saidx_t lrs_len = -1;

	calcSALCP32(text, sample_size, sa, lcp);

	for(saidx_t j = 0; j <= sample_size; j++) {
		if(lcp[j] > lrs_len) lrs_len = lcp[j];
	}

	return(lrs_len);
}

long int len_LRS64(const uint8_t text[], const int sample_size){
	vector <saidx64_t> sa(sample_size+1, -1);
	vector <saidx64_t> lcp(sample_size+1, -1);
	saidx64_t lrs_len = -1;

	calcSALCP64(text, sample_size, sa, lcp);

	for(saidx64_t j = 0; j <= sample_size; j++) {
		if(lcp[j] > lrs_len) lrs_len = lcp[j];
	}

	return(lrs_len);
}
/*
* ---------------------------------------------
*			 HELPER FUNCTIONS
* ---------------------------------------------
*/

void calc_collision_proportion(const vector<double> &p, long double &p_col){

	p_col = 0.0L;
	
	for(unsigned int i = 0; i < p.size(); i++){
		p_col += powl((long double)(p[i]), 2.0L);
	}
}

/*
* ---------------------------------------------
* 		  			 TEST
* ---------------------------------------------
*/

bool len_LRS_test(const uint8_t data[], const int L, const int k, const int verbose, const char *label) {
	// p_col is the probability of collision on a per-symbol basis under an IID assumption (this is related to the collision entropy).
	// p_col >= 1/k, which bounds this.
	// Note, for SP 800-90B k<=256, so we can bound p_col >= 2^-8.
	vector<double> p(k, 0.0);
	calc_proportions(data, p, L);
	long double p_col = 0.0;
	calc_collision_proportion(p, p_col);

	assert(p_col >= 1.0L / ((long double) k));
	assert(p_col <= 1.0L);

	// It is possible for p_col to be exactly 1 (e.g., if the input data is all one symbol)
	// In this instance, a collision of any length up to L-1 has probability 1.
	if(p_col > 1.0L - LDBL_EPSILON) {
		if(verbose == 2) {
			printf("\tPr(X >= 1) = 1.0\n");
		} else if(verbose > 2) {
			printf("%s Longest Repeated Substring: P_col = 1.0\n", label);
			printf("%s Longest Repeated Substring: Pr(X >= 1) = 1.0\n", label);
		} 
		return true;
	}
	assert(p_col < 1.0L);

	// The length of the longest repeated substring (LRS) for the supplied data is W.
	long int W;
	if(L<SAINDEX_MAX) {
		W = len_LRS32(data, L);
	} else {
		W = len_LRS64(data, L);
	}

	// p_col^W is the probability of collision of a W-length string under an IID assumption;
	// this may be quite close to 0.
	// We know that p_col >= 2^-8, but if W is sufficiently large, then the result will be smaller than LDBL_MIN 
	// (2^-16382 on modern Intel platforms); if this happens then we can't represent p_col^W directly as a long double.
	// There isn't anything we can do about this error condition without moving to an arbitrary precision calculation, so
	// for now, we'll just detect this condition and abort the calculation.
	long double p_colPower = powl(p_col, (long double)W);
	assert(p_colPower >= LDBL_MIN);
	assert(p_colPower <= 1.0L-LDBL_EPSILON);

	// There is some delicacy in calculating Pr(X>=1) as some of the intermediary values may be quite close to 0 or 1.
	// We first want to calculate the probability of not having a collision of length W for a single pair of independent 
	// W-length strings.  This quantity is (1-p_col^W).
	// We know that p_col >= 2^-8, but if W is sufficiently large, then the result will be smaller than 
	// LDBL_EPSILON (2^-63 on modern Intel platforms); if this happens then we can't represent (1-p_col^W) directly as a long double.
	// We instead calculate the log of the probability of not having a collision of length W for a single pair of independent W-length strings.
	// Recall that log1p(x) = log(1+x); this form is useful when |x| is small. We are particularly concerned with the case where p_col^W is a small 
	// positive value, which would make log1p(-p_col^W) a negative value quite close to 0.
	long double logProbNoColsPerPair = log1pl(-p_colPower);
	assert(logProbNoColsPerPair < 0.0L);

	//(L - W + 1) is the number of overlapping contiguous substrings of length W in a string of length L.
	// The number of pairs of such overlapping substrings is N = (L - W + 1) choose 2.
	// This is the number of ways of choosing 2 substrings of length W from a string of length L.
	long int N = n_choose_2(L - W + 1);

	if(verbose > 1) {
		if(verbose > 2) {
			printf("%s Longest Repeated Substring results: P_col = %.22Lg\n", label, p_col);
			printf("%s Longest Repeated Substring results: W = %ld\n", label, W);
		} else {
			printf("%s Longest Repeated Substring results\n", label);
			printf("\tP_col: %Lf\n", p_col);
			printf("\tLength of LRS: %ld\n", W);
		}

		// Calculate the probability of not encountering a collision after N sets of independent pairs;
		// this is an application of the Binomial Distribution.
		// Using the CDF at X=0 (or equivalently, the PDF at X=0), we find the probability of seeing at least 1 collision,
		// under the assumption that each substring is independent of every other substring (i.e., under the per-round independence 
		// assumption required by the Binomial Distribution), is:
		// probNoCols = 1 - Pr (X = 0) = 1 - (1 - p_col^W)^N. 
		// Unfortunately, if W>1, many of these substrings are not independent! In the case that two strings overlap, then
		// there is clearly a dependency between the trials. As such, this test has some statistical construction problems.  
		// Empirical testing has shown that, for ideal looking data, the observed failure rate is under the desired 
		// threshold of 1/1000.
		// Note, this N is O(L^2), so use of this value as an exponent tends to cause underflows here;
		// in this case, this probability isn't accurately representable using the precision that we have to work with, but it is expected to
		// round reasonably.
		long double probNoCols = expl(((long double)N)*logProbNoColsPerPair);
		if(verbose > 2) {
			printf("%s Longest Repeated Substring results: Pr(X >= 1) = %.22Lg\n",  label, 1.0L - probNoCols);
		} else {
			printf("\tPr(X >= 1): %Lf\n", 1.0L - probNoCols);
		}
	}

	//We don't need the above to have worked to come to a conclusion on the test, however:
	// The LRS test is considered a "Pass"
	// iff Pr(X>=1) >= 1/1000
	// iff 1 - Pr(X=0) >= 1/1000
	// iff 1 - (1-p_col^W)^N >= 1/1000
	// iff 0.999 >= (1-p_col^W)^N
	// iff log(0.999) >= N*log(1-p_col^W)
	// iff log(0.999) >= N*log1p(-p_col^W)
	return logl(0.999L) >= ((long double)N)* logProbNoColsPerPair;
}
