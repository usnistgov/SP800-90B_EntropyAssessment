#pragma once

#include "utils.h"
#include <divsufsort.h>

#define SAINDEX_MAX INT32_MAX

//Using the Kasai (et. al) O(n) time "13n space" algorithm.
//"Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications", by Kasai, Lee, Arimura, Arikawa, and Park
//https://doi.org/10.1007/3-540-48194-X_17
//http://web.cs.iastate.edu/~cs548/references/linear_lcp.pdf
//The default implementation uses 4 byte indexes
//Note that indexes should be signed, so the next natural size is int64_t
static void sa2lcp(const byte text[], long int n, const vector<saidx_t> &sa, vector<saidx_t> &lcp) {
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

void calcSALCP(const byte text[], long int n, vector<saidx_t> &sa, vector<saidx_t> &lcp) {
	int32_t res;

	assert(n < INT32_MAX); //This is the default type, but it can be compiled to use 64 bit indexes (and then this should be INT64_MAX)
	assert(n > 0); //This is the default type, but it can be compiled to use 64 bit indexes (and then this should be INT64_MAX)
	assert(sa.size() == (size_t)(n+1));
	assert(sa.size() == (size_t)(n+1));

	sa[0] = (saidx_t)n;

	res=divsufsort((const sauchar_t *)text, (saidx_t *)(sa.data()+1), (saidx_t)n);
	assert(res==0);
   	sa2lcp(text, n, sa, lcp);
}

/* Based on the algorithm outlined by Aaron Kaufer
 * This is described here:
 * http://www.untruth.org/~josh/sp80090b/Kaufer%20Further%20Improvements%20for%20SP%20800-90B%20Tuple%20Counts.pdf
 */
void SAalgs(const byte text[], long int n, int k, double &t_tuple_res, double &lrs_res, bool verbose) {
	vector <saidx_t> sa(n+1, -1); //each value is at most n-1
	vector <saidx_t> L(n+2, -1); //each value is at most n-1

   	long int u; //The length of a string: 1 <= u <= v+1 <= n
   	long int v; //The length of the LRS. 1 <= v <= n-1
	long int c; //contains a count from A
	long int j; //0 <= j <= v+1 <= n
	saidx_t t; //Takes values from LCP array. 0 <= t < n

	double curPmax;
	double Pmax;
	double pu;

	assert(n>0);
	assert(k>0);
	assert(n <= SAINDEX_MAX - 1);

	calcSALCP(text, n, sa, L);

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
		double curP = ((double)(Q[i]))/((double)(n-i+1));
		double curPMax = pow(curP, 1.0/(double)i);
		 //fprintf(stderr, "t-Tuple Estimate: P[%ld] = %.17g ( %d / %ld )\n", i, curP, Q[i], n-i+1);
		 //fprintf(stderr, "t-Tuple Estimate: P_max[%ld] = %.17g\n", i, curPMax);

		if(curPMax > Pmax) {
			Pmax=curPMax;
		}
	}

	//finalize the t-tuple estimate
	if(Pmax > 0.0) {
		//We encountered a valid t, so we can run the test
		pu = Pmax + ZALPHA*sqrt(Pmax*(1.0 - Pmax)/((double)(n - 1)));
		if(verbose) printf("t-Tuple Estimate: t = %ld, p-hat_max = %.17g, p_u = %.17g\n", u-1, Pmax, pu);
		if(pu > 1.0) {
			pu = 1.0;
		}

		t_tuple_res = -log2(pu);
	} else {
		if(verbose) printf("t-Tuple Estimate: No strings are repeated 35 times. t-Tuple estimate failed.\n");
		t_tuple_res = -1.0;
	}

	//calculate the LRS estimate
	if(v>=u) {
		vector <long int> S(v+1, 0);
		memset(A.data(), 0, sizeof(saidx_t)*sizeof((v+2)));

		for(long int i = 1; i <= n; i++) {
			if((L[i-1] >= u) && (L[i] < L[i-1])) {
				saidx_t b = L[i];

				//A[u] stores the number of u-length tuples. We need to eventually clear down to A[u]=A[b+1].
				if(b < u) b = u-1;

				for(t = L[i-1]; t > b; t--) {
					A[t] += A[t+1];
					A[t+1] = 0;

					assert(A[t] >= 0);
					//update sum
					//Note that (c choose 2) is just (c)(c-1)/2.
					//The numerator of this expression is necessarily even
					//Dividing an even quantity by 2 is the same as right shifting by 1.
					S[t] += ((((uint64_t)(A[t]+1) * (uint64_t)(A[t]))))>>1; /* update sum */ 
				}

				if(b >= u) A[b] += A[b+1]; /* carry over count for t = L[i] */
				A[b+1] = 0;
			}

			if(L[i] >= u) A[L[i]]++; /* update count for t = L[i] */
		}

		//We now have a complete set of numerators in S
		Pmax = 0.0;
		for(long int i=u; i<=v; i++) {
			double curP = ((double)S[i]) / (double)(((n-i)*(n-i+1))>>1);
			double curPMax = pow(curP, 1.0/((double)i));
			 //fprintf(stderr, "LRS Estimate: P_%ld = %.17g ( %zu / %zu )\n", i, curP, S[i], ((n-i)*(n-i+1))>>1);
			 //fprintf(stderr, "LRS Estimate: P_{max,%ld} = %.17g\n", i, curPMax);


			if(Pmax < curPMax) {
				Pmax = curPMax;
			}
		}

		pu = Pmax + ZALPHA*sqrt(Pmax*(1.0 - Pmax)/((double)(n - 1)));
		if(pu > 1.0) {
			pu = 1.0;
		}
		if(verbose) printf("LRS Estimate: u = %ld, v = %ld, p-hat = %.17g, p_u = %.17g\n", u, v, Pmax, pu);

		lrs_res = -log2(pu);
	} else {
		printf("LRS Estimate: v<u. Can't Run LRS Test.\n");
		lrs_res = -1.0;
		return;
	}

	return;
}

int len_LRS(const byte text[], const int sample_size){
	vector <saidx_t> sa(sample_size+1, -1);
	vector <saidx_t> lcp(sample_size+1, -1);
	saidx_t lrs_len = -1;

	calcSALCP(text, sample_size, sa, lcp);

	for(saidx_t j = 0; j <= sample_size; j++) {
		if(lcp[j] > lrs_len) lrs_len = lcp[j];
	}

	return(lrs_len);
}

/*
* ---------------------------------------------
*			 HELPER FUNCTIONS
* ---------------------------------------------
*/

void calc_collision_proportion(const vector<double> &p, double &p_col){
	
	for(unsigned int i = 0; i < p.size(); i++){
		p_col += pow(p[i], 2);
	}
}

/*
* ---------------------------------------------
* 		  			 TEST
* ---------------------------------------------
*/

bool len_LRS_test(const byte data[], const int sample_size, const int alphabet_size, const bool verbose){

	vector<double> p(alphabet_size, 0.0);
	calc_proportions(data, p, sample_size);

	double p_col = 0.0;
	calc_collision_proportion(p, p_col);

	int lrs = len_LRS(data, sample_size);
	int n = sample_size - lrs + 1;
	long int overlap = n_choose_2(n);

	double pr_x = 1 - pow(1 - pow(p_col, lrs), overlap);

	if(verbose){
		cout << "Longest Repeated Substring results" << endl;
		cout << "\tP_col: " << p_col << endl;
		cout << "\tLength of LRS: " << lrs << endl;
		cout << "\tPr(X >= 1): " << pr_x << endl;
	}

	return (pr_x >= 0.001);
}
