#pragma once

#include "utils.h"
#include <divsufsort.h>

//Using the Kasai (et. al) O(n) time "13n space" algorithm.
//"Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications", by Kasai, Lee, Arimura, Arikawa, and Park
//https://doi.org/10.1007/3-540-48194-X_17
//http://web.cs.iastate.edu/~cs548/references/linear_lcp.pdf
//The default implementation uses 4 byte indexes
//Note that indexes should be signed, so the next natural size is int64_t
static void sa2lcp(const byte text[], long int n, const vector<saidx_t> &sa, vector<saidx_t> &lcp){
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

void calcSALCP(const byte text[], long int n, int k, vector<saidx_t> &sa, vector<saidx_t> &lcp){
	int32_t res;

	assert(n < INT32_MAX); //This is the default type, but it can be compiled to use 64 bit indexes (and then this should be INT64_MAX)
	assert(sa.size() == n+1);
	assert(sa.size() == n+1);

	sa[0] = (saidx_t)n;

	res=divsufsort((const sauchar_t *)text, (saidx_t *)(sa.data()+1), (saidx_t)n);
	assert(res==0);
   	sa2lcp(text, n, sa, lcp);
}


static void countTuples(const vector<saidx_t> &sa, const vector<saidx_t> &lcp, long int L, long int tupleSize, double &Pmax, double &PWmax){
	long int h;
	long int Q=0;
	double PW;
	double P;
	uint64_t PWnum=0;

	for(long int m=0; m<=L; m+=h) {
		h=1;
		if(sa[m] + tupleSize <= L) {
			//This is a tuple to count
			while((m+h <= L) && ((saidx_t)tupleSize <= lcp[m+h])) {
				h++;
			}
			//h now contains the number of times this particular tuple occurred in the text

			//Account for the tuple
			//For t-tuple test
			if(Q < h) Q=h;

			//For the LRS test
			//The numerator must necessarily be even, so we can divide by 2 without problems
			//Check for overflows when adding to Psum element (unsigned 64 bit integers)
			PWnum += (((uint64_t)h)*((uint64_t)h-1)) >> 1;
		} 
	}

	//Now we're done with this tuple size. Finalize the outputs
	//t-tuple test
	if(Q >= 35) {
		//tupleSize <= t
		P = ((double)Q)/((double)(L-tupleSize+1));
		Pmax = pow(P, 1.0/(double)tupleSize);

		//tupleSize < u; set this value to flag nonsense
		PWmax = -1.0;
	} else {
		//Otherwise, tupleSize > t
		//a.k.a. tupleSize >= u
		PW = ((double)PWnum) / (double)(((L-tupleSize+1)*(L-tupleSize))>>1);
		PWmax = pow(PW, 1.0/((double)tupleSize));
		//tupleSize > t; set this value to flag nonsense
		Pmax = -1.0;
	}
}

void SAalgs(const byte text[], long int L, int k, double &t_tuple_res, double &lrs_res, bool verbose){
	vector <saidx_t> sa(L+1, -1);
	vector <saidx_t> lcp(L+1, -1);
   	long int u;
   	long int v;
   	long int lrs_len;
	double curPmax;
	double curPWmax;
	double Pmax;
	double PWmax;
	double pu;

	assert(L>0);
	assert(k>0);
	assert(L < INT32_MAX); //This is the default type, but it can be compiled to use 64 bit indexes (and then this should be INT64_MAX)

	calcSALCP(text, L, k, sa, lcp);

	//Find the length of the LRS, v
	lrs_len = -1;
	for(long int j=0; j<=L; j++) {
		if(lcp[j]>lrs_len) lrs_len = lcp[j];
	}

	assert(lrs_len>0);
	//lrs_len is now set correctly
   	v = lrs_len;

	//This is work factor on order O(L * v)
	Pmax = -1.0;
	PWmax = -1.0;
	u = 0;
	for(long int j=1; j<=v; j++) {
		countTuples(sa, lcp, L, j, curPmax, curPWmax);
		if(u == 0) {
			if(curPmax < 0.0) {
				u = j;
            			PWmax = curPWmax;
         		} else {
				// curPmax >= 0
				assert(curPWmax < 0.0);
				if(curPmax > Pmax) Pmax = curPmax;
			}
		} else {
			if(curPWmax > PWmax) PWmax = curPWmax;
		}
	}

	//finalize the t-tuple estimate
	if(Pmax > 0.0) {
		//We encountered a valid t, so we can run the test
		pu = Pmax + ZALPHA*sqrt(Pmax*(1.0 - Pmax)/((double)(L - 1)));
		if(verbose) printf("t-Tuple Estimate: t = %ld, p-hat_max = %.17g, p_u = %.17g\n", u-1, Pmax, pu);
		if(pu > 1.0) {
			pu = 1.0;
		}

		t_tuple_res = -log2(pu);
	} else {
		if(verbose) printf("t-Tuple Estimate: No strings are repeated 35 times. t-Tuple estimate failed.\n");
		t_tuple_res = -1.0;
	}

	//finalize the LRS estimate
	if(v>=u) {
		pu = PWmax + ZALPHA*sqrt(PWmax*(1.0 - PWmax)/((double)(L - 1)));
		if(pu > 1.0) {
			pu = 1.0;
		}
		if(verbose) printf("LRS Estimate: u = %ld, v = %ld, P_{max,W} = %.17g, p_u = %.17g\n", u, v, PWmax, pu);

		lrs_res = -log2(pu);
	} else {
		fprintf(stderr, "LRS Estimate: v<u. Can't Run LRS Test.\n");
		lrs_res = -1.0;
		return;
	}

	return;
}

int len_LRS(const byte text[], const int sample_size, int k){
	vector <saidx_t> sa(sample_size+1, -1);
	vector <saidx_t> lcp(sample_size+1, -1);
	saidx_t lrs_len = -1;

	calcSALCP(text, sample_size, k, sa, lcp);

	for(saidx_t j = 0; j <= sample_size; j++) {
		if(lcp[j] > lrs_len) lrs_len = lcp[j];
	}

	return(lrs_len);
}

void count_tuples(const byte data[], const int length, map<vector<byte>, int> &tuples, const int sample_size){

	for(int i = 0; i < sample_size-length; i++){
		vector<byte> substring = substr(data, i, length, sample_size);
		if(tuples.find(substring) == tuples.end()){
			tuples[substring] = 1;
		}else{
			tuples[substring]++;
		}
	}
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

	int lrs = len_LRS(data, sample_size, alphabet_size);
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
