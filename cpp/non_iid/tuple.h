#pragma once
#include "../shared/utils.h"

// Section 6.3.5 - t-Tuple Estimate
double t_tuple_test(byte *data, long len, int alph_size, long *u_p){
	long i, pos, c, tup_size, Q, thresh = 35;
	double t_tuple_min_entropy, min_entropy, pmax;
	vector<byte> tuple;
	map<vector<byte>, vector<long>> tuples, prev_tuples;
	map<vector<byte>, vector<long>>::iterator itr;
	
	t_tuple_min_entropy = log2(alph_size);

	// find smallest tuple occurrences
	Q = 0;
	tup_size = 1;
	tuple.push_back(0); // allocate space for tuple
	for(pos = 0; pos < len; pos++){ 
		memcpy(tuple.data(), data+pos, tup_size);
		tuples[tuple].push_back(pos);
		if(tuples[tuple].size() > Q) Q = tuples[tuple].size();
	}

	while(Q >= thresh){
		min_entropy = -log2(Q/(double)(len-tup_size+1)) / tup_size;
		if(min_entropy <= t_tuple_min_entropy) t_tuple_min_entropy = min_entropy; 

		Q = 0;
		tuples.swap(prev_tuples);
		tup_size++;
		tuple.push_back(0); // allocate additional space for tuple
		// find larger tuple occurrences from smaller ones
		for(itr = prev_tuples.begin(); itr != prev_tuples.end(); itr++){
			c = (itr->second).size();	
			if(c >= thresh){
				for(i = 0; i < c; i++){
					pos = (itr->second)[i];
					if(pos+tup_size-1 < len){
						memcpy(tuple.data(), data+pos, tup_size);
						tuples[tuple].push_back(pos);
						if(tuples[tuple].size() > Q) Q = tuples[tuple].size();
					}
				}
			}
		}
		prev_tuples.clear(); // destroy previous tuples
	}

	*u_p = tup_size; // set start tuple size 'u' for LRS test

	pmax = pow(2.0, -t_tuple_min_entropy);
	return -log2(min(1.0, pmax + 2.576*sqrt(pmax*(1.0-pmax)/(len-1.0))));
}

// Section 6.3.6 - Longest Repeated Substring (LRS) Estimate
double lrs_test(byte *data, long len, int alph_size, long u){
	long i, pos, c, tup_size, Q;
	double lrs_min_entropy, min_entropy, pmax, sum;
	vector<byte> tuple;
	map<vector<byte>, vector<long>> tuples, prev_tuples;
	map<vector<byte>, vector<long>>::iterator itr;

	// Since 1/(n*n) < sum(C(C-1))/[n(n-1)], 
	// it follows that LRS min-entropy < 2*log2(n)
	lrs_min_entropy = 2*log2(alph_size);

	// find smallest tuple occurrences
	Q = 0;
	if(u > 0) tup_size = u;
	else tup_size = 1;
	for(i = 0; i < tup_size; i++) tuple.push_back(0); // allocate space for tuple
	for(pos = 0; pos < len; pos++){ 
		if(pos+tup_size-1 < len){
			memcpy(tuple.data(), data+pos, tup_size);
			tuples[tuple].push_back(pos);
			if(tuples[tuple].size() > Q) Q = tuples[tuple].size();
		}
	}

	while(Q >= 2){
		sum = 0;
		for(itr = tuples.begin(); itr != tuples.end(); itr++){
			c = (itr->second).size();	
			if(c >= 2) sum += c*(c-1);
		}
		min_entropy = (-log2(sum)+log2((double)len-tup_size+1)+log2((double)len-tup_size)) / tup_size;
		if(min_entropy <= lrs_min_entropy) lrs_min_entropy = min_entropy; 

		Q = 0;
		tuples.swap(prev_tuples);
		tup_size++;
		tuple.push_back(0); // allocate additional space for tuple
		// find larger tuple occurrences from smaller ones
		for(itr = prev_tuples.begin(); itr != prev_tuples.end(); itr++){
			c = (itr->second).size();	
			if(c >= 2){
				for(i = 0; i < c; i++){
					pos = (itr->second)[i];
					if(pos+tup_size-1 < len){
						memcpy(tuple.data(), data+pos, tup_size);
						tuples[tuple].push_back(pos);
						if(tuples[tuple].size() > Q) Q = tuples[tuple].size();
					}
				}
			}
		}
		prev_tuples.clear(); // destroy previous tuples
	}

	pmax = pow(2.0, -lrs_min_entropy);
	return -log2(min(1.0, pmax + 2.576*sqrt(pmax*(1.0-pmax)/(len-1.0))));
}
