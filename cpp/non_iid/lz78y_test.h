#pragma once
#include "../shared/utils.h"

#define B 16
#define MAX_DICTIONARY_SIZE 65536

// Section 6.3.10 - LZ78Y Prediction Estimate
double LZ78Y_test(byte *data, long len, int alph_size, const bool verbose){
	int dict_size;
	long i, j, N, C, count, max_count, run_len, max_run_len;
	byte y, prediction;
	array<byte, B> prev;
	double p_global, p_local;
	bool found_prev, have_prediction;
	// j             prev          y          D[x,y]
	array<map<array<byte, B>, map<byte, long>>, B> D;

	if(len < B+2){	
		printf("\t*** Warning: not enough samples to run LZ78Y test (need more than %d) ***\n", B+2);
		return -1.0;
	}

	N = len-B-1;
	C = 0;
	run_len = 0;
	max_run_len = 0;

	// initialize dictionary counts
	dict_size = 0;
	memset(prev.data(), 0, B);
	for(j = 0; j < B; j++){
		memcpy(prev.data(), data+B-j-1, j+1);
		D[j][prev][data[B]] = 1;
		dict_size++;
	}

	// perform predictions
	for(i = B+1; i < len; i++){
		max_count = 0;
		have_prediction = false;
		memset(prev.data(), 0, B);
		for(j = 0; j < B; j++){
			// check if prev has been previously seen. If (prev,y) has not occurred, 
			// then do not make a prediction for current d and larger d's 
			// as well, since it will not occur for them either. In other words,
			// prediction is NULL, so do not update the scoreboard.
			if((j == 0) || found_prev){
				memcpy(prev.data(), data+i-j-1, j+1);
				if(D[j].find(prev) == D[j].end()) found_prev = false;
				else found_prev = true;
			}

			if(found_prev){
				y = max_map(D[j][prev]);
				count = D[j][prev][y];
				if(count >= max_count){
					max_count = count;
					prediction = y;
					have_prediction = true;
				}
				D[j][prev][data[i]]++;
			}
			else if(dict_size < MAX_DICTIONARY_SIZE){
				memcpy(prev.data(), data+i-j-1, j+1);
				D[j][prev][data[i]] = 1;
				dict_size++;
			}
		}
	
		// test	prediction of maximum (prev,y) pair
		if(have_prediction && (prediction == data[i])){
			C++;
			if(++run_len > max_run_len) max_run_len = run_len;
		}
		else run_len = 0;
	}

	p_global = calc_p_global(C, N);
	p_local = calc_p_local(max_run_len, N);

	if(verbose) printf("LZ78Y Prediction Estimate: N = %ld, Pglobal' = %.17g (C = %ld) Plocal = %.17g (r = %ld)\n", N, p_global, C, p_local, max_run_len+1);

	return -log2(max(max(p_global, p_local), 1/(double)alph_size));
}
