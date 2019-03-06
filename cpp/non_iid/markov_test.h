#pragma once
#include "../shared/utils.h"

// Section 6.3.3 - Markov Estimate
// data is assumed to be binary (e.g., bit string)
double markov_test(byte* data, long len, const bool verbose){
	long i, C_0, C_1, C_00, C_01, C_10, C_11;
	double H_min, tmp_min_entropy, P_0, P_1, P_00, P_01, P_10, P_11;

	C_0 = 0.0;
	C_00 = 0.0;
	C_10 = 0.0;

	// get counts for unconditional and transition probabilities
	for(i = 0; i < len-1; i++){
		if(data[i] == 0){
			C_0++;
			if(data[i+1] == 0) C_00++;
		}
		else if(data[i+1] == 0) C_10++;
	}

	C_1 = len - 1 - C_0;
	C_01 = C_0 - C_00;
	C_11 = C_1 - C_10;

	P_00 = C_00 / (double)(C_00 + C_01);
	P_10 = C_10 / (double)(C_10 + C_11);
	P_01 = 1.0 - P_00;
	P_11 = 1.0 - P_10;

	// account for last symbol
	if(data[len-1] == 0) C_0++;

	P_0 = C_0 / (double)len;
	P_1 = 1.0 - P_0;

	if(verbose) printf("Markov Estimate: P_0 = %.17g, P_1 = %.17g, P_0,0 = %.17g, P_0,1 = %.17g, P_1,0 = %.17g, P_1,1 = %.17g, ", P_0, P_1, P_00, P_01, P_10, P_11);

	H_min = 128.0;

	// Sequence 00...0
	if(P_00 > 0){
		tmp_min_entropy = -log2(P_0) - 127*log2(P_00); 
		if(tmp_min_entropy < H_min) H_min = tmp_min_entropy;
	}

	// Sequence 0101...01
	if((P_01 > 0) && (P_10 > 0)){
		tmp_min_entropy = -log2(P_0) - 64*log2(P_01) - 63*log2(P_10);
		if(tmp_min_entropy < H_min) H_min = tmp_min_entropy;
	}

        // Sequence 011...1
	if((P_01 > 0) && (P_11 > 0)){
        	tmp_min_entropy = -log2(P_0) - log2(P_01) - 126*log2(P_11);
       		if(tmp_min_entropy < H_min) H_min = tmp_min_entropy;	
	}

        // Sequence 100...0
	if((P_10 > 0) && (P_00 > 0)){
        	tmp_min_entropy = -log2(P_1) - log2(P_10) - 126*log2(P_00);
       		if(tmp_min_entropy < H_min) H_min = tmp_min_entropy;
	}

        // Sequence 1010...10
	if((P_10 > 0) && (P_01 > 0)){
       		tmp_min_entropy = -log2(P_1) - 64*log2(P_10) - 63*log2(P_01);
        	if(tmp_min_entropy < H_min) H_min = tmp_min_entropy;
	}

        // Sequence 11...1
	if(P_11 > 0){
        	tmp_min_entropy = -log2(P_1) - 127*log2(P_11);
       		if(tmp_min_entropy < H_min) H_min = tmp_min_entropy;
	}

	if(verbose) printf("p_max = %.17g\n", pow(2.0, -H_min));

	return fmin(H_min/128.0, 1.0);
}
