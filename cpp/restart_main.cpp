#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "non_iid/collision_test.h"
#include "non_iid/lz78y_test.h"
#include "non_iid/multi_mmc_test.h"
#include "non_iid/lag_test.h"
#include "non_iid/multi_mcw_test.h"
#include "non_iid/compression_test.h"
#include "non_iid/markov_test.h"

//Each test has a targeted chance of roughly 0.000005, and we need to witness at least 5 failures, so this should be no less than 1000000
#define SIMULATION_ROUNDS 5000000

void print_usage(){
	printf("Usage is: ea_restart <file_name> <bits_per_word> <H_I> <-i|-n> [-v]\n\n");
	printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (words).\n");
	printf("\t <bits_per_word>: Must be between 1-8, inclusive.\n");
	printf("\t <H_I>: Initial entropy estimate.\n");
	printf("\t <-i|-n>: '-i' for IID data, '-n' for non-IID data.\n");
	printf("\t -v: Optional verbosity flag for more output.\n");
	printf("\n");
	printf("\t Restart samples are assumed to be packed into 8-bit values, where the rightmost 'bits_per_word'\n");
	printf("\t bits constitute the sample. For example, if 'bits_per_word' is 3, then the four samples\n"); 
	printf("\t 0x6F, 0xA4, 0x39, 0x58, would be truncated to 0x07, 0x04, 0x01, 0x00.\n");
	printf("\n");
        printf("\t If there are less than 2^{bits_per_word} symbols observed in the data, the alphabet is\n");
        printf("\t mapped down to 0, 1, 2, ..., alph_size-1 in ascending numeric order of the symbols.\n");
        printf("\t For example, given 'bits_per_word' is 4, if the data consists of the three unique symbols\n");
        printf("\t 0x7, 0x3, 0xA, they would be mapped down to 0x3 => 0x0, 0x7 => 0x1, 0xA => 0x2.\n");
	printf("\n");
	printf("\t This program performs restart testing as described in Restart Tests (Section 3.1.4). The data\n"); 
	printf("\t consists of 1000 restarts, each with 1000 samples. The data is converted to rows and columns\n");
	printf("\t as described Section 3.1.4.1. The sanity check (Section 3.1.4.3) and the validation test\n");
	printf("\t (Section 3.1.4.2) are performed on this data.\n");
	printf("\n");
	printf("\t If the restart data passes the sanity check and validation test, this program returns\n");
	printf("\t min(H_r, H_c, H_I), which is either the validated entropy assessment or used to derive\n");
	printf("\t 'h_in' if conditioning is used (Section 3.1.5).\n"); 
	printf("\n");
}

//Here, we simulate a sort of "worst case" for this test, where there are a maximal number of symbols with maximal probability,
//and the rest is distributed to the other symbols
long int simulateCount(int k, double H_I, uint64_t *xoshiro256starstarState) {
	long int counts[256];
	int current_symbol;
	int k_max;
	long int max_count=0;
	double p, max_cutoff, p_min, cur_rand;

	assert(k<=256);

	p = pow(2.0, -H_I);

	k_max = floor(1.0/p);

	assert(k_max <= k);


	if(k>k_max) {
		max_cutoff = p * k_max;
		p_min = (1.0-max_cutoff)/(k-k_max);
	} else {
		max_cutoff = 1.0;
		p_min = 0.0;
	}

	for(int j=0; j<k; j++) counts[j] = 0;

	for(int j=0; j<1000; j++) {
		cur_rand = randomUnit(xoshiro256starstarState);
		if(cur_rand < max_cutoff) {
			current_symbol = floor(cur_rand / p);
			assert((current_symbol >= 0) && (current_symbol < k_max));
		} else {
			current_symbol = floor((cur_rand-max_cutoff) / p_min) + k_max;
			assert((current_symbol >= k_max) && (current_symbol < k));
		}
		counts[current_symbol]++;
	}

	for(int j=0; j<k; j++) {
		if(max_count < counts[j]) max_count = counts[j];
	}

	return max_count;
}

//This returns the bound (cutoff) for the test. Counts equal to this value should pass.
//Larger values should fail.
long int simulateBound(double alpha, int k, double H_I){
	uint64_t xoshiro256starstarMainSeed[4];
	uint64_t xoshiro256starstarSeed[4];
	vector<long int> results(SIMULATION_ROUNDS, -1);
	long int returnIndex;

	assert((k>1) && (k<=256));

	seed(xoshiro256starstarMainSeed);

        #pragma omp parallel private(xoshiro256starstarSeed)
	{
		memcpy(xoshiro256starstarSeed, xoshiro256starstarMainSeed, sizeof(xoshiro256starstarMainSeed));
		//Cause the RNG to jump omp_get_thread_num() * 2^128 calls
		xoshiro_jump(omp_get_thread_num(), xoshiro256starstarSeed);

		#pragma omp for
		for(int i = 0; i < SIMULATION_ROUNDS; i++){
			results[i] = simulateCount(k, H_I, xoshiro256starstarSeed);
		}
	}

	sort(results.begin(), results.end());
	assert((results[0]>=(1000/k)) && (results[0] <= 1000));
	assert((results[SIMULATION_ROUNDS-1]>=(1000/k)) && (results[SIMULATION_ROUNDS-1] <= 1000));
	assert(results[0] <= results[SIMULATION_ROUNDS-1]);

	returnIndex = ((size_t)floor((1.0 - alpha) * ((double)SIMULATION_ROUNDS))) - 1;
	assert((returnIndex >= 0) && (returnIndex < SIMULATION_ROUNDS));

	return(results[returnIndex]);
}

int main(int argc, char* argv[]){
	bool iid, verbose = false;
	const char verbose_flag = 'v';
	char *file_path;
	int r = 1000, c = 1000;
	int counts[256];
	long int X_cutoff;
	long i, j, X_i, X_r, X_c, X_max;
	double H_I, H_r, H_c, alpha, ret_min_entropy; 
	byte *rdata, *cdata;
	data_t data;

	// Parse args
	if(argc != 5 && argc != 6){
		printf("Incorrect usage.\n");
		print_usage();
		exit(-1);
	}
	else{
		// get filename
		file_path = argv[1];

		// get bits per word
		data.word_size = atoi(argv[2]);
		if(data.word_size < 1 || data.word_size > 8){
			printf("Invalid bits per word.\n");
			print_usage();
			exit(-1);
		}

		// get H_I	
		H_I = atof(argv[3]);
		if((H_I < 0) || (H_I > data.word_size)){
			printf("H_I must be nonnegative and at most 'bits_per_word'.\n");
			print_usage();
			exit(-1);
		}

		// get IID or non-IID
		if(argv[4][1] == 'i') iid = true;
		else if(argv[4][1] == 'n') iid = false;
		else{
			printf("Must specify whether data is IID or non-IID.\n");
			print_usage();
			exit(-1);
		}

		if(argc == 6) verbose = (argv[5][1] == verbose_flag);
	}

	if(verbose) printf("Opening file: '%s'\n", file_path);

	if(!read_file(file_path, &data)){
		printf("Error reading file.\n");
		print_usage();
		exit(-1);
	}

        if(data.alph_size == 1){
                printf("Symbol alphabet consists of 1 symbol. No entropy awarded...\n");
                free_data(&data);
                exit(-1);
        }

	if(data.len != MIN_SIZE){ 
		printf("\n*** Error: data does not contain %d samples ***\n\n", MIN_SIZE);
		exit(-1);
	}
	if(verbose) printf("Number of Symbols: %ld\n", data.len);
	if(verbose){
		if(data.alph_size < (1 << data.word_size)) printf("\nSymbols have been mapped down to an alphabet size of %d unique symbols\n\n", data.alph_size);
		else printf("\nSymbol alphabet consists of %d unique symbols\n\n", data.alph_size);
	}

	rdata = data.symbols;
	cdata = (byte*)malloc(data.len);
	if(cdata == NULL){
		printf("Error: failure to initialize memory for columns\n");
		exit(-1);
	}

	printf("H_I: %f\n", H_I);

	alpha = 1 - exp(log(0.99)/(r + c));
	X_cutoff = simulateBound(alpha, data.alph_size, H_I);
	printf("ALPHA: %.17g, X_cutoff: %ld\n", alpha, X_cutoff);

	// get maximum row count
	X_r = 0;
	for(i = 0; i < r; i++){ //row
		memset(counts, 0, 256*sizeof(int));
		X_i = 0;
		for(j = 0; j < c; j++){//column
			//[i*r+j] is row i, column j
			//So, we're fixing a row, and then iterate through various columns
			if(++counts[rdata[i*r+j]] > X_i) X_i = counts[rdata[i*r+j]];
		}
		if(X_i > X_r) X_r = X_i;
	}

	// construct column data from row data and get maximum column count
	X_c = 0;
	for(j = 0; j < c; j++){ //columns
		memset(counts, 0, 256*sizeof(int));
		X_i = 0;
		for(i = 0; i < r; i++){
			//[i*r+j] is row i, column j
			//So, we're fixing a column and iterating through various rows
			cdata[j*c+i] = rdata[i*r+j];
			if(++counts[cdata[j*c+i]] > X_i) X_i = counts[cdata[j*c+i]];
		}
		if(X_i > X_c) X_c = X_i;
	}

	// perform sanity check on rows and columns of restart data (Section 3.1.4.3)
	X_max = max(X_r, X_c);
	printf("X_max: %ld\n", X_max);
	if(X_max > X_cutoff){
		printf("\n*** Restart Sanity Check Failed ***\n");
		exit(-1);
	}
	else if(verbose) printf("\nRestart Sanity Check Passed...\n");

	// The maximum min-entropy is -log2(1/2^word_size) = word_size
	H_c = data.word_size;
	H_r = data.word_size;

	if(iid)	printf("\nRunning IID tests...\n\n");
	else printf("\nRunning non-IID tests...\n\n");

	printf("Running Most Common Value Estimate...\n");

	// Section 6.3.1 - Estimate entropy with Most Common Value
	ret_min_entropy = most_common(rdata, data.len, data.alph_size, verbose);
	if(verbose) printf("\tMost Common Value Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
	H_r = min(ret_min_entropy, H_r);
	ret_min_entropy = most_common(cdata, data.len, data.alph_size, verbose);
	if(verbose) printf("\tMost Common Value Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
	H_c = min(ret_min_entropy, H_c);

	if(!iid){
		if(data.word_size == 1){
			printf("\nRunning Entropic Statistic Estimates (bit strings only)...\n");

			// Section 6.3.2 - Estimate entropy with Collision Test (for bit strings only)
			ret_min_entropy = collision_test(rdata, data.len, verbose);
			if(verbose) printf("\tCollision Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy); 
			H_r = min(ret_min_entropy, H_r);
			ret_min_entropy = collision_test(cdata, data.len, verbose);
			if(verbose) printf("\tCollision Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy); 
			H_c = min(ret_min_entropy, H_c);

			// Section 6.3.3 - Estimate entropy with Markov Test (for bit strings only)
			ret_min_entropy = markov_test(rdata, data.len, verbose);
			if(verbose) printf("\tMarkov Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy); 
			H_r = min(ret_min_entropy, H_r);
			ret_min_entropy = markov_test(cdata, data.len, verbose);
			if(verbose) printf("\tMarkov Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy); 
			H_c = min(ret_min_entropy, H_c);

			// Section 6.3.4 - Estimate entropy with Compression Test (for bit strings only)
			ret_min_entropy = compression_test(rdata, data.len, verbose);
			if(ret_min_entropy >= 0){
				if(verbose) printf("\tCompression Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy); 
				H_r = min(ret_min_entropy, H_r);
			}
			ret_min_entropy = compression_test(cdata, data.len, verbose);
			if(ret_min_entropy >= 0){
				if(verbose) printf("\tCompression Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy); 
				H_c = min(ret_min_entropy, H_c);
			}
		}

		printf("\nRunning Tuple Estimates...\n");

		// Section 6.3.5 - Estimate entropy with t-Tuple Test
		double row_t_tuple_res, row_lrs_res;
		double col_t_tuple_res, col_lrs_res;
		SAalgs(rdata, data.len, data.alph_size, row_t_tuple_res, row_lrs_res, verbose);
		SAalgs(cdata, data.len, data.alph_size, col_t_tuple_res, col_lrs_res, verbose);

		if(verbose) printf("\tT-Tuple Test Estimate (Rows) = %f / %d bit(s)\n", row_t_tuple_res, data.word_size);
		H_r = min(row_t_tuple_res, H_r);

		if(verbose) printf("\tT-Tuple Test Estimate (Cols) = %f / %d bit(s)\n", col_t_tuple_res, data.word_size);
		H_c = min(col_t_tuple_res, H_c);

		// Section 6.3.6 - Estimate entropy with LRS Test
		if(verbose) printf("\tLRS Test Estimate (Rows) = %f / %d bit(s)\n", row_lrs_res, data.word_size);
		H_r = min(row_lrs_res, H_r);
		if(verbose) printf("\tLRS Test Estimate (Cols) = %f / %d bit(s)\n", col_lrs_res, data.word_size);
		H_c = min(col_lrs_res, H_c);

		printf("\nRunning Predictor Estimates...\n");

		// Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test
		ret_min_entropy = multi_mcw_test(rdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_r = min(ret_min_entropy, H_r);
		}
		ret_min_entropy = multi_mcw_test(cdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_c = min(ret_min_entropy, H_c);
		}

		// Section 6.3.8 - Estimate entropy with Lag Prediction Test
		ret_min_entropy = lag_test(rdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tLag Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_r = min(ret_min_entropy, H_r);
		}
		ret_min_entropy = lag_test(cdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tLag Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_c = min(ret_min_entropy, H_c);
		}

		// Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)
		ret_min_entropy = multi_mmc_test(rdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_r = min(ret_min_entropy, H_r);
		}
		ret_min_entropy = multi_mmc_test(cdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_c = min(ret_min_entropy, H_c);
		}

		// Section 6.3.10 - Estimate entropy with LZ78Y Test
		ret_min_entropy = LZ78Y_test(rdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tLZ78Y Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_r = min(ret_min_entropy, H_r);
		}
		ret_min_entropy = LZ78Y_test(cdata, data.len, data.alph_size, verbose);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tLZ78Y Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_c = min(ret_min_entropy, H_c);
		}
	}

	printf("\nH_r: %f\n", H_r);
	printf("H_c: %f\n", H_c);
	printf("H_I: %f\n\n", H_I);

	if(min(H_r, H_c) < H_I/2.0) printf("*** min(H_r, H_c) < H_I/2, Validation Testing Failed ***\n");
	else{
		printf("Validation Test Passed...\n\n");
		printf("min(H_r, H_c, H_I): %f\n\n", min(min(H_r, H_c), H_I));
	}

	free(cdata);
	free_data(&data);
}
