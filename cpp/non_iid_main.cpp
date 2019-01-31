#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "non_iid/collision_test.h"
#include "non_iid/lz78y_test.h"
#include "non_iid/multi_mmc_test.h"
#include "non_iid/lag_test.h"
#include "non_iid/multi_mcw_test.h"
#include "non_iid/tuple.h"
#include "non_iid/compression_test.h"
#include "non_iid/markov_test.h"


void print_usage(){
	printf("Usage is: ea_non_iid <file_name> <bits_per_word> <-i|-c> <-a|-t> [-v]\n\n");
	printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (words).\n");
	printf("\t <bits_per_word>: Must be between 1-8, inclusive.\n");
	printf("\t <-i|-c>: '-i' for initial entropy estimate, '-c' for conditioned sequential dataset entropy estimate.\n");
	printf("\t <-a|-t>: '-a' tests all bits in bitstring, '-t' truncates bitstring to %d bits.\n", MIN_SIZE);
	printf("\t          ('-a' is forced if bits_per_word is 1)\n");
	printf("\t -v: Optional verbosity flag for more output.\n");
	printf("\n");
	printf("\t Samples are assumed to be packed into 8-bit values, where the rightmost 'bits_per_word'\n");
	printf("\t bits constitute the sample. For example, if 'bits_per_word' is 3, then the four samples\n"); 
	printf("\t 0x6F, 0xA4, 0x39, 0x58, would be truncated to 0x07, 0x04, 0x01, 0x00.\n");
	printf("\n");
	printf("\t If there are less than 2^{bits_per_word} symbols observed in the data, the alphabet is\n");
	printf("\t mapped down to 0, 1, 2, ..., alph_size-1 in ascending numeric order of the symbols.\n");
	printf("\t For example, given 'bits_per_word' is 4, if the data consists of the three unique symbols\n");
	printf("\t 0x7, 0x3, 0xA, they would be mapped down to 0x3 => 0x0, 0x7 => 0x1, 0xA => 0x2.\n");
	printf("\n");
	printf("\t -i: Initial Entropy Estimate (Section 3.1.3)\n");
	printf("\n");
	printf("\t\t Computes the initial entropy estimate H_I as described in Section 3.1.3\n");
	printf("\t\t (not accounting for H_submitter) using the ten entropy estimators specified in\n");
	printf("\t\t Section 6.3.  If 'bits_per_word' is greater than 1, the samples are also\n");
	printf("\t\t converted to bitstrings. Two entropy estimates are computed: H_original and H_bitstring.\n");
	printf("\t\t Note that if 'bits_per_word' is 1, only H_bitstring is computed.\n"); 
	printf("\t\t Returns min(H_original, bits_per_word X H_bitstring). The initial entropy\n");
	printf("\t\t estimate H_I = min(H_submitter, H_original, bits_per_word X H_bitstring).\n");
	printf("\n");
	printf("\t -c: Conditioned Sequential Dataset Entropy Estimate (Section 3.1.5.2)\n");
	printf("\n");
	printf("\t\t Computes the entropy estimate per bit h' for the conditioned sequential dataset if the\n");
	printf("\t\t conditioning function is non-vetted. The samples are converted to a bitstring.\n");
	printf("\t\t Returns h' = min(H_bitstring).\n");
	printf("\n");
}

int main(int argc, char* argv[]){
	bool initial_entropy, all_bits, verbose = false;
	const char verbose_flag = 'v';
	char *file_path;
	long u, bu;
	double H_original, H_bitstring, ret_min_entropy; 
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

		// get initial entropy estimate or conditioned dataset flag
		if(argv[3][1] == 'i') initial_entropy = true;
		else if(argv[3][1] == 'c') initial_entropy = false;
		else{
			printf("Must specify whether computing initial entropy estimate or conditioned dataset.\n");
			print_usage();
			exit(-1);
		}

		// get bitstring length flag
		if(data.word_size > 1){
			if(argv[4][1] == 'a') all_bits = true;
			else if(argv[4][1] == 't') all_bits = false;
			else{
				printf("Must specify whether or not to truncate bitsring to %d bits.\n", MIN_SIZE);
				print_usage();
				exit(-1);
			}
		}
		else all_bits = true;

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

	if(!all_bits && (data.blen > MIN_SIZE)) data.blen = MIN_SIZE;

	if(verbose && initial_entropy && (data.word_size > 1)) printf("Number of Symbols: %ld\n", data.len);
	if(verbose) printf("Number of Binary Symbols: %ld\n", data.blen);
	if(data.len < MIN_SIZE) printf("\n*** Warning: data contains less than %d samples ***\n\n", MIN_SIZE);
	if(verbose){
		if(data.alph_size < (1 << data.word_size)) printf("\nSymbols have been mapped down to an alphabet size of %d unique symbols\n", data.alph_size);
		else printf("\nSymbol alphabet consists of %d unique symbols\n", data.alph_size);
	}

	// The maximum min-entropy is -log2(1/2^word_size) = word_size
	// The maximum bit string min-entropy is 1.0
	H_original = data.word_size;
	H_bitstring = 1.0;

	printf("\nRunning non-IID tests...\n\n");
	printf("Running Most Common Value Estimate...\n");

	// Section 6.3.1 - Estimate entropy with Most Common Value
	ret_min_entropy = most_common(data.bsymbols, data.blen, 2);
	if(verbose) printf("\tMost Common Value Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
	H_bitstring = min(ret_min_entropy, H_bitstring);
	if(initial_entropy && (data.word_size > 1)){
		ret_min_entropy = most_common(data.symbols, data.len, data.alph_size);
		if(verbose) printf("\tMost Common Value Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
		H_original = min(ret_min_entropy, H_original);
	}

	printf("\nRunning Entropic Statistic Estimates (bit strings only)...\n");

	// Section 6.3.2 - Estimate entropy with Collision Test (for bit strings only)
	ret_min_entropy = collision_test(data.bsymbols, data.blen);
	if(verbose) printf("\tCollision Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy); 
	H_bitstring = min(ret_min_entropy, H_bitstring);

	// Section 6.3.3 - Estimate entropy with Markov Test (for bit strings only)
	ret_min_entropy = markov_test(data.bsymbols, data.blen);
	if(verbose) printf("\tMarkov Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy); 
	H_bitstring = min(ret_min_entropy, H_bitstring);


	// Section 6.3.4 - Estimate entropy with Compression Test (for bit strings only)
	ret_min_entropy = compression_test(data.bsymbols, data.blen);
	if(ret_min_entropy >= 0){
		if(verbose) printf("\tCompression Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy); 
		H_bitstring = min(ret_min_entropy, H_bitstring);
	}

	printf("\nRunning Tuple Estimates...\n");

	// Section 6.3.5 - Estimate entropy with t-Tuple Test
	ret_min_entropy = t_tuple_test(data.bsymbols, data.blen, 2, &bu);
	if(verbose) printf("\tT-Tuple Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
	H_bitstring = min(ret_min_entropy, H_bitstring);
	if(initial_entropy && (data.word_size > 1)){
		ret_min_entropy = t_tuple_test(data.symbols, data.len, data.alph_size, &u);
		if(verbose) printf("\tT-Tuple Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
		H_original = min(ret_min_entropy, H_original);
	}

	// Section 6.3.6 - Estimate entropy with LRS Test
	ret_min_entropy = lrs_test(data.bsymbols, data.blen, 2, bu);
	if(verbose) printf("\tLRS Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
	H_bitstring = min(ret_min_entropy, H_bitstring);
	if(initial_entropy && (data.word_size > 1)){
		ret_min_entropy = lrs_test(data.symbols, data.len, data.alph_size, u);
		if(verbose) printf("\tLRS Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
		H_original = min(ret_min_entropy, H_original);
	}

	printf("\nRunning Predictor Estimates...\n");

	// Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test
	ret_min_entropy = multi_mcw_test(data.bsymbols, data.blen, 2);
	if(ret_min_entropy >= 0){
		if(verbose) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
		H_bitstring = min(ret_min_entropy, H_bitstring);
	}
	if(initial_entropy && (data.word_size > 1)){
		ret_min_entropy = multi_mcw_test(data.symbols, data.len, data.alph_size);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_original = min(ret_min_entropy, H_original);
		}
	}

	// Section 6.3.8 - Estimate entropy with Lag Prediction Test
	ret_min_entropy = lag_test(data.bsymbols, data.blen, 2);
	if(ret_min_entropy >= 0){
		if(verbose) printf("\tLag Prediction Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
		H_bitstring = min(ret_min_entropy, H_bitstring);
	}
	if(initial_entropy && (data.word_size > 1)){
		ret_min_entropy = lag_test(data.symbols, data.len, data.alph_size);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tLag Prediction Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_original = min(ret_min_entropy, H_original);
		}
	}

	// Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)
	ret_min_entropy = multi_mmc_test(data.bsymbols, data.blen, 2);
	if(ret_min_entropy >= 0){
		if(verbose) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
		H_bitstring = min(ret_min_entropy, H_bitstring);
	}
	if(initial_entropy && (data.word_size > 1)){
		ret_min_entropy = multi_mmc_test(data.symbols, data.len, data.alph_size);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_original = min(ret_min_entropy, H_original);
		}
	}

	// Section 6.3.10 - Estimate entropy with LZ78Y Test
	ret_min_entropy = LZ78Y_test(data.bsymbols, data.blen, 2);
	if(ret_min_entropy >= 0){
		if(verbose) printf("\tLZ78Y Prediction Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
		H_bitstring = min(ret_min_entropy, H_bitstring);
	}
	if(initial_entropy && (data.word_size > 1)){
		ret_min_entropy = LZ78Y_test(data.symbols, data.len, data.alph_size);
		if(ret_min_entropy >= 0){
			if(verbose) printf("\tLZ78Y Prediction Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
			H_original = min(ret_min_entropy, H_original);
		}
	}

	if(initial_entropy){
		printf("\nH_original: %f\n", H_original);
		printf("H_bitstring: %f\n\n", H_bitstring);
		printf("min(H_original, %d X H_bitstring): %f\n\n", data.word_size, min(H_original, data.word_size*H_bitstring));
	}
	else printf("\nh': %f\n", H_bitstring);

	free_data(&data);
}
