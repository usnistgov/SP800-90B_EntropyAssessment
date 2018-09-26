#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "iid/permutation_tests.h"
#include "iid/chi_square_tests.h"

void print_usage(){
	printf("Usage is: ea_iid <file_name> <bits_per_word> <-i|-c> <-a|-t> [-v]\n\n");
	printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (words).\n");
	printf("\t <bits_per_word>: Must be between 1-8, inclusive.\n");
	printf("\t <-i|-c>: '-i' for initial entropy estimate, '-c' for conditioned sequential dataset entropy estimate.\n");
	printf("\t <-a|-t>: '-a' tests all bits in bitstring, '-t' truncates bitstring to %d bits.\n", MIN_SIZE);
	printf("\t\t ('-a' is forced if bits_per_word is 1)\n");
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
	double rawmean, median;
	char* file_path;
	int num_threads = 4;
	data_t data;

	// Parse args
	if(argc < 6){
		printf("Incorrect usage.\n");
		print_usage();
		exit(-1);
	}else{

		// Gather filename
		file_path = argv[1];

		// Gather bits per word
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
				printf("Must specify whether or not to truncate bitstring to %d bits.\n", MIN_SIZE);
				print_usage();
				exit(-1);
			}
		}
		else all_bits = true;

		if(argc == 6) verbose = (argv[5][1] == verbose_flag);
	}

	if(verbose){
		printf("Opening file: '%s'\n", file_path);
	}

	//byte* dataset = new byte[sample_size];
	//if(!read_file(file_path, dataset, word_size, sample_size)){
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

	// Calculate baseline statistics
	//alphabet_size = pow(2, word_size);
	int alphabet_size = data.alph_size;
	int sample_size = data.len;

	printf("Calculating baseline statistics...\n");
	calc_stats(&data, rawmean, median);

	if(verbose){
		printf("\tRaw Mean: %f\n", rawmean);
		printf("\tMedian: %f\n", median);
		printf("\tBinary: %s\n\n", (alphabet_size == 2 ? "true" : "false"));
	}

	double start_time = omp_get_wtime();

	// Compute the min-entropy of the dataset
	double H_min = most_common(data.symbols, sample_size, alphabet_size, verbose);
	printf("min-entropy = %f\n\n", H_min);

	// Compute chi square stats
	bool chi_square_test_pass = chi_square_tests(data.symbols, sample_size, alphabet_size, verbose);

	if(chi_square_test_pass){
		printf("** Passed chi square tests\n\n");
	}else{
		printf("** Failed chi square tests\n\n");
		//return -1;
	}

	// Compute length of the longest repeated substring stats
	bool len_LRS_test_pass = len_LRS_test(data.symbols, sample_size, alphabet_size, verbose);

	if(len_LRS_test_pass){
		printf("** Passed length of longest repeated substring test\n\n");
	}else{
		printf("** Failed length of longest repeated substring test\n\n");
		//return -1;
	}

	// Compute permutation stats
	bool perm_test_pass = permutation_tests(&data, rawmean, median, num_threads, verbose);

	if(perm_test_pass){
		printf("** Passed IID permutation tests\n\n");
	}else{
		printf("** Failed IID permutation tests\n\n");
		//return -1;
	}

	double run_time = omp_get_wtime() - start_time;
	printf("Total Time elapsed: %f\n", run_time);

	free_data(&data);
	return 0;
}
