#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "iid/permutation_tests.h"
#include "iid/chi_square_tests.h"
#include <omp.h>
#include <getopt.h>
#include <limits.h>



[[ noreturn ]] void print_usage(){
	printf("Usage is: ea_iid [-i|-c] [-a|-t] [-v] [-l <index>,<samples> ] <file_name> [bits_per_symbol]\n\n");
	printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (samples).\n");
	printf("\t [bits_per_symbol]: Must be between 1-8, inclusive. By default this value is inferred from the data.\n");
	printf("\t [-i|-c]: '-i' for initial entropy estimate, '-c' for conditioned sequential dataset entropy estimate. The initial entropy estimate is the default.\n");
	printf("\t [-a|-t]: '-a' produces the 'H_bitstring' assessment using all read bits, '-t' truncates the bitstring used to produce the `H_bitstring` assessment to %d bits. Test all data by default.\n", MIN_SIZE);
	printf("\t Note: When testing binary data, no `H_bitstring` assessment is produced, so the `-a` and `-t` options produce the same results for the initial assessment of binary data.\n");
	printf("\t -v: Optional verbosity flag for more output. Can be used multiple times.\n");
	printf("\t -l <index>,<samples>\tRead the <index> substring of length <samples>.\n");
	printf("\n");
	printf("\t Samples are assumed to be packed into 8-bit values, where the least significant 'bits_per_symbol'\n");
	printf("\t bits constitute the symbol.\n");
	printf("\n");
	printf("\t -i: Initial Entropy Estimate (Section 3.1.3)\n");
	printf("\n");
	printf("\t\t Computes the initial entropy estimate H_I as described in Section 3.1.3\n");
	printf("\t\t (not accounting for H_submitter) using the entropy estimators specified in\n");
	printf("\t\t Section 6.3.  If 'bits_per_symbol' is greater than 1, the samples are also\n");
	printf("\t\t converted to bitstrings and assessed to create H_bitstring; for multi-bit symbols,\n");
	printf("\t\t two entropy estimates are computed: H_original and H_bitstring.\n");
	printf("\t\t Returns min(H_original, bits_per_symbol X H_bitstring). The initial entropy\n");
	printf("\t\t estimate H_I = min(H_submitter, H_original, bits_per_symbol X H_bitstring).\n");
	printf("\n");
	printf("\t -c: Conditioned Sequential Dataset Entropy Estimate (Section 3.1.5.2)\n");
	printf("\n");
	printf("\t\t Computes the entropy estimate per bit h' for the conditioned sequential dataset if the\n");
	printf("\t\t conditioning function is non-vetted. The samples are converted to a bitstring.\n");
	printf("\t\t Returns h' = min(H_bitstring).\n");
	printf("\n");
	exit(-1);
}

int main(int argc, char* argv[]){

	bool initial_entropy, all_bits;
	int verbose = 0;
	double rawmean, median;
	char* file_path;
	data_t data;
	int opt;
	unsigned long subsetIndex = ULONG_MAX;
	unsigned long subsetSize = 0;
	unsigned long long inint;
	char *nextOption;

	data.word_size = 0;
	initial_entropy = true;
	all_bits = true;

	while ((opt = getopt(argc, argv, "icatvl:")) != -1) {
		switch(opt) {
			case 'i':
				initial_entropy = true;
				break;
			case 'c':
				initial_entropy = false;
				break;
			case 'a':
				all_bits = true;
				break;
			case 't':
				all_bits = false;
				break;
			case 'v':
				verbose++;
				break;
                        case 'l':
                                inint = strtoull(optarg, &nextOption, 0);
                                if((inint > ULONG_MAX) || (errno == EINVAL) || (nextOption == NULL) || (*nextOption != ',')) {
                                        print_usage();
                                }
                                subsetIndex = inint;

                                nextOption++;

                                inint = strtoull(nextOption, NULL, 0);
                                if((inint > ULONG_MAX) || (errno == EINVAL)) {
                                        print_usage();
                                }
                                subsetSize = inint;
                                break;
			default:
				print_usage();
		}
	}

	argc -= optind;
	argv += optind;

	// Parse args
	if((argc != 2) && (argc != 1)){
		printf("Incorrect usage.\n");
		print_usage();
	}

	// get filename
	file_path = argv[0];

	if(argc == 2) {
		// get bits per word
		data.word_size = atoi(argv[1]);
		if(data.word_size < 1 || data.word_size > 8){
			printf("Invalid bits per symbol.\n");
			print_usage();
		}
	}

	if(verbose > 0){
		printf("Opening file: '%s'\n", file_path);
	}

	if(!read_file_subset(file_path, &data, subsetIndex, subsetSize)){
		printf("Error reading file.\n");
		print_usage();
	}
	if(verbose > 0) printf("Loaded %ld samples of %d distinct %d-bit-wide symbols\n", data.len, data.alph_size, data.word_size);

	if(data.alph_size <= 1){
		printf("Symbol alphabet consists of 1 symbol. No entropy awarded...\n");
		free_data(&data);
		exit(-1);
	}

	if(!all_bits && (data.blen > MIN_SIZE)) data.blen = MIN_SIZE;

	if((verbose > 0) && ((data.alph_size > 2) || !initial_entropy)) printf("Number of Binary samples: %ld\n", data.blen);
	if(data.len < MIN_SIZE) printf("\n*** Warning: data contains less than %d samples ***\n\n", MIN_SIZE);
	if(verbose > 0){
		if(data.alph_size < (1 << data.word_size)) printf("\nSamples have been translated\n");
	}

	// Calculate baseline statistics
	int alphabet_size = data.alph_size;
	int sample_size = data.len;

	printf("Calculating baseline statistics...\n");
	calc_stats(&data, rawmean, median);

	if(verbose > 0){
		printf("\tRaw Mean: %f\n", rawmean);
		printf("\tMedian: %f\n", median);
		printf("\tBinary: %s\n\n", (alphabet_size == 2 ? "true" : "false"));
	}

	double H_original = data.word_size;
	double H_bitstring = 1.0;

	// Compute the min-entropy of the dataset
	if(initial_entropy) {
		H_original = most_common(data.symbols, sample_size, alphabet_size, verbose, "Literal");
	}

	if(((data.alph_size > 2) || !initial_entropy)) {
		H_bitstring = most_common(data.bsymbols, data.blen, 2, verbose, "Bitstring");
	}

        if(verbose <= 1) {
                if(initial_entropy){
                        printf("H_original: %f\n", H_original);
                        if(data.alph_size > 2) {
                                printf("H_bitstring: %f\n", H_bitstring);
                                printf("min(H_original, %d X H_bitstring): %f\n", data.word_size, min(H_original, data.word_size*H_bitstring));
                        }
                } else printf("h': %f\n", H_bitstring);
        } else {
                double h_assessed = data.word_size;

                if((data.alph_size > 2) || !initial_entropy) {
                        h_assessed = min(h_assessed, H_bitstring * data.word_size);
                        printf("H_bitstring = %.17g\n", H_bitstring);
                }

                if(initial_entropy) {
                        h_assessed = min(h_assessed, H_original);
                        printf("H_original: %.17g\n", H_original);
                }

                printf("Assessed min entropy: %.17g\n", h_assessed);
        }

	printf("\n");

	// Compute chi square stats
	bool chi_square_test_pass = chi_square_tests(data.symbols, sample_size, alphabet_size, verbose);

	if(chi_square_test_pass){
		printf("** Passed chi square tests\n\n");
	}else{
		printf("** Failed chi square tests\n\n");
	}

	// Compute length of the longest repeated substring stats
	bool len_LRS_test_pass = len_LRS_test(data.symbols, sample_size, alphabet_size, verbose, "Literal");

	if(len_LRS_test_pass){
		printf("** Passed length of longest repeated substring test\n\n");
	}else{
		printf("** Failed length of longest repeated substring test\n\n");
	}

	// Compute permutation stats
	bool perm_test_pass = permutation_tests(&data, rawmean, median, verbose);

	if(perm_test_pass){
		printf("** Passed IID permutation tests\n\n");
	}else{
		printf("** Failed IID permutation tests\n\n");
	}

	free_data(&data);
	return 0;
}
