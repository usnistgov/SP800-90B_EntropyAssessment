#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "iid/permutation_tests.h"
#include "iid/chi_square_tests.h"


void print_usage(){
	cout << "Usage is: ./iid_main.out <file_name> <bits_per_word> [-v]" << endl;
	cout << "\t <file_name>: Must be relative path to a binary file with at least 1 million entries (words)." << endl;
	cout << "\t <bits_per_word>: Must be between 1-8, inclusive." << endl;
	cout << "\t -v: Optional verbosity flag for more output." << endl;
	cout << endl;
}

int main(int argc, char* argv[]){

	byte dataset[SIZE];
	bool verbose = false, is_binary;
	double mean, median;
	char* file_path;

	// Parse args
	if(argc != 3 && argc != 4){
		cout << "Incorrect usage." << endl;
		print_usage();
		exit(-1);
	}else{

		// Gather filename
		file_path = argv[1];

		// Gather bits per word
		const int word_size = atoi(argv[2]);
		if(word_size < 1 || word_size > 8){
			cout << "Invalid bits per word." << endl;
			print_usage();
			exit(-1);
		}

		const char verbose_flag = 'v';
		if(argc == 4){
			verbose = (argv[3][1] == verbose_flag);
		}
	}

	if(verbose){
		cout << "Opening file: " << file_path << endl;
	}

	if(!read_file(file_path, dataset)){
		cout << "Error reading file. Need 1 million entries for the tests to work." << endl;
		print_usage();
		exit(-1);
	}

	// Calculate baseline statistics
	cout << "Calculating baseline statistics..." << endl;
	calc_stats(dataset, mean, median, is_binary);

	if(verbose){
		cout << "\tMean: " << mean << endl;
		cout << "\tMedian: " << median << endl;
		cout << "\tBinary: " << (is_binary ? "true" : "false") << endl << endl;
	}

	// Compute permutation stats
	bool perm_test_pass = permutation_tests(dataset, mean, median, is_binary, verbose);

	if(perm_test_pass){
		cout << "** Passed IID permutation tests" << endl;
	}else{
		cout << "** Failed IID permutation tests" << endl;
		return -1;
	}

	// Compute chi square stats
	bool chi_square_test_pass = chi_square_tests(dataset, mean, median, is_binary, verbose);

	if(chi_square_test_pass){
		cout << "** Passed chi square tests" << endl;
	}else{
		cout << "** Failed chi square tests" << endl;
		return -1;
	}

	// Compute length of the longest repeated substring stats
	bool len_LRS_test_pass = len_LRS_test(dataset, verbose);

	if(len_LRS_test_pass){
		cout << "** Passed length of longest repeated substring test" << endl;
	}else{
		cout << "** Failed length of longest repeated substring test" << endl;
		return -1;
	}

	// Compute the min-entropy of the dataset
	double H_min = most_common(dataset);
	cout << "min-entropy = " << H_min << endl;

	return 0;
}
