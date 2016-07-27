#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "iid/permutation_tests.h"
#include "iid/chi_square_tests.h"


byte dataset[SIZE];

double mean = 0.0;
double median = 0.0;
bool is_binary = false;

int main(){

	// Read in file
	const char* file_path = "../bin/truerand_8bit.bin";
	read_file(file_path, dataset);

	// Calculate baseline statistics
	cout << "Calculating baseline statistics..." << endl;
	calc_stats(dataset, mean, median, is_binary);

	#ifdef VERBOSE
	cout << "Mean: " << mean << endl;
	cout << "Median: " << median << endl;
	cout << "Binary: " << (is_binary ? "true" : "false") << endl;
	#endif

	// Compute permutation stats
	bool perm_test_pass = permutation_tests(dataset, mean, median, is_binary);

	if(perm_test_pass){
		cout << "** Passed IID permutation tests" << endl;
	}else{
		cout << "** Failed IID permutation tests" << endl;
		return 0;
	}

	// Compute chi square stats
	bool chi_square_test_pass = chi_square_tests(dataset, mean, median, is_binary);

	if(chi_square_test_pass){
		cout << "** Passed chi square tests" << endl;
	}else{
		cout << "** Failed chi square tests" << endl;
		return 0;
	}

	// Compute length of the longest repeated substring stats
	bool len_LRS_test_pass = len_LRS_test(dataset);

	if(len_LRS_test_pass){
		cout << "** Passed length of longest repeated substring test" << endl;
	}else{
		cout << "** Failed length of longest repeated substring test" << endl;
		return 0;
	}

	// Compute the min-entropy of the dataset
	double H_min = most_common(dataset);
	cout << "min-entropy = " << H_min << endl;

	return 0;
}