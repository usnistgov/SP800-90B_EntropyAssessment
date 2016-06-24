#include "utils.h"
#include "permutation_tests.h"


byte dataset[SIZE];

double mean = 0.0;
double median = 0.0;
bool is_binary = false;

int main(){

	// Read in file
	const char* file_path = "../bin/beacon_rand.bin";
	read_file(file_path, dataset);

	// Calculate baseline statistics
	cout << "Calculating baseline statistics..." << endl;
	calc_stats(dataset, mean, median, is_binary);

	#ifdef VERBOSE
	cout << "Mean: " << mean << endl;
	cout << "Median: " << median << endl;
	cout << "Binary: " << (is_binary ? "true" : "false") << endl;
	#endif

	permutation_tests(dataset, mean, median, is_binary);

	return 0;
}