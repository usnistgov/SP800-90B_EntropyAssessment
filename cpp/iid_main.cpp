#include "utils.h"
#include "permutation_tests.h"


using namespace std;

typedef unsigned char byte;

#define SIZE 1000000
#define PERMS 10

byte dataset[SIZE];

double mean = 0.0;
double median = 0.0;
bool is_binary = false;

int main(){

	// Read in file
	const char* file_path = "../bin/beacon_rand.bin";
	read_file(file_path);

	// Set variable
	for(int i = 0; i < SIZE; i++){
		data[i] = s[i];
	}

	// Calculate baseline statistics
	cout << "Calculating baseline statistics..." << endl;
	calc_stats();

	#ifdef VERBOSE
	cout << "Mean: " << mean << endl;
	cout << "Median: " << median << endl;
	cout << "Binary: " << (is_binary ? "true" : "false") << endl;
	#endif

	permutation_tests(s);

	return 0;
}