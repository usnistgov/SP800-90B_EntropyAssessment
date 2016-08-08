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


byte dataset[SIZE];

int main(){

	// Read in file
	const char* file_path = "../bin/truerand_8bit.bin";
	read_file(file_path, dataset);

	// The maximum min-entropy is -log2(1/SIZE)
	double min_entropy = log2(SIZE);
	double H_min;

	#ifdef VERBOSE
	cout << endl << "Running non-IID tests..." << endl;
	#endif

	// Section 6.3.1 - Estimate entropy with Most Common Value
	// H_min = most_common(dataset);
	// #ifdef VERBOSE
	// cout << "Most Common Value Estimate = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);

	// Section 6.3.2 - Estimate entropy with Collision Test
	// H_min = collision_test(dataset, 8);
	// #ifdef VERBOSE
	// cout << "Collision Test Estimate = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);

	// Section 6.3.3 - Estimate entropy with Markov Test

	// Section 6.3.4 - Estimate entropy with Compression Test
	H_min = compression_test(dataset);
	#ifdef VERBOSE
	cout << "Compression Test Estimate = " << H_min << endl;
	#endif
	min_entropy = min(min_entropy, H_min);

	// Section 6.3.5 - Estimate entropy with t-Tuple Test
	// H_min = t_tuple_test(dataset);
	// #ifdef VERBOSE
	// cout << "t-Tuple Test Estimate = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);

	// Section 6.3.6 - Estimate entropy with Longest Repeated Substring Test (LRS)
	// Slow, needs speedup, wrap LRS functions into a single call
	// H_min = LRS_test_noniid(dataset);
	// #ifdef VERBOSE
	// cout << "Longest Reapeated Substring Test Estimate = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);

	// Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test
	// H_min = multi_mcw_test(dataset);
	// #ifdef VERBOSE
	// cout << "Multi Most Common in Window (Multi MCW) Test = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);

	// Section 6.3.8 - Estimate entropy with Lag Prediction Test
	// H_min = lag_test(dataset);
	// #ifdef VERBOSE
	// cout << "Lag Prediction Test = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);

	// Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)
	// H_min = multi_mmc_test(dataset);
	// #ifdef VERBOSE
	// cout << "Multi Markov Model with Counting (MultiMMC) Prediction Test = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);

	// Section 6.3.10 - Estimate entropy with LZ78Y Test
	// Not super fast, just a touch longer than the python but with -O3 is super fast
	// H_min = LZ78Y_test(dataset);
	// #ifdef VERBOSE
	// cout << "LZ78Y Prediction Test = " << H_min << endl;
	// #endif
	// min_entropy = min(min_entropy, H_min);
}