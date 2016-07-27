#include "shared/utils.h"
#include "shared/most_common.h"
#include "non_iid/collision_test.h"


byte dataset[SIZE];

int main(){

	// Read in file
	const char* file_path = "../bin/truerand_8bit.bin";
	read_file(file_path, dataset);

	// The maximum min-entropy is -log2(1/SIZE)
	double min_entropy = log2(SIZE);

	// Section 6.3.1 - Estimate entropy with Most Common Value
	double H_min = most_common(dataset);
	
	#ifdef VERBOSE
	cout << "Most Common Value Estimate = " << H_min << endl;
	#endif

	min_entropy = min(min_entropy, H_min);

	// Section 6.3.2 - Estimate entropy with Collision Test
	// involves remapping to smaller word size
	H_min = collision_test(dataset);

	#ifdef VERBOSE
	cout << "Collision Test Estimate = " << H_min << endl;
	#endif

	min_entropy = min(min_entropy, H_min);

	// Section 6.3.3 - Estimate entropy with Markov Test

	// Section 6.3.4 - Estimate entropy with Compression Test

	// Section 6.3.5 - Estimate entropy with t-Tuple Test

	// Section 6.3.6 - Estimate entropy with Longest Repeated Substring Test (LRS)

	// Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test

	// Section 6.3.8 - Estimate entropy with Lag Prediction Test

	// Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)

	// Section 6.3.10 - Estimate entropy with LZ78Y Test
}