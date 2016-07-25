#include "shared/utils.h"
#include "shared/most_common.h"


byte dataset[SIZE];

int main(){

	// Read in file
	const char* file_path = "../bin/truerand_8bit.bin";
	read_file(file_path, dataset);

	// The maximum min-entropy is -log2(1/SIZE)
	double min_entropy = log2(SIZE);

	// Estimate entropy with Most Common Value
	double H_min = most_common(dataset);
	
	#ifdef VERBOSE
	cout << "Most Common Value Estimate = " << H_min << endl;
	#endif

	min_entropy = min(min_entropy, H_min);

	// Estimate entropy with Collision Test

	// Estimate entropy with Markov Test

	// Estimate entropy with Compression Test

	// Estimate entropy with t-Tuple Test

	// Estimate entropy with Longest Repeated Substring Test (LRS)

	// Estimate entropy with Multi Most Common in Window Test

	// Estimate entropy with Lag Prediction Test

	// Estimate entropy with Multi MMC Test

	// Estimate entropy with LZ78Y Test
}