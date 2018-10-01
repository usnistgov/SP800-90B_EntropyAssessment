#pragma once

#include <stdlib.h>
#include <bzlib.h> // sudo apt-get install libbz2-dev
#include "../shared/utils.h"

// The tests used
const unsigned int num_tests = 19;
const string test_names[] = {"excursion","numDirectionalRuns","lenDirectionalRuns","numIncreasesDecreases","numRunsMedian","lenRunsMedian","avgCollision","maxCollision","periodicity(1)","periodicity(2)","periodicity(8)","periodicity(16)","periodicity(32)","covariance(1)","covariance(2)","covariance(8)","covariance(16)","covariance(32)","compression"};

using namespace std;

/*
* ---------------------------------------------
* 	  TASKS FOR PERMUTATION TESTS
* ---------------------------------------------
*/

// 5.1 Conversion I
// Takes a binary sequence and partitions it into 8-bit blocks
// Blocks have the number of 1's counted and totaled
//
// Requires binary data
vector<byte> conversion1(const byte data[], const int sample_size){
	vector<byte> ret((sample_size / 8) + 1, 0);

	for(unsigned int i = 0; i < sample_size; ++i){
		ret[i/8] += data[i];	// integer division to ensure the size of ret is sample_size / 8
	}

	return ret;
}

// 5.1 Conversion II
// Takes a binary sequence and partitions it into 8-bit blocks
// Blocks are then converted to decimal
//
// Requires binary data
vector<byte> conversion2(const byte data[], const int sample_size){
	vector<byte> ret((sample_size / 8) + 1, 0);

	for(unsigned int i = 0; i < sample_size; ++i){
		ret[i/8] += data[i] << (8 - ((i+1)%8));
	}

	return ret;
}

// 5.1.1 Excursion Test
// Measures how far the running sum of values deviates from the
// average value at each point in the set
//
// Requires binary or non-binary data
double excursion(const byte data[], const double mean, const int sample_size){
	double d_i = 0;
	double max = 0;
	double running_sum = 0;

	for(unsigned int i = 0; i < sample_size; ++i){
		running_sum += data[i];
		d_i = abs(running_sum - ((i+1) * mean));

		if(d_i > max){
			max = d_i;
		}
	}

	return max;
}

// Helper for 5.1.2, 5.1.3, and 5.1.4
// Builds a array of the runs of consecutive values
// Pushes -1 to the array if the value is > than the next
// Pushes +1 to the array if the value is <= than the next
//
// Requires non-binary data
vector<int> alt_sequence1(const byte data[], const int sample_size){
	vector<int> ret(sample_size-1, 0);

	for(unsigned int i = 0; i < sample_size-1; ++i){
		ret[i] = ((data[i] > data[i+1]) ? -1 : 1);
	}

	return ret;
}

// Helper for 5.1.5 and 5.1.6
// Builds a array of the runs of values compared to the median
// Pushes +1 to the array if the value is >= the median
// Pushes -1 to the array if the value is < than the median
vector<int> alt_sequence2(const byte data[], const double median, const int sample_size){
	vector<int> ret(sample_size, 0);

	for(unsigned int i = 0; i < sample_size; ++i){
		ret[i] = ((data[i] < median) ? -1 : 1);
	}

	return ret;
}

// 5.1.2 Number of Directional Runs
// Determines the number of runs in the sequence.
// A run is when multiple consecutive values are all >= the prior
// or all < the prior
//
// Requires data from alt_sequence1, binary data needs conversion1 first
//
//
// 5.1.5 Number of Runs Based on the Median
// Determines the number of runs that are constructed with respect
// to the median of the dataset
// This is similar to a normal run, but instead of being compared
// to the previous value, each value is compared to the median
//
// Requires data from alt_sequence2
unsigned int num_directional_runs(const vector<int> &alt_seq){
	unsigned int num_runs = 0;

	//Account for the first run (which always exists for non-empty strings)
	if(alt_seq.size() > 0) num_runs ++;

	// openmp optimization
	#pragma omp parallel for reduction(+:num_runs)
	for(unsigned int i = 1; i < alt_seq.size(); ++i){
		if(alt_seq[i] != alt_seq[i-1]){
			++num_runs;
		}
	}

	return num_runs;
}

// 5.1.3 Length of Directional Runs
// Determines the length of the longest run
//
// Requires data from alt_sequence1, binary data needs conversion1 first
//
//
// 5.1.6 Length of Runs Based on the Median
// Determines the length of the longest run that is constructed
// with respect to the median
//
// Requires data from alt_sequence2
unsigned int len_directional_runs(const vector<int> &alt_seq){
	unsigned int max_run = 0;
	unsigned int run = 1;

	for(unsigned int i = 1; i < alt_seq.size(); ++i){
		
		// Use if-else because if the length of the run increases, then it could still go on
		if(alt_seq[i] == alt_seq[i-1]){
			++run;
		}else{
			if(run > max_run){
				max_run = run;
			}
			run = 1;
		}
	}

	// Handle last run
	if(run > max_run){
		max_run = run;
	}

	return max_run;
}

// 5.1.4 Number of Increases and Decreases
// Determines the maximum number of increases or decreases between
// consecutive values
//
// Requires data from alt_sequence1, binary data needs conversion1 first
unsigned int num_increases_decreases(const vector<int> &alt_seq){
	unsigned int pos = 0;

	// openmp optimization
	#pragma omp parallel for reduction(+:pos)
	for(unsigned int i = 0; i < alt_seq.size(); ++i){
		if(alt_seq[i] == 1)
			++pos;
	}

	unsigned int reverse_pos = alt_seq.size() - pos;
	return max(pos, reverse_pos);
}

// Helper function to prepare for 5.1.7 and 5.1.8
// Should take about 2-3 seconds for a 1mil size list
// Maybe speed up by randomizing the check_size progression
// based on the average amount until the next collision (about 16-19 for a full 256 alphabet size)
// Just advance maybe 4-5 per iteration and back up if collisions are found based on how many collisions are found.
// POC: https://repl.it/C4eI
vector<unsigned int> find_collisions(const byte data[], const unsigned int n){
	vector<unsigned int> ret;
	set<unsigned int> dups;

	unsigned int i = 0;
	unsigned int check_size;

	// Begin from each element
	while(i < n){
		check_size = 0;

		// Progressively increase the number of elements checked
		while(check_size < (n - i)){

			// Toss elements into a set
			dups.insert(data[i+check_size]);

			// If sizes don't match up then a collision exists
			if(dups.size() != check_size+1){

				// Record info on collision and end inner loop
				// Advance outer loop past the collision end
				ret.push_back(check_size+1);
				i += check_size;
				check_size = n;
				dups.clear();
			}

			++check_size;
		}

		++i;
	}

	return ret;
}

// 5.1.7 Average Collision Test
// Counts the number of successive samples until a duplicate is found
//
// Requires non-binary data or binary data from conversion2
double avg_collision(const vector<unsigned int> &col_seq){

	return divide(sum(col_seq), col_seq.size());
}

// 5.1.8 Maximum Collision Test
// Determines the maximum number of samples without a duplicate
//
// Requires non-binary data or binary data from conversion2
unsigned int max_collision(const vector<unsigned int> &col_seq){
	unsigned int max = 0;
	for(unsigned int i = 0; i < col_seq.size(); ++i){
		if(max < col_seq[i]) max = col_seq[i];
	}

	return max;
}

// 5.1.9 Periodicity Test
// Determines the number of periodic structures
// Based on lag parameter p = [1, 2, 8, 16, 32]
//
// Requires non-binary data or binary data from conversion1
unsigned int periodicity(const byte data[], const unsigned int p, const unsigned int n){
	unsigned int T = 0;

	#pragma omp parallel for reduction(+:T)
	for(unsigned int i = 0; i < n-p; ++i){
		if(data[i] == data[i+p]){
			++T;
		}
	}

	return T;
}

// 5.1.10 Covariance Test
// Measures the strength of lagged correlation
// Based on lag parameter p = [1, 2, 8, 16, 32]
//
// Requires non-binary data or binary data from conversion1
unsigned long int covariance(const byte data[], const unsigned int p, const unsigned int n){
	unsigned long int T = 0;

	#pragma omp parallel for reduction(+:T)
	for(unsigned int i = 0; i < n-p; ++i){
		T += data[i] * data[i+p];
	}

	return T;
}

// 5.1.11 Compression Test
// Compresses the data using bzip2 and determines the length
// of the resulting compressed data
//
// Can handle binary and non-binary data
unsigned int compression(const byte data[], const int sample_size){
	
	// Build string of bytes
	string msg;
	char buffer[8];

	//Reserve the necessary size sample_size*(floor(log10(2^8))+1)
	msg.reserve(4*sample_size);

	for(unsigned int i = 0; i < sample_size; ++i){
		sprintf(buffer, "%u ", data[i]);
		msg += buffer;
	}

	// Remove the extra ' ' at the end
	msg.pop_back();

	// Set up structures for compression
	char* source = (char*)msg.c_str();
	unsigned int dest_len = ceil(1.01*strlen(source)) + 600;
	char* dest = new char[dest_len];

	// Compress and capture the size of the compressed data
	int rc = BZ2_bzBuffToBuffCompress(dest, &dest_len, source, strlen(source), 5, 0, 0);

	// Free memory
	delete[](dest);

	// Return with proper return code
	if(rc == BZ_OK){
		return dest_len;
	}else{
		return 0;
	}
}

/*
* ---------------------------------------------
* 	  HELPERS FOR PERMUTATION TEST ITERATION
* ---------------------------------------------
*/

void excursion_test(const byte data[], const double mean, const int sample_size, map<string, long double> &stats){

	stats["excursion"] = excursion(data, mean, sample_size);
}

void directional_tests(const byte data[], const int alphabet_size, const int sample_size, map<string, long double> &stats){

	vector<int> alt_seq;

	if(alphabet_size == 2){
		vector<byte> cs1 = conversion1(data, sample_size);
		alt_seq = alt_sequence1(cs1.data(), sample_size/8);		// conversion1 reduces the total size by a factor of 8
	}else{
		alt_seq = alt_sequence1(data, sample_size);
	}

	stats["numDirectionalRuns"] = num_directional_runs(alt_seq);
	stats["lenDirectionalRuns"] = len_directional_runs(alt_seq);
	stats["numIncreasesDecreases"] = num_increases_decreases(alt_seq);
}

void consecutive_runs_tests(const byte data[], const double median, const int alphabet_size, const int sample_size, map<string, long double> &stats){

	vector<int> alt_seq;

	if(alphabet_size == 2){
		vector<byte> cs2 = conversion2(data, sample_size);
		alt_seq = alt_sequence2(cs2.data(), 0.5, sample_size/8);	// conversion2 reduces the total size by a factor of 8
	}else{
		alt_seq = alt_sequence2(data, median, sample_size);
	}

	stats["numRunsMedian"] = num_directional_runs(alt_seq);
	stats["lenRunsMedian"] = len_directional_runs(alt_seq);
}

void collision_tests(const byte data[], const int alphabet_size, const int sample_size, map<string, long double> &stats){

	vector<unsigned int> col_seq;

	if(alphabet_size == 2){
		vector<byte> cs2 = conversion2(data, sample_size);
		col_seq = find_collisions(cs2.data(), sample_size/8);		// conversion2 reduces the total size by a factor of 8
	}else{
		col_seq = find_collisions(data, sample_size);
	}

	stats["avgCollision"] = avg_collision(col_seq);
	stats["maxCollision"] = max_collision(col_seq);
}

void periodicity_tests(const byte data[], const int alphabet_size, const int sample_size, map<string, long double> &stats){

	if(alphabet_size == 2){
		vector<byte> cs1 = conversion1(data, sample_size);
		stats["periodicity(1)"] = periodicity(cs1.data(), 1, sample_size/8);
		stats["periodicity(2)"] = periodicity(cs1.data(), 2, sample_size/8);
		stats["periodicity(8)"] = periodicity(cs1.data(), 8, sample_size/8);
		stats["periodicity(16)"] = periodicity(cs1.data(), 16, sample_size/8);
		stats["periodicity(32)"] = periodicity(cs1.data(), 32, sample_size/8);
	}else{
		stats["periodicity(1)"] = periodicity(data, 1, sample_size);
		stats["periodicity(2)"] = periodicity(data, 2, sample_size);
		stats["periodicity(8)"] = periodicity(data, 8, sample_size);
		stats["periodicity(16)"] = periodicity(data, 16, sample_size);
		stats["periodicity(32)"] = periodicity(data, 32, sample_size);
	}
}

void covariance_tests(const byte data[], const int alphabet_size, const int sample_size, map<string, long double> &stats){

	if(alphabet_size == 2){
		vector<byte> cs1 = conversion1(data, sample_size);
		stats["covariance(1)"] = covariance(cs1.data(), 1, sample_size/8);		// top should be cs1
		stats["covariance(2)"] = covariance(cs1.data(), 2, sample_size/8);
		stats["covariance(8)"] = covariance(cs1.data(), 8, sample_size/8);
		stats["covariance(16)"] = covariance(cs1.data(), 16, sample_size/8);
		stats["covariance(32)"] = covariance(cs1.data(), 32, sample_size/8);
	}else{
		stats["covariance(1)"] = covariance(data, 1, sample_size);
		stats["covariance(2)"] = covariance(data, 2, sample_size);
		stats["covariance(8)"] = covariance(data, 8, sample_size);
		stats["covariance(16)"] = covariance(data, 16, sample_size);
		stats["covariance(32)"] = covariance(data, 32, sample_size);
	}
}

void compression_test(const byte data[], const int sample_size, map<string, long double> &stats){

	stats["compression"] = compression(data, sample_size);
}

void run_tests(const byte data[], const double mean, const double median, const int alphabet_size, const int sample_size, map<string, long double> &stats){

	// Perform tests
	//double start_time = omp_get_wtime();

	excursion_test(data, mean, sample_size, stats);
	directional_tests(data, alphabet_size, sample_size, stats);
	consecutive_runs_tests(data, median, alphabet_size, sample_size, stats);
	collision_tests(data, alphabet_size, sample_size, stats);
	periodicity_tests(data, alphabet_size, sample_size, stats);
	covariance_tests(data, alphabet_size, sample_size, stats);
	compression_test(data, sample_size, stats);

	//cout << endl << "Iteration time elapsed: " << omp_get_wtime() - start_time << endl;
}

/*
* ---------------------------------------------
* 			  PERMUTATION TEST
* ---------------------------------------------
*/

void print_results(map<int, array<int, 2>> &C){
	cout << endl << endl;
	cout << "                statistic  C[i][0]  C[i][1]" << endl;
	cout << "-------------------------------------------" << endl;
	for(unsigned int i = 0; i < num_tests; ++i){
		if((C[i][0] + C[i][1] <= 5) || C[i][0] >= PERMS-5){
			cout << setw(24) << test_names[i] << "*";
		}else{
			cout << setw(25) << test_names[i];
		}
		cout << setw(8) << C[i][0];
		cout << setw(8) << C[i][1] << endl;
	}
	cout << "(* denotes failed test)" << endl;
	cout << endl;
}

bool permutation_tests(const byte ds[], const double mean, const double median, const int alphabet_size, const int sample_size, const int num_threads, const bool verbose){

	uint64_t xoshiro256starstarSeed[4];

	// We need a copy because the tests take in by reference and modify it
	byte data[sample_size];
	for(int i = 0; i < sample_size; ++i){
		data[i] = ds[i];
	}

	// Counters for the pass/fail of each statistic
	map<int, array<int, 2>> C;

	// Original test results (t) and permuted test results (t' or tp)
	map<string, long double> t, tp;

	// Build map of results
	for(unsigned int i = 0; i < num_tests; ++i){
		C[i] = {0};

		t[test_names[i]] = -1;
		tp[test_names[i]] = -1;
	}

	// Run initial tests
	cout << "Beginning initial tests..." << endl;
	seed(xoshiro256starstarSeed);
	run_tests(data, mean, median, alphabet_size, sample_size, t);

	if(verbose){
		cout << endl << "Initial test results" << endl;
		for(int i = 0; i < num_tests; i++){
			cout << setw(23) << test_names[i] << ": ";
			cout << t[test_names[i]] << endl;
		}
		cout << endl;
	}
	
	cout << "Beginning permutation tests... these may take some time" << endl;

	for(int i = 0; i < PERMS; ++i) {
		if(verbose){
			cout << "\rPermutation Test: " << divide(i, PERMS)*100 << "% complete" << flush;
		}

		shuffle(data, sample_size, xoshiro256starstarSeed);
		run_tests(data, mean, median, alphabet_size, sample_size, tp);

		// Aggregate results into the counters
		for(unsigned int j = 0; j < num_tests; ++j){
			if(tp[test_names[j]] > t[test_names[j]]){
				C[j][0]++;
			}else if(tp[test_names[j]] == t[test_names[j]]){
				C[j][1]++;
			}
		}
	}

	if(verbose) print_results(C);

	for(unsigned int i = 0; i < num_tests; ++i){
		if((C[i][0] + C[i][1] <= 5) || C[i][0] >= PERMS-5){
			return false;
	 	}
	}

	return true;
}
