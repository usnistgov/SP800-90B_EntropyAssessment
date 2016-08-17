#pragma once

#include "bzlib.h" // sudo apt-get install libbz2-dev
#include "ThreadPool.h"

#include "../shared/utils.h"

// Binary partition size
#define BIN_PART_SIZE (SIZE/8)

// The tests used
const unsigned int num_tests = 19;
const string test_names[] = {"excursion","numDirectionalRuns","lenDirectionalRuns","numIncreasesDecreases","numRunsMedian","lenRunsMedian","avgCollision","maxCollision","periodicity(1)","periodicity(2)","periodicity(8)","periodicity(16)","periodicity(32)","covariance(1)","covariance(2)","covariance(8)","covariance(16)","covariance(32)","compression"};

/*
* ---------------------------------------------
* 	  TASKS FOR PERMUTATION TESTS
* ---------------------------------------------
*/

// 5.1 Conversion I
// Takes a binary sequence and partitions it into 8-bit blocks
// Blocks have the number of 1's counted and totaled
array<byte, BIN_PART_SIZE> conversion1(const byte data[]){
	array<byte, BIN_PART_SIZE> ret = {0};

	for(unsigned int i = 0; i < SIZE; i+=8){
		for(unsigned int j = 0; j < 8; j++){
			ret[i/8] += data[i+j];
		}
	}

	return ret;
}

// 5.1 Conversion II
// Takes a binary sequence and partitions it into 8-bit blocks
// Blocks are then converted to decimal
array<byte, BIN_PART_SIZE> conversion2(const byte data[]){
	array<byte, BIN_PART_SIZE> ret = {0};

	for(unsigned int i = 0; i < SIZE; i+=8){
		for(unsigned int j = 0; j < 8; j++){
			ret[i/8] += data[i+j] * pow(2, (7-j));
		}
	}

	return ret;
}

// 5.1.1 Excursion Test
// Measures how far the running sum of values deviates from the
// average value at each point in the set
double excursion(const byte data[], const double mean){
	double d_i = 0;
	double max = 0;
	double running_sum = 0;

	for(unsigned int i = 0; i < SIZE; i++){
		running_sum += data[i];
		d_i = running_sum - ((i+1) * mean);

		if(d_i > max){
			max = d_i;
		}else if(-1 * d_i > max){
			max = -1 * d_i;
		}
	}

	return max;
}

// Helper for 5.1.2, 5.1.3, and 5.1.4
// Builds a array of the runs of consecutive values
// Pushes +1 to the array if the value is >= the previous
// Pushes -1 to the array if the value is < than the previous
vector<int> alt_sequence1(const byte data[], const unsigned int n){
	vector<int> ret(n-1, 0);

	for(unsigned int i = 0; i < n-1; i++){
		ret[i] = ((data[i] > data[i+1]) ? -1 : 1);
	}

	return ret;
}

// Helper for 5.1.5 and 5.1.6
// Builds a array of the runs of values compared to the median
// Pushes +1 to the array if the value is >= the median
// Pushes -1 to the array if the value is < than the median
vector<int> alt_sequence2(const byte data[], const double median, const unsigned int n){
	vector<int> ret(n, 0);

	for(unsigned int i = 0; i < n; i++){
		ret[i] = ((data[i] < median) ? -1 : 1);
	}

	return ret;
}

// 5.1.2 Number of Directional Runs
// Determines the number of runs in the sequence.
// A run is when multiple consecutive values are all >= the prior
// or all < the prior
// 5.1.5 Number of Runs Based on the Median
// Determines the number of runs that are constructed with respect
// to the median of the dataset
// This is similar to a normal run, but instead of being compared
// to the previous value, each value is compared to the median
unsigned int num_directional_runs(const vector<int> &alt_seq){
	unsigned int num_runs = 0;
	for(unsigned int i = 1; i < alt_seq.size(); i++){
		if(alt_seq[i] != alt_seq[i-1]){
			num_runs++;
		}
	}

	return num_runs;
}

// 5.1.3 Length of Directional Runs
// Determines the length of the longest run
// 5.1.6 Length of Runs Based on the Median
// Determines the length of the longest run that is constructed
// with respect to the median
unsigned int len_directional_runs(const vector<int> &alt_seq){
	unsigned int max_run = 0;
	unsigned int run = 1;

	for(unsigned int i = 1; i < alt_seq.size(); i++){
		if(alt_seq[i] == alt_seq[i-1]){
			run++;
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
unsigned int num_increases_decreases(const vector<int> &alt_seq){
	unsigned int pos = 0;
	for(unsigned int i = 0; i < alt_seq.size(); i++){
		if(alt_seq[i] == 1) pos++;
	}

	unsigned reverse_pos = alt_seq.size() - pos;
	return max(pos, reverse_pos);
}

// Helper function to prepare for 5.1.7 and 5.1.8
// Should take about 2-3 seconds for a 1mil size list
// Maybe speed up by randomizing the check_size progression
// based on the average amount until the next collision (about 16-19)
// Just advance maybe 4-5 per iteration and back up if collisions are found based on how many collision are found.
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

			check_size++;
		}

		i++;
	}

	return ret;
}

// 5.1.7 Average Collision Test
// Counts the number of successive samples until a duplicate is found
double avg_collision(const vector<unsigned int> &col_seq){

	return divide(sum(col_seq), col_seq.size());
}

// 5.1.8 Maximum Collision Test
// Determines the maximum number of samples without a duplicate
unsigned int max_collision(const vector<unsigned int> &col_seq){
	unsigned int max = 0;
	for(unsigned int i = 0; i < col_seq.size(); i++){
		if(max < col_seq[i]) max = col_seq[i];
	}

	return max;
}

// 5.1.9 Periodicity Test
// Determines the number of periodic structures
// Based on lag parameter p
unsigned int periodicity(const byte data[], const unsigned int p, const unsigned int n){
	unsigned int T = 0;

	for(unsigned int i = 0; i < n-p; i++){
		if(data[i] == data[i+p]){
			T++;
		}
	}

	return T;
}

// 5.1.10 Covariance Test
// Measures the strength of lagged correlation
// Based on lag parameter p
unsigned long int covariance(const byte data[], const unsigned int p, const unsigned int n){
	unsigned long int T = 0;

	for(unsigned int i = 0; i < n-p; i++){
		T += data[i] * data[i+p];
	}

	return T;
}

// 5.1.11 Compression Test
// Compresses the data using bzip2 and determines the length
// of the resulting compressed data
unsigned int compression(const byte data[]){

	// Build string of bytes
	string msg = "";
	for(unsigned int i = 0; i < SIZE; i++){
		msg += to_string((unsigned int)data[i]);		// Dependent on C++11
		msg += " ";
	}

	// Set up structures for compression
	char* source = (char*)msg.c_str();
	unsigned int dest_len = 2*SIZE;
	char* dest = new char[dest_len];

	// Compress
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

void excursion_test(const byte data[], const double mean, map<string, long double> &stats){

	stats["excursion"] = excursion(data, mean);
}

void directional_tests(const byte data[], const bool is_binary, map<string, long double> &stats){

	vector<int> alt_seq;

	if(is_binary){
		array<byte, BIN_PART_SIZE> cs1 = conversion1(data);
		alt_seq = alt_sequence1(cs1.data(), BIN_PART_SIZE);
	}else{
		alt_seq = alt_sequence1(data, SIZE);
	}

	stats["numDirectionalRuns"] = num_directional_runs(alt_seq);
	stats["lenDirectionalRuns"] = len_directional_runs(alt_seq);
	stats["numIncreasesDecreases"] = num_increases_decreases(alt_seq);
}

void consecutive_runs_tests(const byte data[], const double median, const bool is_binary, map<string, long double> &stats){

	vector<int> alt_seq;

	if(is_binary){
		array<byte, BIN_PART_SIZE> cs2 = conversion2(data);
		alt_seq = alt_sequence2(cs2.data(), 0.5, BIN_PART_SIZE);
	}else{
		alt_seq = alt_sequence2(data, median, SIZE);
	}

	stats["numRunsMedian"] = num_directional_runs(alt_seq);
	stats["lenRunsMedian"] = len_directional_runs(alt_seq);
}

void collision_tests(const byte data[], const bool is_binary, map<string, long double> &stats){

	vector<unsigned int> col_seq;

	if(is_binary){
		array<byte, BIN_PART_SIZE> cs2 = conversion2(data);
		col_seq = find_collisions(cs2.data(), BIN_PART_SIZE);
	}else{
		col_seq = find_collisions(data, SIZE);
	}

	stats["avgCollision"] = avg_collision(col_seq);
	stats["maxCollision"] = max_collision(col_seq);
}

void periodicity_tests(const byte data[], const bool is_binary, map<string, long double> &stats){

	if(is_binary){
		array<byte, BIN_PART_SIZE> cs1 = conversion1(data);
		stats["periodicity(1)"] = periodicity(cs1.data(), 1, BIN_PART_SIZE);
		stats["periodicity(2)"] = periodicity(cs1.data(), 2, BIN_PART_SIZE);
		stats["periodicity(8)"] = periodicity(cs1.data(), 8, BIN_PART_SIZE);
		stats["periodicity(16)"] = periodicity(cs1.data(), 16, BIN_PART_SIZE);
		stats["periodicity(32)"] = periodicity(cs1.data(), 32, BIN_PART_SIZE);
	}else{
		stats["periodicity(1)"] = periodicity(data, 1, SIZE);
		stats["periodicity(2)"] = periodicity(data, 2, SIZE);
		stats["periodicity(8)"] = periodicity(data, 8, SIZE);
		stats["periodicity(16)"] = periodicity(data, 16, SIZE);
		stats["periodicity(32)"] = periodicity(data, 32, SIZE);
	}
}

void covariance_tests(const byte data[], const bool is_binary, map<string, long double> &stats){

	if(is_binary){
		array<byte, BIN_PART_SIZE> cs1 = conversion1(data);
		stats["covariance(1)"] = covariance(cs1.data(), 1, BIN_PART_SIZE);		// top should be cs1
		stats["covariance(2)"] = covariance(cs1.data(), 2, BIN_PART_SIZE);
		stats["covariance(8)"] = covariance(cs1.data(), 8, BIN_PART_SIZE);
		stats["covariance(16)"] = covariance(cs1.data(), 16, BIN_PART_SIZE);
		stats["covariance(32)"] = covariance(cs1.data(), 32, BIN_PART_SIZE);
	}else{
		stats["covariance(1)"] = covariance(data, 1, SIZE);
		stats["covariance(2)"] = covariance(data, 2, SIZE);
		stats["covariance(8)"] = covariance(data, 8, SIZE);
		stats["covariance(16)"] = covariance(data, 16, SIZE);
		stats["covariance(32)"] = covariance(data, 32, SIZE);
	}
}

void compression_test(const byte data[], map<string, long double> &stats){

	stats["compression"] = compression(data);
}

void run_tests(const byte data[], const double mean, const double median, const bool is_binary, map<string, long double> &stats){

	// Perform tests
	excursion_test(data, mean, stats);
	directional_tests(data, is_binary, stats);
	consecutive_runs_tests(data, median, is_binary, stats);	
	collision_tests(data, is_binary, stats);
	periodicity_tests(data, is_binary, stats);
	covariance_tests(data, is_binary, stats);
	compression_test(data, stats);
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
	for(int i = 0; i < num_tests; i++){
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

bool permutation_tests(const byte ds[], const double mean, const double median, const bool is_binary, const int num_threads, const bool verbose){

	// We need a copy because the tests take in by reference and modify it
	byte data[SIZE];
	for(int i = 0; i < SIZE; i++){
		data[i] = ds[i];
	}

	// Counters for the pass/fail of each statistic
	map<int, array<int, 2>> C;

	// Original test results (t) and permuted test results (t' or tp)
	map<string, long double> t;
	// map<string, long double> tp;

	// Build map of results
	for(int i = 0; i < num_tests; i++){
		C[i] = {0};

		t[test_names[i]] = -1;
		// tp[test_names[i]] = -1;
	}

	// Run initial tests
	cout << "Beginning initial tests..." << endl;
	run_tests(data, mean, median, is_binary, t);

	/*
	* if(verbose){
	* 	cout << endl << "Initial test results" << endl;
	* 	for(int i = 0; i < num_tests; i++){
	* 		cout << setw(23) << test_names[i] << ": ";
	* 		cout << t[test_names[i]] << endl;
	* 	}
	* 	cout << endl;
	* }
	*/

	// Permutation tests, shuffle -> run -> aggregate
	ThreadPool pool(num_threads);
	vector<future<map<string, long double>>> results;

	cout << "Beginning permutation tests..." << endl;
	for(int i = 0; i < PERMS; i++){

		if(verbose){
			cout << "\rPermutation Test: " << divide(i, PERMS)*100 << "% complete" << flush;
		}

		auto lambda_test = [&](){
			map<string, long double> tp;
			for(int i = 0; i < num_tests; i++){
				tp[test_names[i]] = -1;
			}
			shuffle(data);
			run_tests(data, mean, median, is_binary, tp);
			return tp;
		};

		results.emplace_back(pool.enqueue(lambda_test));

		// shuffle(data);
		// run_tests(data, mean, median, is_binary, tp);
	}

	for(auto &&result : results){
		auto tmp_result = result.get();
		for(int i = 0; i < num_tests; i++){
			if(tmp_result[test_names[i]] > t[test_names[i]]){
				C[i][0]++;
			}else if(tmp_result[test_names[i]] == t[test_names[i]]){
				C[i][1]++;
			}
		}
	}

	// // Aggregate results into the counters
	// for(int j = 0; j < num_tests; j++){
	// 	if(tp[test_names[j]] > t[test_names[j]]){
	// 		C[j][0]++;
	// 	}else if(tp[test_names[j]] == t[test_names[j]]){
	// 		C[j][1]++;
	// 	}
	// }

	if(verbose) print_results(C);

	for(int i = 0; i < num_tests; i++){
		if((C[i][0] + C[i][1] <= 5) || C[i][0] >= PERMS-5){
			return false;
	 	}
	}

	return true;
}
