#pragma once

#include <stdlib.h>
#include <bzlib.h> // sudo apt-get install libbz2-dev
#include "../shared/utils.h"
#include "../shared/TestCase.h"
#include <assert.h>
#include <unistd.h>

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
vector<uint8_t> conversion1(const uint8_t data[], const int sample_size){
	vector<uint8_t> ret((sample_size / 8) + ((sample_size%8==0)?0:1), 0);

	for(int i = 0; i < sample_size; ++i){
		ret[i/8] += data[i];	// integer division to ensure the size of ret is sample_size / 8
	}

	return ret;
}

// 5.1 Conversion II
// Takes a binary sequence and partitions it into 8-bit blocks
// Blocks are then converted to decimal
//
// Requires binary data
vector<uint8_t> conversion2(const uint8_t data[], const int sample_size){
	vector<uint8_t> ret((sample_size / 8) + ((sample_size%8==0)?0:1), 0);

	for(int i = 0; i < sample_size; ++i) {
		ret[i/8] += data[i] << (7 - i%8);
	}

	return ret;
}

// 5.1.1 Excursion Test
// Measures how far the running sum of values deviates from the
// average value at each point in the set
//
// Requires binary or non-binary data
double excursion(const uint8_t data[], const double rawmean, const int sample_size){
	double d_i = 0;
	double max = 0;
	double running_sum = 0;

	for(int i = 0; i < sample_size; ++i){
		running_sum += data[i];
		d_i = abs(running_sum - ((i+1) * rawmean));

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
vector<int> alt_sequence1(const uint8_t data[], const int sample_size){
	vector<int> ret(sample_size-1, 0);

	for(int i = 0; i < sample_size-1; ++i){
		ret[i] = ((data[i] > data[i+1]) ? -1 : 1);
	}

	return ret;
}

// Helper for 5.1.5 and 5.1.6
// Builds a array of the runs of values compared to the median
// Pushes +1 to the array if the value is >= the median
// Pushes -1 to the array if the value is < than the median
vector<int> alt_sequence2(const uint8_t data[], const double median, const int sample_size){
	vector<int> ret(sample_size, 0);

	for(int i = 0; i < sample_size; ++i){
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
	for(unsigned int i = 0; i < alt_seq.size(); ++i){
		if(alt_seq[i] == 1)
			++pos;
	}

	unsigned int reverse_pos = alt_seq.size() - pos;
	return max(pos, reverse_pos);
}

// Helper function to prepare for 5.1.7 and 5.1.8
vector<unsigned int> find_collisions(const uint8_t data[], const unsigned int n, const unsigned int k){
	vector<unsigned int> ret;
	vector<bool> dups(k, false);

	unsigned long int i=0;
	unsigned long int j=0;

	// Begin at the start
	while(i + j < n){
		for(unsigned int l=0; l<k; l++) dups[l] = false;

		// Progressively increase the number of elements checked
		while(i + j < n) {
			// Check for a collision
			if(dups[data[i+j]]) {
				// Record info on collision and end inner loop
				// Advance outer loop past the collision end
				// Note that j is not the current window size,
				// it is one less that the window size.
				ret.push_back(j+1);
				i += j;
				j=0;
				break;
			} else {
				dups[data[i+j]]=true;
				++j;
			}

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
unsigned int periodicity(const uint8_t data[], const unsigned int p, const unsigned int n){
	unsigned int T = 0;

	assert(n>=p);

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
unsigned long int covariance(const uint8_t data[], const unsigned int p, const unsigned int n){
	unsigned long int T = 0;

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
unsigned int compression(const uint8_t data[], const int sample_size, const uint8_t max_symbol){
	char *msg;
	unsigned int curlen = 0;
	char *curmsg;

	assert(max_symbol > 0);

	// Build string of bytes
	// Reserve the necessary size sample_size*(floor(log10(max_symbol))+2)
	// This is "worst case" and accounts for the space at the end of the number, as well.
	msg = new char[(size_t)(floor(log10(max_symbol))+2.0)*sample_size+1];
	msg[0] = '\0';
	curmsg = msg;

	for(int i = 0; i < sample_size; ++i) {
		int res;
		res = sprintf(curmsg, "%u ", data[i]);
		assert(res >= 2);
		curlen += res;
		curmsg += res;
	}

	if(curlen > 0) {
		// Remove the extra ' ' at the end
		assert(curmsg > msg);
		curmsg--;
		*curmsg = '\0';
		curlen--;
	}

	// Set up structures for compression
	unsigned int dest_len = ceil(1.01*curlen) + 600;
	char* dest = new char[dest_len];

	// Compress and capture the size of the compressed data
	int rc = BZ2_bzBuffToBuffCompress(dest, &dest_len, msg, curlen, 5, 0, 0);

	// Free memory
	delete[](dest);
	delete[](msg);

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

void excursion_test(const uint8_t data[], const double rawmean, const int sample_size, long double* stats, const bool *test_status){

	if(test_status[0]) stats[0] = excursion(data, rawmean, sample_size);
}

void directional_tests(const uint8_t data[], const int alphabet_size, const int sample_size, long double *stats, const bool *test_status){

	vector<int> alt_seq;

	if(test_status[1] || test_status[2] || test_status[3]) {
		if(alphabet_size == 2){
			vector<uint8_t> cs1 = conversion1(data, sample_size);
			alt_seq = alt_sequence1(cs1.data(), cs1.size());		// conversion1 reduces the total size by a factor of 8
		}else{
			alt_seq = alt_sequence1(data, sample_size);
		}

		if(test_status[1]) stats[1] = num_directional_runs(alt_seq);
		if(test_status[2]) stats[2] = len_directional_runs(alt_seq);
		if(test_status[3]) stats[3] = num_increases_decreases(alt_seq);
	}
}

void consecutive_runs_tests(const uint8_t data[], const double median, const int alphabet_size, const int sample_size, long double *stats, const bool *test_status){

	vector<int> alt_seq;

	if(test_status[4] || test_status[5]) {
		if(alphabet_size == 2){
			alt_seq = alt_sequence2(data, 0.5, sample_size);
		}else{
			alt_seq = alt_sequence2(data, median, sample_size);
		}

		if(test_status[4]) stats[4] = num_directional_runs(alt_seq);
		if(test_status[5]) stats[5] = len_directional_runs(alt_seq);
	}
}

void collision_tests(const uint8_t data[], const int alphabet_size, const int sample_size, long double *stats, const bool *test_status){

	vector<unsigned int> col_seq;

	if(test_status[7] || test_status[6]) {
		if(alphabet_size == 2){
			vector<uint8_t> cs2 = conversion2(data, sample_size);
			col_seq = find_collisions(cs2.data(), cs2.size(), 256);		// conversion2 reduces the total size by a factor of 8
		}else{
			col_seq = find_collisions(data, sample_size, alphabet_size);
		}

		if(test_status[6]) stats[6] = avg_collision(col_seq);
		if(test_status[7]) stats[7] = max_collision(col_seq);
	}
}

void periodicity_tests(const uint8_t data[], const int alphabet_size, const int sample_size,long double *stats, const bool *test_status){

	if(test_status[8] || test_status[9] || test_status[10] || test_status[11] || test_status[12]) {
		if(alphabet_size == 2){
			vector<uint8_t> cs1 = conversion1(data, sample_size);
			if(test_status[8]) stats[8] = periodicity(cs1.data(), 1, cs1.size());
			if(test_status[9]) stats[9] = periodicity(cs1.data(), 2, cs1.size());
			if(test_status[10]) stats[10] = periodicity(cs1.data(), 8, cs1.size());
			if(test_status[11]) stats[11] = periodicity(cs1.data(), 16, cs1.size());
			if(test_status[12]) stats[12] = periodicity(cs1.data(), 32, cs1.size());
		}else{
			if(test_status[8]) stats[8] = periodicity(data, 1, sample_size);
			if(test_status[9]) stats[9] = periodicity(data, 2, sample_size);
			if(test_status[10]) stats[10] = periodicity(data, 8, sample_size);
			if(test_status[11]) stats[11] = periodicity(data, 16, sample_size);
			if(test_status[12]) stats[12] = periodicity(data, 32, sample_size);
		}
	}
}

void covariance_tests(const uint8_t data[], const int alphabet_size, const int sample_size, long double *stats, const bool *test_status){

	if(test_status[13] || test_status[14] || test_status[15] || test_status[16] || test_status[17]) {
		if(alphabet_size == 2){
			vector<uint8_t> cs1 = conversion1(data, sample_size);
			if(test_status[13]) stats[13] = covariance(cs1.data(), 1, cs1.size());		// top should be cs1
			if(test_status[14]) stats[14]  = covariance(cs1.data(), 2, cs1.size());
			if(test_status[15]) stats[15]  = covariance(cs1.data(), 8, cs1.size());
			if(test_status[16]) stats[16]  = covariance(cs1.data(), 16, cs1.size());
			if(test_status[17]) stats[17]  = covariance(cs1.data(), 32, cs1.size());
		}else{
			if(test_status[13]) stats[13]  = covariance(data, 1, sample_size);
			if(test_status[14]) stats[14] = covariance(data, 2, sample_size);
			if(test_status[15]) stats[15] = covariance(data, 8, sample_size);
			if(test_status[16]) stats[16] = covariance(data, 16, sample_size);
			if(test_status[17]) stats[17] = covariance(data, 32, sample_size);
		}
	}
}

void compression_test(const uint8_t data[], const int sample_size, long double *stats, const uint8_t max_symbol, const bool *test_status){

	if(test_status[18]) stats[18] = compression(data, sample_size, max_symbol);
}

void run_tests(const data_t *dp, const uint8_t data[], const uint8_t rawdata[], const double rawmean, const double median, long double *stats, const bool *test_status){

	// Perform tests
	excursion_test(rawdata, rawmean, dp->len, stats, test_status);
	directional_tests(data, dp->alph_size, dp->len, stats, test_status);
	consecutive_runs_tests(data, median, dp->alph_size, dp->len, stats, test_status);
	collision_tests(data, dp->alph_size, dp->len, stats, test_status);
	periodicity_tests(data, dp->alph_size, dp->len, stats, test_status);
	if(dp->alph_size == 2) {
		//The two conversions only make sense if the two symbols are 0 and 1.
		covariance_tests(data, dp->alph_size, dp->len, stats, test_status);
	} else {
		covariance_tests(rawdata, dp->alph_size, dp->len, stats, test_status);
	}
	compression_test(rawdata, dp->len, stats, dp->maxsymbol, test_status);
}

/*
 * ---------------------------------------------
 * 			  PERMUTATION TEST
 * ---------------------------------------------
 */

void print_results(int C[][3], const int verbose){
	cout << endl << endl;
	cout << "                statistic  C[i][0]  C[i][1]  C[i][2]" << endl;
	cout << "----------------------------------------------------" << endl;
	for(unsigned int i = 0; i < num_tests; ++i){
		if((C[i][0] + C[i][1] <= 5) || C[i][1] + C[i][2] <= 5){
			cout << setw(24) << test_names[i] << "*";
		}else{
			cout << setw(25) << test_names[i];
		}
		cout << setw(8) << C[i][0];
		cout << setw(8) << C[i][1];
		cout << setw(8) << C[i][2] << endl;
	}
	cout << "(* denotes failed test)" << endl;
	cout << endl;
}

void populateTestCase(IidTestCase &tc, int C[][3]){

    PermutationTestResult tr0;
    PermutationTestResult tr1;
    PermutationTestResult tr2;
    
    tr0.iteration = 0;
    tr1.iteration = 1;
    tr2.iteration = 2;
    
    tr0.excursion = C[0][0];
    tr1.excursion = C[0][1];
    tr2.excursion = C[0][2];

    tr0.numDirectionalRuns = C[1][0];
    tr1.numDirectionalRuns = C[1][1];
    tr2.numDirectionalRuns = C[1][2];
       
    tr0.lenDirectionalRuns = C[2][0];
    tr1.lenDirectionalRuns = C[2][1];
    tr2.lenDirectionalRuns = C[2][2];
    
    tr0.numIncreasesDecreases = C[3][0];
    tr1.numIncreasesDecreases = C[3][1];
    tr2.numIncreasesDecreases = C[3][2];
    
    tr0.numRunsMedian = C[4][0];
    tr1.numRunsMedian = C[4][1];
    tr2.numRunsMedian = C[4][2];
    
    tr0.lenRunsMedian = C[5][0];
    tr1.lenRunsMedian = C[5][1];
    tr2.lenRunsMedian = C[5][2];
    
    tr0.avgCollision = C[6][0];
    tr1.avgCollision = C[6][1];
    tr2.avgCollision = C[6][2];
    
    tr0.maxCollision = C[7][0];
    tr1.maxCollision = C[7][1];
    tr2.maxCollision = C[7][2];
    
    tr0.periodicity1 = C[8][0];
    tr1.periodicity1 = C[8][1];
    tr2.periodicity1 = C[8][2];
    
    tr0.periodicity2 = C[9][0];
    tr1.periodicity2 = C[9][1];
    tr2.periodicity2 = C[9][2];
    
    tr0.periodicity8 = C[10][0];
    tr1.periodicity8 = C[10][1];
    tr2.periodicity8 = C[10][2];

    tr0.periodicity16 = C[11][0];
    tr1.periodicity16 = C[11][1];
    tr2.periodicity16 = C[11][2];

    tr0.periodicity32 = C[12][0];
    tr1.periodicity32 = C[12][1];
    tr2.periodicity32 = C[12][2];

    tr0.covariance1 = C[13][0];
    tr1.covariance1 = C[13][1];
    tr2.covariance1 = C[13][2];

    tr0.covariance2 = C[14][0];
    tr1.covariance2 = C[14][1];
    tr2.covariance2 = C[14][2];

    tr0.covariance8 = C[15][0];
    tr1.covariance8 = C[15][1];
    tr2.covariance8 = C[15][2];
    
    tr0.covariance16 = C[16][0];
    tr1.covariance16 = C[16][1];
    tr2.covariance16 = C[16][2];

    tr0.covariance32 = C[17][0];
    tr1.covariance32 = C[17][1];
    tr2.covariance32 = C[17][2];

    tr0.compression = C[18][0];
    tr1.compression = C[18][1];
    tr2.compression = C[18][2];
    
    tc.testResults.push_back(tr0);
    tc.testResults.push_back(tr1);
    tc.testResults.push_back(tr2);
}

bool permutation_tests(const data_t *dp, const double rawmean, const double median, const int verbose, IidTestCase &tc){
	uint64_t xoshiro256starstarMainSeed[4];
	bool istty;

	// Progress
	size_t completed = 0;

	// Counters for the pass/fail of each statistic
	int C[num_tests][3];

	// Original test results (t) 
	long double t[num_tests];
	bool test_status[num_tests];

	istty = (isatty(STDOUT_FILENO)==1);

	// Build map of results
	for(unsigned int i = 0; i < num_tests; ++i){
		C[i][0] = 0;
		C[i][1] = 0;
		C[i][2] = 0;

		t[i] = -1;
		test_status[i] = true;
	}

	// Run initial tests
	if(verbose == 2) cout << "Beginning initial tests..." << endl;
	seed(xoshiro256starstarMainSeed);

	run_tests(dp, dp->symbols, dp->rawsymbols, rawmean, median, t, test_status);

	if(verbose == 2) {
		cout << endl << "Initial test results" << endl;
		for(unsigned int i = 0; i < num_tests; i++){
			cout << setw(23) << test_names[i] << ": ";
			cout << t[i] << endl;
		}
		cout << endl;
	} else if (verbose >= 3) {
		for(unsigned int i = 0; i < num_tests; i++){
			printf("Permutation testing: Unpermuted result %s = %.22Lg\n", test_names[i].c_str(), t[i]);
		}
	}
	
	if(verbose == 2) cout << "Beginning permutation tests... these may take some time" << endl;

	#pragma omp parallel
	{
		uint8_t *data;
		uint8_t *rawdata;
		uint64_t xoshiro256starstarSeed[4];
		long double tp[num_tests];
		int passed_count;

		data = new uint8_t[dp->len];
		rawdata = new uint8_t[dp->len];

		// Init results
		for(unsigned int i = 0; i < num_tests; ++i){
			tp[i] = -1;
		}

		for(int i = 0; i < dp->len; ++i){
			data[i] = dp->symbols[i];
			rawdata[i] = dp->rawsymbols[i];
		}

		passed_count = 0;
		memcpy(xoshiro256starstarSeed, xoshiro256starstarMainSeed, sizeof(xoshiro256starstarMainSeed));
		//Cause the RNG to jump omp_get_thread_num() * 2^128 calls
		xoshiro_jump(omp_get_thread_num(), xoshiro256starstarSeed);

		#pragma omp for
		for(int i = 0; i < PERMS; ++i) {
			if(passed_count < 19) {
				char statusMessage[1024];
				size_t statusMessageLength = 0;

				FYshuffle(data, rawdata, dp->len, xoshiro256starstarSeed);
				run_tests(dp, data, rawdata, rawmean, median, tp, test_status);

				// Aggregate results into the counters
				#pragma omp critical(resultUpdate)
				{
					for(unsigned int j = 0; j < num_tests; ++j){
						if(test_status[j]) {
							if(tp[j] > t[j]){
								C[j][0]++;
							} else if(tp[j] == t[j]){
								C[j][1]++;
							} else {
								C[j][2]++;
							}
							if((C[j][0] + C[j][1] > 5) && (C[j][1] + C[j][2] > 5)) {
								test_status[j] = false;
							}
						}
					}
					passed_count = 0;
					for(unsigned int j=0; j < num_tests; j++) if(!test_status[j]) passed_count++;
					completed ++;
				} // end resultUpdate

				if(verbose == 2) {
					int res;
					/* Construct pretty output regardless of whether on terminal (tty) or 
					* redirected to another file descriptor (eg. redirect to file).
					* Note that if using something like 'tee' to replicate the output
					* then it might be handy to use 'unbuffer' to fake the call into
					* thinking it is still being sent to a tty.
					*/
					if(istty) {
						statusMessage[0] = '\r';
						statusMessage[1] = '\0';
						statusMessageLength = 1;
					} else {
						statusMessage[0] = '\0';
						statusMessageLength = 0;
					}

					res = snprintf(statusMessage+statusMessageLength, sizeof(statusMessage)-statusMessageLength, "%6.02f%% of Permutation test rounds, %6.02f%% of Permutation tests", (100.0*((float)completed)/((float)PERMS)), (100.0*((float)passed_count)/19.0));
					assert(res>0);
					statusMessageLength += res;
					assert(statusMessageLength < sizeof(statusMessage));

					/* If not displaying to screen, then we can print even more information. Ultimately
					* we want the '\n' however printed when not printing to terminal so that the redirected
					* output looks nicer. 
					*/
					if(!istty)  {
						res = snprintf(statusMessage+statusMessageLength, sizeof(statusMessage)-statusMessageLength, " (Core %d/%d, passed_count %d)\n", omp_get_thread_num(), omp_get_num_threads()-1, passed_count);
						assert(res>0);
						statusMessageLength += res;
						assert(statusMessageLength < sizeof(statusMessage));
					}
					#pragma omp critical(verboseOutput)
					{
						fputs(statusMessage, stdout);
						fflush(stdout);
					}
				}
			} else {
				//We don't have a lock for this branch, so make one to update the completed count.
				#pragma omp atomic
				completed ++;
			}

		}
        	delete[](data);
        	delete[](rawdata);
	} //end parallel

	if(verbose > 1) print_results(C, verbose);
        
    populateTestCase(tc, C);
	
    for(unsigned int i = 0; i < num_tests; ++i){
		if((C[i][0] + C[i][1] <= 5) || (C[i][1] + C[i][2] <= 5)){
			return false;
	 	}
	}

	return true;
}
