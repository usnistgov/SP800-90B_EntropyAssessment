#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <time.h>
#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <string.h>
#include <iomanip>

#include "bzlib.h" // sudo apt-get install libbz2-dev

#define SIZE 1000000
#define PERMS 100

using namespace std;

typedef unsigned char byte;

// Global variables
byte data[SIZE];
double mean = 0.0;
double median = 0.0;
bool is_binary = false;
map<int, int*> C;
map<string, long int> t;
map<string, long int> tp;
const int num_tests = 19;
const string test_names[] = {"excursion","numDirectionalRuns","lenDirectionalRuns","numIncreasesDecreases","numRunsMedian","lenRunsMedian","avgCollision","maxCollision","periodicity(1)","periodicity(2)","periodicity(8)","periodicity(16)","periodicity(32)","covariance(1)","covariance(2)","covariance(8)","covariance(16)","covariance(32)","compression"};

// Read in binary file to test
void read_file(const char* file_path){
	FILE* file = NULL;

	#ifdef VERBOSE
	printf("Opening file %s\n", file_path);
	#endif

	file = fopen(file_path, "rb");
	fread(data, SIZE, 1, file);
	fclose(file);

	#ifdef VERBOSE
	printf("Data read\n");
	#endif
}

// Fisher-Yates Fast (in place) shuffle algorithm
void shuffle(byte arr[]){
	srand(time(NULL));
	long int r;

	for(long int i = SIZE-1; i > 0; i--){
		r = rand() % (i+1);
		swap(arr[r], arr[i]);
	}
}

// Quick sum array
long int sum(byte arr[]){
	long int sum = 0;
	for(long int i = 0; i < SIZE; i++){
		sum += arr[i];
	}

	return sum;
}

// Quick sum vector
long int sum(vector<int> v){
	long int sum = 0;
	for(long int i = 0; i < v.size(); i++){
		sum += v[i];
	}

	return sum;
}

// Calculate baseline statistics
void calc_stats(){

	// Calculate mean
	mean = sum(data) / (long double) SIZE;

	// Sort in a vector for median/min/max
	vector<byte> v(data, data+SIZE);
	sort(v.begin(), v.end());

	byte min = v[0];
	byte max = v[SIZE-1];

	if(min == 0 && max == 1){
		is_binary = true;
		median = 0.5;
	}else{
		long int half = SIZE / 2;
		median = (v[half] + v[half-1]) / 2.0;
	}
}

double excursion(){
	double d_i = 0;
	double max = 0;
	double running_sum = 0;

	for(long int i = 0; i < SIZE; i++){
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

vector<int> alt_sequence1(){
	vector<int> ret(SIZE-1, 0);

	for(int i = 0; i < SIZE-1; i++){
		ret[i] = ((data[i] > data[i+1]) ? -1 : 1);
	}

	return ret;
}

vector<int> alt_sequence2(){
	vector<int> ret(SIZE, 0);

	for(int i = 0; i < SIZE; i++){
		ret[i] = ((data[i] < median) ? -1 : 1);
	}

	return ret;
}

long int num_directional_runs(vector<int> alt_seq){
	long int num_runs = 0;
	for(int i = 1; i < SIZE; i++){
		if(alt_seq[i] != alt_seq[i-1]){
			num_runs++;
		}
	}

	return num_runs;
}

long int len_directional_runs(vector<int> alt_seq){
	long int max_run = 0;
	long int run = 1;

	for(int i = 1; i < SIZE; i++){
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

long int num_increases_decreases(vector<int> alt_seq){
	long int pos = 0;
	for(int i = 0; i < SIZE; i++){
		if(alt_seq[i] == 1) pos++;
	}

	return max(pos, SIZE-pos);
}

long int num_runs_median(vector<int> alt_seq){
	long int num_runs = 1;

	for(int i = 1; i < SIZE; i++){
		if(alt_seq[i] != alt_seq[i-1]){
			num_runs++;
		}
	}

	return num_runs;
}

long int len_runs_median(vector<int> alt_seq){
	long int max_run = 0;
	long int run = 1;

	for(int i = 1; i < SIZE; i++){
		if(alt_seq[i] == alt_seq[i-1]){
			run++;
		}else{
			if(run > max_run){
				max_run = run;
			}
			run = 1;
		}
	}

	// Handle the last run
	if(run > max_run){
		max_run = run;
	}

	return max_run;
}

// Should take about 2-3 seconds for a 1mil size list
// Maybe speed up by randomizing the check_size progression
// based on the average amount until the next collision (about 16-19)
// Just advance maybe 4-5 per iteration and back up if collisions are found based on how many collision are found.
// POC: https://repl.it/C4eI
vector<int> find_collisions(){
	vector<int> ret;
	set<int> dups;

	long int i = 0;
	int check_size;

	// Begin from each element
	while(i < SIZE){

		check_size = 0;

		// Progressively increase the number of elements checked
		while(check_size < (SIZE - i)){

			// Toss elements into a set
			dups.insert(data[i+check_size]);

			// If sizes don't match up then a collision exists
			if(dups.size() != check_size+1){

				// Record info on collision and end inner loop
				// Advance outer loop past the collision end
				ret.push_back(check_size+1);
				i += check_size;
				check_size = SIZE;
				dups.clear();
			}

			check_size++;
		}

		i++;
	}

	return ret;
}

double avg_collision(vector<int> col_seq){

	return sum(col_seq)/float(col_seq.size());
}

long int max_collision(vector<int> col_seq){
	long int max = -1;
	for(int i = 0; i < col_seq.size(); i++){
		if(max < col_seq[i]) max = col_seq[i];
	}

	return max;
}

long int periodicity(int p){
	long int T = 0;

	for(int i = 0; i < SIZE-p; i++){
		if(data[i] == data[i+p]){
			T++;
		}
	}

	return T;
}

long int covariance(int p){
	long int T = 0;

	for(int i = 0; i < SIZE-p; i++){
		T += data[i] * data[i+p];
	}

	return T;
}

unsigned int compression(){

	// Build string of bytes
	string msg = "";
	for(int i = 0; i < SIZE; i++){
		msg += to_string((int)data[i]);		// Dependent on C++11
		msg += " ";
	}

	char* source = (char*)msg.c_str();
	unsigned int dest_len = 2*SIZE;
	char* dest = new char[dest_len];

	int rc = BZ2_bzBuffToBuffCompress(dest, &dest_len, source, strlen(source), 5, 0, 0);

	if(rc == BZ_OK){
		return dest_len;
	}else{
		return 0;
	}
}

int main(){

	// Read in file
	const char* file_path = "../bin/urand.bin";
	read_file(file_path);

	// Calculate baseline statistics
	calc_stats();

	#ifdef VERBOSE
	cout << "Mean: " << mean << endl;
	cout << "Median: " << median << endl;
	cout << (is_binary ? "Is Binary." : "Is Not Binary.") << endl;
	#endif

	// Build map of results
	for(int i = 0; i < num_tests; i++){
		C[i] = new int[2];
		C[i][0] = 0;
		C[i][1] = 0;

		t[test_names[i]] = -1;
		tp[test_names[i]] = -1;
	}

	// Begin initial tests
	t["excursion"] = excursion();

	// Directional tests
	vector<int> alt_seq = alt_sequence1();
	t["numDirectionalRuns"] = num_directional_runs(alt_seq);
	t["lenDirectionalRuns"] = len_directional_runs(alt_seq);
	t["numIncreasesDecreases"] = num_increases_decreases(alt_seq);

	// Runs tests
	alt_seq = alt_sequence2();
	t["numRunsMedian"] = num_runs_median(alt_seq);
	t["lenRunsMedian"] = len_runs_median(alt_seq);

	// Collision tests
	vector<int> col_seq = find_collisions();
	t["avgCollision"] = avg_collision(col_seq);
	t["maxCollision"] = max_collision(col_seq);

	// Periodicity tests
	t["periodicity(1)"] = periodicity(1);
	t["periodicity(2)"] = periodicity(2);
	t["periodicity(8)"] = periodicity(8);
	t["periodicity(16)"] = periodicity(16);
	t["periodicity(32)"] = periodicity(32);

	// Covariance tests
	t["covariance(1)"] = covariance(1);
	t["covariance(2)"] = covariance(2);
	t["covariance(8)"] = covariance(8);
	t["covariance(16)"] = covariance(16);
	t["covariance(32)"] = covariance(32);

	// Compression test
	t["compression"] = compression();

	// Permutation tests
	for(int i = 0; i < PERMS; i++){

		// Permute the data
		shuffle(data);

		// Excursion test
		tp["excursion"] = excursion();

		// Directional tests
		vector<int> alt_seq = alt_sequence1();
		tp["numDirectionalRuns"] = num_directional_runs(alt_seq);
		tp["lenDirectionalRuns"] = len_directional_runs(alt_seq);
		tp["numIncreasesDecreases"] = num_increases_decreases(alt_seq);

		// Runs tests
		alt_seq = alt_sequence2();
		tp["numRunsMedian"] = num_runs_median(alt_seq);
		tp["lenRunsMedian"] = len_runs_median(alt_seq);

		// Collision tests
		vector<int> col_seq = find_collisions();
		tp["avgCollision"] = avg_collision(col_seq);
		tp["maxCollision"] = max_collision(col_seq);

		// Periodicity tests
		tp["periodicity(1)"] = periodicity(1);
		tp["periodicity(2)"] = periodicity(2);
		tp["periodicity(8)"] = periodicity(8);
		tp["periodicity(16)"] = periodicity(16);
		tp["periodicity(32)"] = periodicity(32);

		// Covariance tests
		tp["covariance(1)"] = covariance(1);
		tp["covariance(2)"] = covariance(2);
		tp["covariance(8)"] = covariance(8);
		tp["covariance(16)"] = covariance(16);
		tp["covariance(32)"] = covariance(32);

		// Compression test
		tp["compression"] = compression();

		for(int j = 0; j < num_tests; j++){
			if(tp[test_names[j]] > t[test_names[j]]){
				C[j][0]++;
			}else if(tp[test_names[j]] == tp[test_names[j]]){
				C[j][1]++;
			}
		}
	}

	#ifdef VERBOSE
	cout << endl;
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
	#endif

	// for(int i = 0; i < num_tests; i++){
	// 	if((C[i][0] + C[i][1] <= 5) || C[i][0] >= PERMS-5){
	// 		exit(-1);
	// 	}
	// }

	#ifdef VERBOSE
	for(int i = 0; i < num_tests; i++)
		cout << test_names[i] << ": " << t[test_names[i]] << endl;
	#endif

	cout << endl;

	return 0;
}
