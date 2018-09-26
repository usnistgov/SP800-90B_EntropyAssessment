#pragma once
#include <iostream>		// std::cout
#include <string>		// std::string
#include <map>			// std::map
#include <set>			// std::set
#include <string.h>		// strlen
#include <iomanip>		// setw / setfill
#include <stdio.h>
//#include <stdlib.h>
#include <cstdlib>
#include <vector>		// std::vector
#include <time.h>		// time
#include <algorithm>	// std::sort
#include <cmath>		// pow, log2
#include <array>		// std::array
#include <omp.h>		// openmp 4.0 with gcc 4.9
#include <bitset>
#include <mutex>		// std::mutex

#define SWAP(x, y) do { int s = x; x = y; y = s; } while(0)

#define MIN_SIZE 1000000
#define PERMS 10000

typedef unsigned char byte;

typedef struct data_t data_t;

struct data_t{
	int word_size; 		// bits per symbol
	int alph_size; 		// symbol alphabet size
	byte *symbols; 		// data words
	byte *bsymbols; 	// data words as binary string
	long len; 			// number of words in data
	long blen; 			// number of bits in data
};

using namespace std;

void free_data(data_t *dp){
	if(dp->symbols != NULL) free(dp->symbols);
	if((dp->word_size > 1) && (dp->bsymbols != NULL)) free(dp->bsymbols);
} 

// Read in binary file to test
bool read_file(const char *file_path, data_t *dp){
	FILE *file; 
	int mask, j, max_symbols = 1 << dp->word_size;
	int symbol_map_down_table[max_symbols];
	long rc, i;

	file = fopen(file_path, "rb");
	if(!file){
		printf("Error: could not open '%s'\n", file_path);
		return false;
	}

	rc = (long)fseek(file, 0, SEEK_END);
	if(rc < 0){
		printf("Error: fseek failed\n");
		fclose(file);
		return false;
	}

	dp->len = ftell(file);
	if(dp->len < 0){
		printf("Error: ftell failed\n");
		fclose(file);
		return false;
	}

	rewind(file);

	if(dp->len == 0){
		printf("Error: '%s' is empty\n", file_path);
		fclose(file);
		return false;
	}

	dp->symbols = (byte*)malloc(dp->len);
	if(dp->symbols == NULL){
		printf("Error: failure to initialize memory for symbols\n");
		fclose(file);
		return false;
	}

    rc = fread(dp->symbols, sizeof(byte), dp->len, file);
	if(rc != dp->len){
		printf("Error: file read failure\n");
		fclose(file);
		free(dp->symbols);
		dp->symbols = NULL;
		return false;
	}

	fclose(file);

	// create symbols (samples) and check if they need to be mapped down
	dp->alph_size = 0;
	memset(symbol_map_down_table, 0, max_symbols*sizeof(int));
	mask = max_symbols-1;
	for(i = 0; i < dp->len; i++){ 
		dp->symbols[i] &= mask;
		if(symbol_map_down_table[dp->symbols[i]] == 0) symbol_map_down_table[dp->symbols[i]] = 1;
	}

	for(i = 0; i < max_symbols; i++){
		if(symbol_map_down_table[i] != 0) symbol_map_down_table[i] = (byte)dp->alph_size++;
	}

	// create bsymbols (bitstring) using the non-mapped data
	dp->blen = dp->len * dp->word_size;
	if(dp->word_size == 1) dp->bsymbols = dp->symbols;
	else{
		dp->bsymbols = (byte*)malloc(dp->blen);
		if(dp->bsymbols == NULL){
			printf("Error: failure to initialize memory for bsymbols\n");
			free(dp->symbols);
			dp->symbols = NULL;
			return false;
		}

		for(i = 0; i < dp->len; i++){
			for(j = 0; j < dp->word_size; j++){
				dp->bsymbols[i*dp->word_size+j] = (dp->symbols[i] >> (dp->word_size-1-j)) & 0x1;
			}
		}
	}

	// map down symbols if less than 2^bits_per_word unique symbols
	if(dp->alph_size < max_symbols){
		for(i = 0; i < dp->len; i++) dp->symbols[i] = (byte)symbol_map_down_table[dp->symbols[i]];
	} 

	return true;
}

void seed(){
	srand(time(NULL));
}

// Fisher-Yates Fast (in place) shuffle algorithm
void shuffle(byte arr[], const int sample_size) {
	long int r;
	static mutex shuffle_mutex;
	unique_lock<mutex> lock(shuffle_mutex);

	for (long int i = sample_size - 1; i > 0; --i) {

		r = (rand() / (float)RAND_MAX) * (i + 1); 	// Proven faster than using % to cast random values
		SWAP(arr[r], arr[i]);
	}
}

// Quick sum array  // TODO
long int sum(const byte arr[], const int sample_size) {
	long int sum = 0;
	for (long int i = 0; i < sample_size; ++i) {
		sum += arr[i];
	}

	return sum;
}

// Quick sum std::array // TODO
template<size_t LENGTH>
int sum(const array<int, LENGTH> &arr) {
	int sum = 0;
	for (int i = 0; i < LENGTH; ++i) {
		sum += arr[i];
	}

	return sum;
}

// Quick sum vector
template<typename T>
T sum(const vector<T> &v) {
	T sum = 0;
	for (unsigned long int i = 0; i < v.size(); ++i) {
		sum += v[i];
	}

	return sum;
}

// Calculate baseline statistics
// Finds mean, median, and whether or not the data is binary
void calc_stats(const byte data[], double &mean, double &median, const int sample_size, const int alphabet_size) {

	// Calculate mean
	mean = sum(data, sample_size) / (double)sample_size;

	// Sort in a vector for median/min/max
	vector<byte> v(data, data + sample_size);
	sort(v.begin(), v.end());

	long int half = sample_size / 2;
	median = (v[half] + v[half - 1]) / 2.0;
}


// Map initialization for integers
void map_init(map<byte, int> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0;
	}
}

// Map initialization for doubles
void map_init(map<byte, double> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0.0;
	}
}

// Map initialization for pair<byte, byte> to int
void map_init(map<pair<byte, byte>, int> &m) {
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			m[pair<byte, byte>(i, j)] = 0;
		}
	}
}

// Calculates proportions of each value as an index
void calc_proportions(const byte data[], vector<double> &p, const int sample_size) {
	for (int i = 0; i < sample_size; i++) {
		p[data[i]] += (1.0 / sample_size);
	}
}

// Determines the standard deviation of a dataset
double std_dev(const vector<int> x, const double x_mean) {
	double sum = 0.0;

	for (unsigned int i = 0; i < x.size(); i++) {
		sum += pow(x[i] - x_mean, 2);
	}

	return sqrt(sum / x.size());
}

// Quick formula for n choose 2 (which can be simplified to [n^2 - n] / 2)
long int n_choose_2(const long int n) {
	return ((n*n) - n) / 2;
}

vector<byte> substr(const byte text[], const int pos, const int len, const int sample_size) {
	int substr_len = len;

	if (pos + len > sample_size) {
		substr_len = sample_size - pos;
	}

	vector<byte> substring;

	for (int i = 0; i < substr_len; i++) {
		substring.push_back(text[pos + i]);
	}

	return substring;
}

// Fast substring with no bounds checking
array<byte, 16> fast_substr(const byte text[], const int pos, const int len) {
	array<byte, 16> substring = { 0 };

	for (int i = 0; i < len; i++) {
		substring[i] = text[pos + i];
	}

	return substring;
}

// Return the key that leads to the maximum value
byte max_map(const map<byte, long>& m){
	long max_cnt = 0;
	byte max_key;
	map<byte, long>::const_iterator itr;

	for(itr = m.begin(); itr != m.end(); ++itr){
		if((itr->second > max_cnt) || ((itr->second == max_cnt) && (itr->first > max_key))){
			max_key = itr->first;
			max_cnt = itr->second;
		}
	}

	return max_key;
}

template<typename T>
T max_vector(const vector<T> &vals) {
	T max = vals[0];
	for (unsigned int i = 0; i < vals.size(); i++) {
		if (vals[i] > max) {
			max = vals[i];
		}
	}

	return max;
}

template<typename T>
T max_arr(const T* vals, const unsigned int k){
	T max = vals[0];
	for (unsigned int i = 0; i < k; i++){
		if (vals[i] > max) {
			max = vals[i];
		}
	}

	return max;
}

double divide(const int a, const int b) {
	return ((double)a / (double)b);
}

double calc_p_global(long C, long N){
	double p = C/(double)N;

	if(p > 0) p = min(1.0, p + 2.576*sqrt((p*(1.0-p))/(N-1.0)));
	else p = 1 - pow(0.01, 1.0/(double)N);
	return p;
}

double calc_p_local(long max_run_len, long N){
	int i;
	double p, q, r, x, p_lo, p_hi, log_alpha, exp, eps;

	// binary search for p_local
	r = (double)max_run_len+1;
	log_alpha = log(0.99);
	eps = 1.0 / (1 << 20); // 2^-20
	p_lo = 0.0 + eps; // avoid division by zero
	p_hi = 1.0 - eps; // avoid division by zero
	do{
		p = (p_lo + p_hi) / 2.0;
		q = 1.0-p;
		
		x = 1.0;
		for(i = 0; i < 10; i++) x = 1.0 + q*pow(p, r)*pow(x, r+1.0);
		exp = log(1.0-p*x) - log((r+1.0-r*x)*q) - (N+1.0)*log(x);

		if(log_alpha < exp) p_lo = p;
		else p_hi = p;
	}while(fabs(p_hi - p_lo) > eps);

	return p;
}
