#pragma once

#include <iostream>		// std::cout
#include <string>		// std::string
#include <map>			// std::map
#include <set>			// std::set
#include <string.h>		// strlen
#include <iomanip>		// setw / setfill
#include <stdio.h>
#include <stdlib.h>
#include <vector>		// std::vector
#include <time.h>		// time
#include <algorithm>	// std::sort
#include <math.h>		// pow, log2
#include <array>		// std::array

#define SWAP(x, y) do { int s = x; x = y; y = s; } while(0)

#define SIZE 1000000
#define PERMS 10000

typedef unsigned char byte;

using namespace std;

// Read in binary file to test
void read_file(const char* file_path, byte data[]){
	FILE* file = NULL;

	#ifdef VERBOSE
	cout << "Opening file: " << file_path << endl;
	#endif

	file = fopen(file_path, "rb");
	auto rc = fread(data, 1, SIZE, file);
	fclose(file);

	#ifdef VERBOSE
	cout << "Data read complete" << endl;
	#endif
}

// Fisher-Yates Fast (in place) shuffle algorithm
void shuffle(byte arr[]){
	srand(time(NULL));
	long int r;

	for(long int i = SIZE-1; i > 0; --i){
		r = (rand() / (float)RAND_MAX) * (i+1); 	// Proven faster than using % to cast random values
		SWAP(arr[r], arr[i]);
	}
}

// Quick sum array
long int sum(const byte arr[]){
	long int sum = 0;
	for(long int i = 0; i < SIZE; ++i){
		sum += arr[i];
	}

	return sum;
}

// Quick sum std::array
template<size_t LENGTH>
int sum(const array<int, LENGTH> &arr){
	int sum = 0;
	for(int i = 0; i < LENGTH; ++i){
		sum += arr[i];
	}

	return sum;
}

// Quick sum vector
template<typename T>
long double sum(const vector<T> &v){
	long double sum = 0;
	for(unsigned long int i = 0; i < v.size(); ++i){
		sum += v[i];
	}

	return sum;
}

// Calculate baseline statistics
// Finds mean, median, and whether or not the data is binary
void calc_stats(const byte data[], double &mean, double &median, bool &is_binary){

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

// Map initialization for integers
void map_init(map<byte, int> &m){
	for(int i = 0; i < 256; i++){
		m[i] = 0;
	}
}

// Map initialization for doubles
void map_init(map<byte, double> &m){
	for(int i = 0; i < 256; i++){
		m[i] = 0.0;
	}
}

// Map initialization for pair<byte, byte> to int
void map_init(map<pair<byte, byte>, int> &m){
	for(int i = 0; i < 256; i++){
		for(int j = 0; j < 256; j++){
			m[pair<byte, byte>(i,j)] = 0;
		}
	}
}

// Calculates proportions of each value as an index
void calc_proportions(const byte data[], vector<double> &p){
	for(int i = 0; i < SIZE; i++){
		p[data[i]] += (1.0 / SIZE);
	}
}

// Determines the standard deviation of a dataset
double std_dev(const vector<int> x, const double x_mean){
	double sum = 0.0;

	for(int i = 0; i < x.size(); i++){
		sum += pow(x[i] - x_mean, 2);
	}

	return sqrt(sum / x.size());
}

// Quick formula for n choose 2 (which can be simplified to [n^2 - n] / 2)
long int n_choose_2(const long int n){
	return ((n*n) - n) / 2;
}

vector<byte> substr(const byte text[], const int pos, const int len){
	int substr_len = len;

	if(pos+len > SIZE){
		substr_len = SIZE - pos;
	}

	vector<byte> substring;

	for(int i = 0; i < substr_len; i++){
		substring.push_back(text[pos+i]);
	}

	return substring;
}

// Fast substring with no bounds checking
array<byte, 16> fast_substr(const byte text[], const int pos, const int len){
	array<byte, 16> substring;

	for(int i = 0; i < len; i++){
		substring[i] = text[pos+i];
	}

	return substring;
}

// Return the key that leads to the maximum value
byte max_map(const map<byte, int> m){
	int max = -1;
	byte key;

	map<byte, int>::const_iterator itr;
	for(itr = m.begin(); itr != m.end(); ++itr){
		if(itr->second > max){
			key = itr->first;
			max = itr->second;
		}
	}

	return key;
}

template<typename T>
T max_vector(const vector<T> &vals){
	T max = vals[0];
	for(int i = 0; i < vals.size(); i++){
		if(vals[i] > max){
			max = vals[i];
		}
	}

	return max;
}

template<typename T>
T max_arr(const T* vals, const unsigned int k){
	T max = vals[0];
	for(int i = 0; i < k; i++){
		if(vals[i] > max){
			max = vals[i];
		}
	}

	return max;
}

double divide(const int a, const int b){
	return ((double)a / (double)b);
}

double calc_p_avg(const int C, const int N){
	double p_global = ((double)C / (double)N);
	return p_global + 2.576*sqrt((p_global - p_global*p_global)/(double)N);
}

double find_root(const double p, const int r){
	double q = 1-p;
	double s = 1;

	for(int i = 0; i < 10; i++){
		s = 1+q*(pow(p, r)*pow(s, r+1));
	}

	return s;
}

double calc_qn(const double p, const int r, const int N){
	double q = 1-p;
	double x = find_root(p, r);

	double qn = (1-p*x) / ((r+1-r*x)*q);
	qn /= pow(x, N+1);

	return qn;
}

int find_max_run(const vector<int> &correct){
	int run = 0; 
	int max_run = 0;
	int prev = -1;

	for(int i = 0; i < correct.size(); i++){
		run = (correct[i] ? run+1 : 0);
		
		if(run > max_run){
			max_run = run;
		}

		prev = correct[i];
	}

	return max_run;
}

double calc_run(const vector<int> &correct){
	int N = correct.size();
	double alpha = .99;
	int r = find_max_run(correct);

	// Do a binary search for p
	double p = .5;
	double adj = .5;
	double qn = calc_qn(p, r+1, N);

	for(int i = 0; i < 20; i++){
		adj /= 2;
		if(qn > alpha){
			p += adj;
		}else{
			p -= adj;
		}

		// Find probability there is no run of length r+1
		qn = calc_qn(p, r+1, N);
		if(abs(qn-alpha) <= .0001) break;
	}

	return p;
}