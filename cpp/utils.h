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
#include <math.h>		// pow

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
	fread(data, 1, SIZE, file);
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

// Quick sum vector
long int sum(const vector<int> &v){
	long int sum = 0;
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