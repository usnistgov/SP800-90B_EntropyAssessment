#pragma once

#include <utility>

#include "utils.h"

double critical_value[] = {10.828,13.816,16.266,18.467,20.515,22.458,24.322,26.125,27.877,29.588,31.264,32.91,34.528,36.123,37.697,39.252,40.79,42.312,43.82,45.315,46.797,48.268,49.728,51.179,52.62,54.052,55.476,56.892,58.301,59.703,61.098,62.487,63.87,65.247,66.619,67.985,69.347,70.703,72.055,73.402,74.745,76.084,77.419,78.75,80.077,81.4,82.72,84.037,85.351,86.661,87.968,89.272,90.573,91.872,93.168,94.461,95.751,97.039,98.324,99.607,100.888,102.166,103.442,104.716,105.988,107.258,108.526,109.791,111.055,112.317,113.577,114.835,116.092,117.346,118.599,119.85,121.1,122.348,123.594,124.839,126.083,127.324,128.565,129.804,131.041,132.277,133.512,134.746,135.978,137.208,138.438,139.666,140.893,142.119,143.344,144.567,145.789,147.01,148.23,149.449};

double calc_chi_square_cutoff(int df){
	double x_p = 3.090;
	double h60 = 0.0048;
	double h_v = (60.0/df) * h60;
	double term = 2.0/(9.0 * df);

	return (df * pow(1.0 - term + (x_p - h_v) * sqrt(term), 3));
}

void binary_chi_square_independence(byte data[], double &score, int &df){

}

void chi_square_independence(byte data[], double &score, int &df){

	/*
	* 5.2.1
	*	1. Find the proportion p_i of each x_i in S (data). Calculate
	*	   the expected number of occurances of each possible pair
	*/

	// Proportion of each element to the entire set
	vector<double> p(256, 0.0);

	// Calculate proportions of the number of each element over the whole
	for(int i = 0; i < SIZE; i++){
		p[data[i]] += (1.0 / SIZE);
	}

	// Expected values of each possible pair of values
	vector<pair<double, pair<byte, byte>>> e;

	// Calculate the expected number of occurances for each possible pair
	for(int i = 0; i < 256; i++){
		for(int j = 0; j < 256; j++){
			pair<byte, byte> key(i, j);
			pair<double, pair<byte, byte>> e_pair(p[i]*p[j]*(SIZE-1), key);
			e.push_back(e_pair);
		}
	}

	sort(e.begin(), e.end());

	/*
	* 5.2.1
	*	2. Allocate the possible pairs, starting from the smallest pair
	* 	   frequency, into bins such that the expected value of each bins
	* 	   is at least 5. The expected value of a bin is equal to the sum
	*      of the values of the pairs that are included in the bin.
	*	   After allocating all pairs, if the expected value of the last
	*	   bin is less than 5, merge the last two bins. Let q be the number
	* 	   of bins constructed after this procedure.
	*/

	// Allocate bins
	vector<vector<pair<byte, byte>>> bins;
	vector<double> bin_value(1, 0.0);

	bool first = true;
	for(int i = 0; i < e.size(); i++){

		pair<double, pair<byte, byte>> current_e = e[i];

		double smallest = current_e.first;
		pair<byte, byte> index = current_e.second;

		// If the previous value is large enough, or if there aren't any values yet
		if(bin_value.back() >= 5.0 || first){

			if(first){
				bin_value.pop_back();
			}

			// Add the value to the first bin
			bin_value.push_back(smallest);

			// Add the pair to the parallel bin
			vector<pair<byte, byte>> bin(1, index);
			bins.push_back(bin);

			first = false;

		// If the previous value is not greater than 5.0, add to it
		}else{
			
			// Add the value to the existing bin value
			double val = bin_value.back();
			bin_value.pop_back();

			val += smallest;
			bin_value.push_back(val);

			// Add the pair to the parallel bin
			vector<pair<byte, byte>> bin = bins.back();
			bins.pop_back();

			bin.push_back(index);
			bins.push_back(bin);
		}
	}

	// If the last bin is under 5, combine it with the previous
	if(bin_value.back() < 5){

		// Pop back the last value
		double val = bin_value.back();
		bin_value.pop_back();

		// Add it to the new last value
		vector<double>::reverse_iterator rit = bin_value.rbegin();
		*rit += val;
		bin_value.push_back(*rit);

		// Pop back the last vector of pairs that form the last bin
		vector<pair<byte, byte>> bin = bins.back();
		bins.pop_back();

		// Add them to the last value
		vector<vector<pair<byte, byte>>>::reverse_iterator rit2 = bins.rbegin();
		for(int i = 0; i < bin.size(); i++){
			rit2->push_back(bin[i]);
		}
	}

	// Calculate observed frequency of each pair
	map<pair<byte, byte>, int> observed;

	for(int j = 0; j < 2; j++){
		for(int i = 0; i < SIZE-1; i++){

			pair<byte, byte> key(data[i], data[i+1]);
			if(j == 0){
				observed[key] = 0;
			}else{
				observed[key]++;
			}
		}
	}

	// Accumulate T
	double T = 0.0;
	for(int i = 0; i < bins.size(); i++){

		// Get current bin
		vector<pair<byte, byte>> bin = bins[i];
		
		// Get expected value for bin
		double expected_value = bin_value[i];

		// Calculate observed value for the bin (sum the observed values for each element in the bin)
		double observed_value = 0.0;

		for(int j = 0; j < bin.size(); j++){
			observed_value += observed[bin[j]];
		}

		// Increment T
		T += pow((observed_value - expected_value), 2) / expected_value;
	}

	// Return score and degrees of freedom
	score = T;
	df = bins.size()-1;
}

double chi_square_cutoff(int df){
	if(df < 101){
		return critical_value[df];
	}else{
		return calc_chi_square_cutoff(df);
	}
}

bool chi_square_tests(byte data[], double mean, double median, bool is_binary){

	double score = 0.0;
	int df = 0;

	// Chi Square independence test
	if(is_binary){
		binary_chi_square_independence(data, score, df);
	}else{
		chi_square_independence(data, score, df);
	}

	double cutoff = chi_square_cutoff(df);

	// Print results
	#ifdef VERBOSE
	cout << "Chi square independence" << endl;
	cout << "    score = " << score << endl;
	cout << "    degrees of freedom = " << df << endl;
	cout << "    cutoff = " << cutoff << endl;
	#endif

	// Chi Square goodness of fit test


	return true;
}
