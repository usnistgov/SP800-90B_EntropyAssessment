#pragma once

#include <utility>

#include "utils.h"

/*
* ---------------------------------------------
* 		  HELPER FUNCTIONS / VARIABLES
* ---------------------------------------------
*/

double critical_value[] = {10.828,13.816,16.266,18.467,20.515,22.458,24.322,26.125,27.877,29.588,31.264,32.91,34.528,36.123,37.697,39.252,40.79,42.312,43.82,45.315,46.797,48.268,49.728,51.179,52.62,54.052,55.476,56.892,58.301,59.703,61.098,62.487,63.87,65.247,66.619,67.985,69.347,70.703,72.055,73.402,74.745,76.084,77.419,78.75,80.077,81.4,82.72,84.037,85.351,86.661,87.968,89.272,90.573,91.872,93.168,94.461,95.751,97.039,98.324,99.607,100.888,102.166,103.442,104.716,105.988,107.258,108.526,109.791,111.055,112.317,113.577,114.835,116.092,117.346,118.599,119.85,121.1,122.348,123.594,124.839,126.083,127.324,128.565,129.804,131.041,132.277,133.512,134.746,135.978,137.208,138.438,139.666,140.893,142.119,143.344,144.567,145.789,147.01,148.23,149.449};

double calc_chi_square_cutoff(const int df){
	double x_p = 3.090;
	double h60 = 0.0048;
	double h_v = (60.0/df) * h60;
	double term = 2.0/(9.0 * df);

	return (df * pow(1.0 - term + (x_p - h_v) * sqrt(term), 3));
}

double chi_square_cutoff(int df){
	if(df < 101){
		return critical_value[df];
	}else{
		return calc_chi_square_cutoff(df);
	}
}

void calc_proportions(const byte data[], vector<double> &p){
	
	for(int i = 0; i < SIZE; i++){
		p[data[i]] += (1.0 /SIZE);
	}
}

void calc_expectations(const vector<double> &p, vector<pair<double, pair<byte, byte>>> &e){
	
	for(int i = 0; i < 256; i++){
		for(int j = 0; j < 256; j++){
			e.push_back(pair<double, pair<byte, byte>>(p[i]*p[j]*(SIZE-1), pair<byte, byte>(i, j)));
		}
	}
}

void allocate_bins(const vector<pair<double, pair<byte, byte>>> &e, vector<vector<pair<byte, byte>>> &bins, vector<double> &bin_value){
	
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
}

void check_last_bin(vector<vector<pair<byte, byte>>> &bins, vector<double> &bin_value){
	
	if(bin_value.back() < 5){

		// Pop back the last two values
		double last_val = bin_value.back();
		bin_value.pop_back();

		double second_last_val = bin_value.back();
		bin_value.pop_back();

		// Increment values and add back to the stack
		second_last_val += last_val;
		bin_value.push_back(second_last_val);

		// Pop back the last two vector of pairs that form the last bin
		vector<pair<byte, byte>> last_bin = bins.back();
		bins.pop_back();

		vector<pair<byte, byte>> second_last_bin = bins.back();
		bins.pop_back();

		// Increment last bin and add back to the stack
		for(int i = 0; i < last_bin.size(); i++){
			second_last_bin.push_back(last_bin[i]);
		}

		bins.push_back(second_last_bin);
	}
}

void calc_observed(const byte data[], map<pair<byte, byte>, int> &o){
	
	map_init(o);
	for(int i = 0; i < SIZE-1; i++){
		o[pair<byte, byte>(data[i], data[i+1])]++;
	}
}

void calc_T(const vector<vector<pair<byte, byte>>> &bins, const vector<double> &bin_value, const map<pair<byte, byte>, int> &o, double &T){

	for(int i = 0; i < bins.size(); i++){

		// Calculate observed value for the bin (sum the observed values for each element in the bin)
		double observed_value = 0.0;
		for(int j = 0; j < bins[i].size(); j++){
			observed_value += o.at(bins[i][j]);
		}

		// Increment T
		T += pow((observed_value - bin_value[i]), 2) / bin_value[i];
	}
}

/*
* ---------------------------------------------
* 		  			 TESTS
* ---------------------------------------------
*/

void binary_chi_square_independence(byte data[], double &score, int &df){

}

void chi_square_independence(byte data[], double &score, int &df){

	// Proportion of each element to the entire set
	vector<double> p(256, 0.0);
	calc_proportions(data, p);

	// Calculate the expected number of occurances for each possible pair of values
	vector<pair<double, pair<byte, byte>>> e;
	calc_expectations(p, e);
	sort(e.begin(), e.end());

	// Allocate sorted expected values into bins and accumulate corresponding expected values of entire bins
	vector<vector<pair<byte, byte>>> bins;
	vector<double> bin_value(1, 0.0);
	allocate_bins(e, bins, bin_value);
	
	// Check the last bin to see if it needs to be combined with the previous
	check_last_bin(bins, bin_value);

	// Calculate observed frequency of each pair
	map<pair<byte, byte>, int> o;
	calc_observed(data, o);

	// Calcualte T 
	double T = 0.0;
	calc_T(bins, bin_value, o, T);

	// Return score and degrees of freedom
	score = T;
	df = bins.size()-1;
}

void binary_goodness_of_fit(byte data[10][SIZE/10], double &score, int &df){

}

void goodness_of_fit(byte data[10][SIZE/10], double &score, int &df){

	int sublength = SIZE/10;

	// Get the expected number of each symbol in each subset
	map<byte, double> e;
	map_init(e);

	for(int i = 0; i < 10; i++){
		for(int j = 0; j < sublength; j++){
			e[data[i][j]] += .1; // would be 1 but we divide by 10 later anyways
		}
	}
	
	// Sort the expected values
	vector<pair<double, byte>> e_sorted;
	map<byte, double>::iterator e_itr;
	for(e_itr = e.begin(); e_itr != e.end(); ++e_itr){
		e_sorted.push_back(pair<double, byte>(e_itr->second, e_itr->first));
	}

	sort(e_sorted.begin(), e_sorted.end());

	/*
	* Allocate the sample values into bins. Assign consecutive values 
	* to a bin until the sum of the e_i for those binned items is at
	* least five, then begin assigning the following values to the 
	* next bin. If the expected value of the last bin is less than 
	* five, merge the last two bins. Let q be the number of bins
    * constructed after this procedure.
	*/

    vector<vector<byte>> bins;
    vector<double> E(1, 0.0);
    bool first = true;
    for(int i = 0; i < e_sorted.size(); i++){

    	vector<byte> bin;
    	pair<double, byte> e_pair = e_sorted[i];

    	// Check the previous value 
		// If the value is greater than 5, or the first value
		// Then add it to its own bin
		if(E.back() >= 5.0 || first){

			// Get rid of dummy value
			if(first){
				E.pop_back();
			}

			// Add the value to the bin
			E.push_back(e_pair.first);

			// Add the pair to the bin
			bin.push_back(e_pair.second);
			bins.push_back(bin);

			// Reset flag
			first = false;

		// If the previous value is not greater than 5
		}else{

			// Get previous value and remove from the stack
			double prev_val = E.back();
			E.pop_back();

			// Increment with current value and add it back to the stack
			prev_val += e_pair.first;
			E.push_back(prev_val);

			// Get previous bin and remove from the stack
			bin = bins.back();
			bins.pop_back();

			// Increment with current value and add it back to the stack
			bin.push_back(e_pair.second);
			bins.push_back(bin);
		}
    }

    // If the last bin is less than 5
    // Then combine it with the previous
    if(E.back() < 5.0){

    	// Get previous value and remove it from the stack
    	double last_val = E.back();
    	E.pop_back();

    	// Get second to last value and remove it from the stack
    	double second_last_val = E.back();
    	E.pop_back();

    	// Increment with last value and add it back to the stack
    	second_last_val += last_val;
    	E.push_back(second_last_val);

    	// Get previous bin and remove from the stack
    	vector<byte> last_bin = bins.back();
    	bins.pop_back();

    	// Get second to last bin and remove it from the stack
    	vector<byte> second_last_bin = bins.back();
    	bins.pop_back();

    	// Increment the bin with the content of the last bin
    	for(int i = 0; i < last_bin.size(); i++){
    		second_last_bin.push_back(last_bin[i]);
    	}

    	// Add it back to the stack
		bins.push_back(second_last_bin);
    }

    // Determine observed values from each subset individually
    vector<map<byte, int>> observed;
    map<byte, int> obs_i;
    for(int i = 0; i < 10; i++){

    	// Reset map for this iteration
    	map_init(obs_i);

		for(int j = 0; j < sublength; j++){

			// Increment value
			obs_i[data[i][j]]++;
		}

		// Store in the list
		observed.push_back(obs_i);
	}

	// Begin summing up T
    double T = 0.0;
    for(int i = 0; i < 10; i++){

    	for(int j = 0; j < bins.size(); j++){
    		
    		// Record times a value was observed in a subset 
    		int o = 0;
    		for(int k = 0; k < bins[j].size(); k++){

    			o += observed[i][bins[j][k]];
    		}

    		// Increment T
    		T += pow(o - E[j], 2) / E[j];
    	}
    }

    // Return values
    score = T;
    df = 9*(bins.size()-1);
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

	// Check result to return if test failed
	if(score > cutoff){
		return false;
	}

	// Divide dataset into 10 equal subgroups
	int sublength = SIZE / 10;
	byte data_subsets[10][SIZE/10];

	for(int i = 0; i < 10; i++){
		for(int j = 0; j < sublength; j++){
			data_subsets[i][j] = data[i*sublength + j];
		}
	}

	// Chi Square goodness of fit test
	if(is_binary){
		binary_goodness_of_fit(data_subsets, score, df);
	}else{
		goodness_of_fit(data_subsets, score, df);
	}

	cutoff = chi_square_cutoff(df);

	// Print results
	#ifdef VERBOSE
	cout << "Chi square goodness of fit" << endl;
	cout << "    score = " << score << endl;
	cout << "    degrees of freedom = " << df << endl;
	cout << "    cutoff = " << cutoff << endl;
	#endif

	// Check result to return if test failed
	if(score > cutoff){
		return false;
	}

	return true;
}
