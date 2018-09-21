#pragma once

#include <utility>
#include "../shared/utils.h"

/*
* ---------------------------------------------
* 		  HELPER FUNCTIONS / VARIABLES
* ---------------------------------------------
*/

double critical_value[] = {10.828,13.816,16.266,18.467,20.515,22.458,24.322,26.125,27.877,29.588,31.264,32.91,34.528,36.123,37.697,39.252,40.79,42.312,43.82,45.315,46.797,48.268,49.728,51.179,52.62,54.052,55.476,56.892,58.301,59.703,61.098,62.487,63.87,65.247,66.619,67.985,69.347,70.703,72.055,73.402,74.745,76.084,77.419,78.75,80.077,81.4,82.72,84.037,85.351,86.661,87.968,89.272,90.573,91.872,93.168,94.461,95.751,97.039,98.324,99.607,100.888,102.166,103.442,104.716,105.988,107.258,108.526,109.791,111.055,112.317,113.577,114.835,116.092,117.346,118.599,119.85,121.1,122.348,123.594,124.839,126.083,127.324,128.565,129.804,131.041,132.277,133.512,134.746,135.978,137.208,138.438,139.666,140.893,142.119,143.344,144.567,145.789,147.01,148.23,149.449};

// TODO not clear if this needs to be updated
double calc_chi_square_cutoff(const int df){
	double x_p = 3.090;
	double h60 = 0.0048;
	double h_v = (60.0/df) * h60;
	double term = 2.0/(9.0 * df);

	return (df * pow(1.0 - term + (x_p - h_v) * sqrt(term), 3));
}

// TODO not clear if this needs to be updated
double chi_square_cutoff(const int df){
	if(df < 101){
		return critical_value[df-1];
	}else{
		return calc_chi_square_cutoff(df);
	}
}

/*
* ---------------------------------------------
* 	  HELPERS FOR CHI_SQUARE_INDEPENDENCE
* ---------------------------------------------
*/

void calc_expectations(const vector<double> &p, vector<pair<double, pair<byte, byte>>> &e, const int sample_size){
	for(int i = 0; i < p.size(); i++){
		for(int j = 0; j < p.size(); j++){
			double exp = p[i] * p[j] * sample_size * 0.5;
			e.push_back(pair<double, pair<byte, byte>>(exp, pair<byte, byte>(i, j)));
		}
	}
}

void allocate_bins(const vector<pair<double, pair<byte, byte>>> &e, vector<vector<pair<byte, byte>>> &bins, vector<double> &bin_value){
	
	bool first = true;
	for(unsigned int i = 0; i < e.size(); i++){

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
	
	if(bin_value.back() < 5.0){

		// Pop back the last two values
		double last_val = bin_value.back();
		bin_value.pop_back();

		double second_last_val = bin_value.back();
		bin_value.pop_back();

		// Increment values and add back to the stack
		second_last_val += last_val;
		bin_value.push_back(second_last_val);

		// Pop back the last two vectors of pairs
		vector<pair<byte, byte>> last_bin = bins.back();
		bins.pop_back();

		vector<pair<byte, byte>> second_last_bin = bins.back();
		bins.pop_back();

		// Increment last bin and add back to the stack
		for(unsigned int i = 0; i < last_bin.size(); i++){
			second_last_bin.push_back(last_bin[i]);
		}

		bins.push_back(second_last_bin);
	}
}

// TODO slow method
void calc_observed(const byte data[], const vector<vector<pair<byte, byte>>> &bins, vector<int> &o, const int sample_size){

	o = vector<int>(bins.size(), 0);

	for(int j = 0; j < sample_size-1; j+=2){

		// Search each bin for (data_i, data_i+1)
		for(int i = 0; i < bins.size(); i++){

			// If the pair exists in the bin, increment occurances by 1
			for(int k = 0; k < bins[i].size(); k++){

				if (bins[i][k].first == data[j] && bins[i][k].second == data[j+1]){
					o[i]++;

					// Each pair can only be found once, so once we find it, move on to the next pair
					k = bins[i].size();
					i = bins.size();
				}
			}
		}
	}
}

void calc_T(const vector<vector<pair<byte, byte>>> &bins, const vector<double> &bin_value, const vector<int> &o, double &T){

	for (unsigned int i = 0; i < bins.size(); i++){

		T += divide(pow((o[i] - bin_value[i]), 2), bin_value[i]);
	}
}

/*
* ---------------------------------------------
* 		  HELPERS FOR GOODNESS_OF_FIT
* ---------------------------------------------
*/

void calc_expectations(byte** data, vector<pair<double, byte>> &e, const int sample_size, const int alphabet_size){

	int sublength = sample_size / 10;

	// Build map of expected values
	// Using a map for the access structure because we need the byte later as well.
	// If we used a vector where the index was the byte, (as it is 0-255) we lose the unification of the expected value
	// and the byte when we sort later on.
	map<byte, double> e_map;
	map_init(e_map);

	for(int i = 0; i < 10; i++){
		for(int j = 0; j < sublength; j++){
			e_map[data[i][j]] += .1; 		// could be 1 but we would divide by 10 later anyways
		}
	}

	// Put the content of the map into a vector (a quickly sortable structure)
	map<byte, double>::iterator e_itr;
	for(e_itr = e_map.begin(); e_itr != e_map.end(); ++e_itr){
		e.push_back(pair<double, byte>(e_itr->second, e_itr->first));
	}
}

void allocate_bins(const vector<pair<double, byte>> &e, vector<vector<byte>> &bins, vector<double> &bin_value){
	
	bool first = true;
	for(unsigned int i = 0; i < e.size(); i++){

		vector<byte> bin;
		pair<double, byte> e_pair = e[i];

		// If the previous value is greater than 5, or is the first value then add it to its own bin
		if(bin_value.back() >= 5.0 || first){

			// Get rid of dummy value and set the flag
			if(first){
				bin_value.pop_back();
				first = false;
			}

			bin_value.push_back(e_pair.first);

			bin.push_back(e_pair.second);
			bins.push_back(bin);

		// If the previous value is not greater than 5, add the current value to the previous bin
		}else{

			// Get previous value, add to it, and place it back on the stack
			double prev_val = bin_value.back() + e_pair.first;
			bin_value.pop_back();
			bin_value.push_back(prev_val);

			// Get previous bin, add to it, and place it back on the stack
			bin = bins.back();
			bin.push_back(e_pair.second);
			bins.pop_back();
			bins.push_back(bin);
		}
	}
}

void check_last_bin(vector<vector<byte>> &bins, vector<double> &bin_value){

	// If the last bin is less than 5
	// Then combine it with the previous
	if(bin_value.back() < 5.0){

		// Pop back the last two values
		double last_val = bin_value.back();
		bin_value.pop_back();

		double second_last_val = bin_value.back();
		bin_value.pop_back();

		// Increment values and add back to the stack
		second_last_val += last_val;
		bin_value.push_back(second_last_val);

		// Pop back the last two vectors
		vector<byte> last_bin = bins.back();
		bins.pop_back();

		vector<byte> second_last_bin = bins.back();
		bins.pop_back();

		// Increment last bin and add it back to the stack
		for(unsigned int i = 0; i < last_bin.size(); i++){
			second_last_bin.push_back(last_bin[i]);
		}

		bins.push_back(second_last_bin);
	}
}

void calc_observed(byte** data, vector<map<byte, int>> &o, const int sample_size){

	int sublength = sample_size / 10;
	map<byte, int> obs_i;
	for(int i = 0; i < 10; i++){

		// Reset map for this iteration
		map_init(obs_i);

		for(int j = 0; j < sublength; j++){
			obs_i[data[i][j]]++;
		}

		o.push_back(obs_i);
	}
}

void calc_T(const vector<vector<byte>> &bins, const vector<double> &bin_value, const vector<map<byte, int>> &o, double &T){
	
	for(int i = 0; i < 10; i++){
		for(unsigned int j = 0; j < bins.size(); j++){
			
			// Record times a value was observed in a subset 
			int o_i = 0;
			for(unsigned int k = 0; k < bins[j].size(); k++){
				o_i += o[i].at(bins[j][k]);
			}

			// Increment T
			T += pow(o_i - bin_value[j], 2) / bin_value[j];
		}
	}
}

/*
* ---------------------------------------------
* 		  			 TESTS
* ---------------------------------------------
*/

void binary_chi_square_independence(const byte data[], double &score, int &df, const int sample_size){

	// Compute proportion of 0s and 1s
	double p0 = 0.0, p1 = 0.0;

	for(int i = 0; i < sample_size; i++){
		p1 += data[i];
	}

	p1 /= sample_size;
	p0 = 1.0 - p1;

	// Compute m
	double min_p = min(p0, p1);
	int m = 11;
	int threshhold = 5;
	while(m > 1){
		if (pow(min_p, m) * (sample_size / m) >= threshhold){
			break;
		}else{
			m--;
		}
	}

	if (m < 2){
		score = 0.0;
		df = 0;
		return;
	}

	// Test is only run if m >= 2
	double T = 0;

	// Count occurances of m-bit tuples by converting to decimal and using as index in a vector
	vector<int> occ(pow(2, m), 0);
	for(int i = 0; i < occ.size(); i++){

		int decimal = 0;
		for(int j = 0; j < m; j++){
			decimal += (pow(2, j) * data[i*m + j]);
		}

		occ[decimal]++;
	}

	for(unsigned int i = 0; i < occ.size(); i++){

		// GCC only, counts the number of 1s in an integer
		int w = __builtin_popcount(i);

		double e = pow(p1, w) * pow(p0, m - w) * (sample_size / m);

		T += pow(occ[i] - e, 2) / e;
	}

	score = T;
	df = pow(2, m) - 2;
}

void chi_square_independence(const byte data[], double &score, int &df,  const int sample_size, const int alphabet_size){

	// Proportion of each element to the entire set
	vector<double> p(alphabet_size, 0.0);
	calc_proportions(data, p, sample_size);

	//for (int i = 0; i < p.size(); i++){
	//	cout << "p[" << i << "] = " << p[i] << endl;
	//}

	// Calculate the expected number of occurances for each possible pair of symbols
	vector<pair<double, pair<byte, byte>>> e;
	calc_expectations(p, e, sample_size);
	sort(e.begin(), e.end());

	//for(int i = 0; i < e.size(); i++){
	//	cout << "e" << (int)e[i].second.first + 1 << (int)e[i].second.second + 1 << " = " << e[i].first << endl;
	//}

	// Allocate sorted expected values into bins and accumulate corresponding expected values of entire bins
	vector<vector<pair<byte, byte>>> bins;
	vector<double> bin_value(1, 0.0);		// Needs a starting value because of the 5.0 bin value limit
	allocate_bins(e, bins, bin_value);
	
	// Check the last bin to see if it needs to be combined with the previous
	check_last_bin(bins, bin_value);

	//for(int i = 0; i < bins.size(); i++){
	//	if (i % 1000 != 0) continue;
	//	cout << "Bin " << i+1 << ": ";
	//	for(int j = 0; j < bins[i].size(); j++){
	//		cout << "(" << (int)bins[i][j].first + 1 << ", " << (int)bins[i][j].second + 1 << "),";
	//	}
	//	cout << " bin_value: " << bin_value[i] << endl; 
	//}

	//cout << "Calculating observed" << endl;

	// Calculate the observed frequency of each pair of symbols
	vector<int> o;
	calc_observed(data, bins, o, sample_size);

	//for(int i = 0; i < bins.size(); i++){
	//	cout << "Bin observed: " << o[i] << endl;
	//}

	// Calcualte T 
	double T = 0.0;
	calc_T(bins, bin_value, o, T);

	//cout << "T: " << T << endl;

	// Return score and degrees of freedom
	score = T;
	df = bins.size() - alphabet_size;
}

void binary_goodness_of_fit(byte** data, double &score, int &df, const int sample_size){

	// Find proportion of 1s to the whole data set
	int sublength = sample_size / 10;
	int ones = 0;

	for(int i = 0; i < 10; i++){
		for(int j = 0; j < sublength; j++){
			ones += data[i][j];
		}
	}

	double p = divide(ones, sample_size);
	double T = 0;

	// Compute expected 0s and 1s in each sub-sequence
	double e0 = (1.0 - p) * sublength;
	double e1 = p * sublength;

	for(int i = 0; i < 10; i++){

		// Count actual 0s and 1s in each sub-sequence
		int o0 = 0, o1 = 0;

		for(int j = 0; j < sublength; j++){
			o1 += data[i][j];
		}

		o0 = sublength - o1;

		// Compute T
		T += (pow(o0 - e0, 2) / e0) + (pow(o1 - e1, 2) / e1);
	}

	score = T;
	df = 9;
}

void goodness_of_fit(byte** subset_data, double &score, int &df, const int sample_size, const int alphabet_size){

	// Get the expected number of each symbol in each subset and sort in ascending order
	vector<pair<double, byte>> e;
	calc_expectations(subset_data, e, sample_size, alphabet_size);
	sort(e.begin(), e.end());

	// Allocate sorted expected values into bins and accumulate corresponding expected values of entire bins
	vector<vector<byte>> bins;
	vector<double> bin_value(1, 0.0);		// Needs a starting value because of the 5.0 bin value limit
	allocate_bins(e, bins, bin_value);
	
	// Check the last bin to see if it needs to be combined with the previous
	check_last_bin(bins, bin_value);

	// Calculate the observed frequency of each symbol in each subset
	vector<map<byte, int>> o;
	calc_observed(subset_data, o, sample_size);

	// Calculate T
	double T = 0.0;
	calc_T(bins, bin_value, o, T);

	// Return score and degrees of freedom
	score = T;
	df = 9*(bins.size()-1);
}

bool chi_square_tests(const byte data[], const double mean, const double median, const int sample_size, const int alphabet_size, const bool verbose){

	double score = 0.0;
	int df = 0;

	// Chi Square independence test
	if(alphabet_size == 2){
		binary_chi_square_independence(data, score, df, sample_size);
	}else{
		chi_square_independence(data, score, df, sample_size, alphabet_size);
	}

	double cutoff = chi_square_cutoff(df);

	// Print results
	if(verbose){
		printf("Chi square independence\n");
		printf("\tscore = %f\n", score);
		printf("\tdegrees of freedom = %d\n", df);
		printf("\tcutoff = %f\n\n", cutoff);
	}

	// Check result to return if test failed
	if(score > cutoff){
		return false;
	}

	// Reset score and df
	score = 0.0;
	df = 0;

	// Divide dataset into 10 equal subgroups
	int sublength = sample_size / 10;
	byte* data_subsets[10];

	for(int i = 0; i < 10; i++){
		data_subsets[i] = new byte[sublength];
		for(int j = 0; j < sublength; j++){
			data_subsets[i][j] = data[i*sublength + j];
		}
	}

	// Chi Square goodness of fit test
	if(alphabet_size == 2){
		binary_goodness_of_fit(data_subsets, score, df, sample_size);
	}else{
		goodness_of_fit(data_subsets, score, df, sample_size, alphabet_size);
	}

	cutoff = chi_square_cutoff(df);

	// Print results
	if(verbose){
		printf("Chi square goodness of fit\n");
		printf("\tscore = %f\n", score);
		printf("\tdegrees of freedom = %d\n", df);
		printf("\tcutoff = %f\n\n", cutoff);
	}

	// Check result to return if test failed
	if(score > cutoff){
		return false;
	}

	return true;
}
