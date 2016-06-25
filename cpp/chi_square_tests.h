#pragma once

#include "utils.h"

double critical_value[] = {10.828,13.816,16.266,18.467,20.515,22.458,24.322,26.125,27.877,29.588,31.264,32.91,34.528,36.123,37.697,39.252,40.79,42.312,43.82,45.315,46.797,48.268,49.728,51.179,52.62,54.052,55.476,56.892,58.301,59.703,61.098,62.487,63.87,65.247,66.619,67.985,69.347,70.703,72.055,73.402,74.745,76.084,77.419,78.75,80.077,81.4,82.72,84.037,85.351,86.661,87.968,89.272,90.573,91.872,93.168,94.461,95.751,97.039,98.324,99.607,100.888,102.166,103.442,104.716,105.988,107.258,108.526,109.791,111.055,112.317,113.577,114.835,116.092,117.346,118.599,119.85,121.1,122.348,123.594,124.839,126.083,127.324,128.565,129.804,131.041,132.277,133.512,134.746,135.978,137.208,138.438,139.666,140.893,142.119,143.344,144.567,145.789,147.01,148.23,149.449};

double calc_chi_square_cutoff(int df){
	double x_p = 3.090;
	double h60 = 0.0048;
	double h_v = (60.0/df) * h60;
	double term = 2.0/(9.0 * df);

	double chi_sqr = df * pow(1.0 - term + (x_p - h_v) * sqrt(term), 3);
}

void binary_chi_square_independence(byte data[], double &score, int &df){

}

void chi_square_independence(byte data[], double &score, int &df){

}

double chi_square_cutoff(int df){
	if(df < 101){
		return critical_value[df];
	}else{
		return calc_chi_square_cutoff(df);
	}
}

bool chi_square_tests(byte data[], double mean, double median, bool is_binary){

	double score = 0;
	int df = 0;

	if(is_binary){
		binary_chi_square_independence(data, &score, &df);
	}else{
		chi_square_independence(data, &score, &df);
	}

	double cutoff = chi_square_cutoff(df);


	return true;
}