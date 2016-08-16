#include "../shared/utils.h"

int d = 1000;
int v = SIZE - d;

double G(const double p){
	
	double G_sum = 0.0;
	double in_sum = 0.0;
	for(int i = 1; i < d+1; i++){
		in_sum += log2(i) * pow(1.0-p, i-1);
	}

	G_sum = pow(p, 2) * in_sum * v;

	vector<double> st;
	for(int i = d+1; i < SIZE+1; i++){
		st.push_back(log2(i) * pow(1.0-p, i-1));
	}

	in_sum = 0;
	for(unsigned int i = 0; i < st.size(); i++){
		in_sum += ((SIZE-i-(d+1)) * st[i]);
	}

	G_sum += (pow(p, 2) * in_sum);

	G_sum += (p * sum(st));

	return (G_sum / (double)v);
}

double EppM(const double p, const int k){
	double q = (1.0 - p) / (k - 1.0);
	return (G(p) + ((k - 1) * G(q)));
}

bool solve_for_p(const double mu_bar, const int k, double &p){

	double min_p = 1.0 / (double)k;
	p = (1 - min_p) / 2.0 + min_p;
	double adj = 1 - min_p;

	double E_p_maxvalid = EppM(1.0 / (double)k, k);
	if(mu_bar > E_p_maxvalid){
		return false;
	}

	double Ep = EppM(p, k);
	while(abs(mu_bar - Ep) > 0.00001){
		adj /= 2.0;
		
		if(mu_bar > Ep){
			p -= adj;
		}else{
			p += adj;
		}

		Ep = EppM(p, k);
	}

	return true;
}

double compression_test(const byte data[]){
	
	// Step 1-3
	//   Partition dataset into two disjoint groups of size d and v
	//   Create a dictionary from the first d observations
	//   Initialize as 0, then have each byte correspond to the last occurance
	//   If dict is nonzero, then assign D and update the dictionary
	//   Otherwise, add value to dictionary and assign D
	vector<int> dict(256, 0);
	vector<double> D;
	for(int i = 0; i < SIZE; i++){
		if(i >= d){
			D.push_back(log2(i-dict[data[i]]));
		}

		dict[data[i]] = i;
	}

	// Step 4
	//   Let b be the number of bits needed to represent the largest possible character
	//   Calculate the mean and standard deviation
	double b = floor(log2(256-1)) + 1;
	double mu = sum(D) / (double)D.size();
	double c = 0.7 - (0.8/b) + ((4 + (32/b) * pow(v, -3.0/b)) / 15);
	double sigma = 0.0;

	for(unsigned int i = 0; i < D.size(); i++){
		sigma += pow(D[i], 2);
	}

	sigma = c * sqrt((sigma / (double)v) - pow(mu, 2));

	// Step 5
	//   Compute lower-bound based on normal distribution
	double mu_bar = mu - (2.576*sigma)/sqrt(v);

	// Step 6
	//   Using a binary search solve for p
	double p = 0.0;
	bool valid = solve_for_p(mu_bar, 256, p);

	// Step 7
	//   Return based on if the search was successful
	return (valid ? -log2(p) : log2(256));
}