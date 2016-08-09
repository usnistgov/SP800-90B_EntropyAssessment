#include "../shared/utils.h"

double markov_test(const byte data[], const int k, double alpha){
	
	// Step 1
	//   Re-define the confidence level to alpha
	int d = 128;
	int alpha_exp = max((int)pow(k, 2), d);
	alpha = pow(alpha, alpha_exp);

	// Step 2
	//   Estimate the initial state probability distribution P
	int occurance[k] = {0};
	for(int i = 0; i < SIZE; i++){
		occurance[data[i]]++;
	}

	double epsilon_term = log2(1.0 / (1.0-alpha));
	double epsilon = sqrt(epsilon_term / (2*SIZE));

	double P[k];
	for(int i = 0; i < k; i++){
		P[i] = min(1.0, divide(occurance[i], SIZE) + epsilon);
	}

	// Step 3
	//   Remove 1 from the last occurance
	//   Estimate probabilities of the transition matrix T
	occurance[k-1]--;

	int trans_occurance[k][k] = {0};
	int o_i = data[0];
	for(int i = 1; i < SIZE; i++){
		int o_j = data[i];
		trans_occurance[o_i][o_j]++;
		o_i = o_j;
	}

	double epsilon_i[k];
	for(int i = 0; i < k; i++){
		if(occurance[i] == 0){
			epsilon_i[i] = 0.0;
		}else{
			sqrt(epsilon_term / (2*occurance[i]));
		}
	}

	double T[k][k] = {0.0};
	for(int i = 0; i < k; i++){
		for(int j = 0; j < k; j++){
			if(occurance[i] == 0){
				T[i][j] = 1.0;
			}else{
				T[i][j] = min(1.0, divide(trans_occurance[i][j], occurance[i]) + epsilon_i[i]);
			}
		}
	}

	// Step 4
	//   Using transition matrix T to find p_max
	for(int j = 1; j < d; j++){
		double h[k];
		for(int c = 0; c < k; c++){
			double Pp[k];
			for(int i = 0; i < k; i++){
				Pp[i] = P[i]*T[i][c];
			}
			h[c] = max_arr(Pp, k);
		}
		copy(h, h+k, P);
	}

	double p_max = max_arr(P, k);

	return -log2(p_max) / d;
}