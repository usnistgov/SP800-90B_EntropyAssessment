#pragma once
#include "../shared/utils.h"

#define NUM_WINS 4

// Section 6.3.7 - Multi Most Common in Window (MCW) Prediction Estimate
double multi_mcw_test(uint8_t *data, long len, int alph_size, const int verbose, const char *label){
	int winner;
	int W[NUM_WINS] = {63, 255, 1023, 4095};
	long i, j, k, N, C, run_len, max_run_len, max_pos; 
	long scoreboard[NUM_WINS] = {0};
	long max_cnts[NUM_WINS] = {0};
	long win_cnts[NUM_WINS][alph_size], win_poses[NUM_WINS][alph_size]; 
	uint8_t frequent[NUM_WINS];
	
	if(len < W[NUM_WINS-1]+1){	
		printf("\t*** Warning: not enough samples to run multiMCW test (need more than %d) ***\n", W[NUM_WINS-1]+1);
		return -1.0;
	}

	N = len-W[0];
	winner = 0;
	C = 0;
	run_len = 0;
	max_run_len = 0;
	for(i = 0; i < NUM_WINS; i++){
		for(j = 0; j < alph_size; j++){
			win_cnts[i][j] = 0;
			win_poses[i][j] = 0;
		}
	}

	// compute initial window counts
	for(i = 0; i < W[NUM_WINS-1]; i++){
		for(j = 0; j < NUM_WINS; j++){
			if(i < W[j]){
				if(max_cnts[j] <= ++win_cnts[j][data[i]]){
					max_cnts[j] = win_cnts[j][data[i]];
					frequent[j] = data[i];
				}
				win_poses[j][data[i]] = i;
			}
		}
	}

	// perform predictions
	for (i = W[0]; i < len; i++){
		// test prediction of winner
		if(frequent[winner] == data[i]){
			C++;
			if(++run_len > max_run_len) max_run_len = run_len;
		}
		else run_len = 0;

		// update scoreboard and select new winner
		for(j = 0; j < NUM_WINS; j++){
			if((i >= W[j]) && (frequent[j] == data[i])){
				if(++scoreboard[j] >= scoreboard[winner]) winner = j;
			}
		}
	
		// update window counts and select new frequents
		for(j = 0; j < NUM_WINS; j++){
			if(i >= W[j]){
				win_cnts[j][data[i-W[j]]]--;
				win_cnts[j][data[i]]++;
				win_poses[j][data[i]] = i;
				if((data[i-W[j]] != frequent[j]) && (max_cnts[j] <= win_cnts[j][data[i]])){
					max_cnts[j] = win_cnts[j][data[i]];
					frequent[j] = data[i];
				}
				else if(data[i-W[j]] == frequent[j]){
					max_cnts[j]--;
					// search for possible new frequent
					max_pos = i-W[j];
					for(k = 0; k < alph_size; k++){
						if((max_cnts[j] < win_cnts[j][k]) || ((max_cnts[j] == win_cnts[j][k]) && (max_pos <= win_poses[j][k]))){
							max_cnts[j] = win_cnts[j][k];
							frequent[j] = k;
							max_pos = win_poses[j][k];
						}
					}
				}
			}
		}
	}

	return(predictionEstimate(C, N, max_run_len, alph_size, "MultiMCW", verbose, label));
}
