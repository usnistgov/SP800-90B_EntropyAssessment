#pragma once
#include "../shared/utils.h"

#define D_MMC 16
#define MAX_ENTRIES 100000

// Section 6.3.9 - MultiMMC Prediction Estimate
double multi_mmc_test(byte *data, long len, int alph_size, const bool verbose){
	int winner, cur_winner;
	int entries[D_MMC];
	long i, d, N, C, run_len, max_run_len;
	long scoreboard[D_MMC] = {0};
	array<byte, D_MMC> x;
	bool found_x;
	// d             x                 y          M_d[x,y]
	array<map<array<byte, D_MMC>, map<byte, long>>, D_MMC> M;

	if(len < 3){	
		printf("\t*** Warning: not enough samples to run multiMMC test (need more than %d) ***\n", 3);
		return -1.0;
	}

	N = len-2;
	winner = 0;
	C = 0;
	run_len = 0;
	max_run_len = 0;

	// initialize MMC counts
	memset(x.data(), 0, D_MMC);
	for(d = 0; d < D_MMC; d++){
		if(d < N){
			memcpy(x.data(), data, d+1);
			M[d][x][data[d+1]] = 1;
		}
		entries[d] = 1;
	}

	// perform predictions
	for (i = 2; i < len; i++){
		cur_winner = winner;
		memset(x.data(), 0, D_MMC);
		for(d = 0; d < D_MMC; d++){
			if(i-2 >= d){
				// check if x has been previously seen. If (x,y) has not occurred, 
				// then do not make a prediction for current d and larger d's 
				// as well, since it will not occur for them either. In other words,
				// prediction is NULL, so do not update the scoreboard.
				if((d == 0) || found_x){
					memcpy(x.data(), data+i-d-1, d+1);
					if(M[d].find(x) == M[d].end()) found_x = false;
					else found_x = true;
				}

				if(found_x){
					// x has occurred, find max (x,y) pair across all y's
					if(max_map(M[d][x]) == data[i]){
						// prediction is correct, udpate scoreboard and winner
						if(++scoreboard[d] >= scoreboard[winner]) winner = d;
						if(d == cur_winner){
							C++;
							if(++run_len > max_run_len) max_run_len = run_len;
						}
					}
					else if(d == cur_winner) run_len = 0;
					M[d][x][data[i]]++;
				}
				else if(entries[d] < MAX_ENTRIES){
					memcpy(x.data(), data+i-d-1, d+1);
					M[d][x][data[i]] = 1;
					entries[d]++;
				}
			}
		}
	}

	return(predictionEstimate(C, N, max_run_len, alph_size, "MultiMMC", verbose));
}
