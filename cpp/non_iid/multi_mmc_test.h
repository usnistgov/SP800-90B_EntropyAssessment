#pragma once
#include "../shared/utils.h"

#define D_MMC 16
#define MAX_ENTRIES 100000

static double binaryMultiMMCPredictionEstimate(const byte *S, long L, const bool verbose)
{
   byte prediction[D_MMC];
   long predictionCount[D_MMC];
   long scoreboard[D_MMC];
   long *binaryDict[D_MMC];
   long winner;
   long curRunOfCorrects;
   long maxRunOfCorrects;
   long correctCount;
   long j, d, i;
   uint32_t curPattern=0;
   bool makeBranches;
   long dictElems[D_MMC];

   assert(L>3);
   assert(D_MMC < 32);

   //Initialize the predictors
   for(j=0; j< D_MMC; j++) {
      prediction[j] = 0;
      scoreboard[j] = 0;
      predictionCount[j] = 0;
      dictElems[j]=0;
      binaryDict[j] = new long[2*(1U<<(j+1))];
      for(i=0; i<(1U<<(j+1))*2; i++) {
         (binaryDict[j])[i] = 0;
      }
   }

   winner = 0;
   maxRunOfCorrects = 0;
   curRunOfCorrects = 0;
   correctCount = 0;

   curPattern = S[0];

   //In C, arrays are 0 indexed.
   //i is the index of the new symbol to be predicted
   //This is all rather confusing, but it helps to run the first and last few cycles, to verify that it works...
   for(i=2; i<L; i++) {
      //4a

      //currentPattern should contain the D_MMC-tuple (S[i-D_MMC-1] ... S[i-2])

      //d is the number of symbols used by the predictor
      for(d=1; (d<=D_MMC) && (d<=(i-1)); d++) {
         //update the state to reflect last round's new value (add S[i-1] to the predictor)
         //We need the d-tuple prior to S[i-1], that is (S[i-d-1], ..., S[i-2])
         makeBranches = dictElems[d-1] < MAX_ENTRIES;
         //This tuple is stored in curPattern. Take the lower d bits.
         if(incrementBinaryDict(binaryDict, d, curPattern, S[i-1], makeBranches, true) && makeBranches) {
            dictElems[d-1]++;
         }
      }

      //4b
      //Add S[i-1] to the curPattern
      curPattern = ((curPattern << 1) | S[i-1]) & ((1 << D_MMC) - 1);
      //curPattern should contain the D_MMC-tuple (S[i-D_MMC] ... S[i-1])

      for(d=1; d<=D_MMC; d++) {
         //Get the predictions
         //predict S[i] by using the prior d bits and the current state
         //We need the d-tuple prior to S[i], that is (S[i-d], ..., S[i-1])
         if(d <= i) {
            predictionCount[d-1] = predictBinaryDict(binaryDict, d, curPattern, prediction+(d-1));
         } else {
            //We can't form a string to query the predictor
            predictionCount[d-1] = 0;
         }
      }

      if((predictionCount[winner]!=0) && (S[i] == prediction[winner])) {
         correctCount++;
         curRunOfCorrects++;
      } else {
         curRunOfCorrects = 0;
      }

      if(curRunOfCorrects > maxRunOfCorrects) {
         maxRunOfCorrects = curRunOfCorrects;
      }

      for(j=0; j<D_MMC; j++) {
         if(predictionCount[j] != 0) {
            if(prediction[j] == S[i]) {
               scoreboard[j]++;
               if(scoreboard[j] >= scoreboard[winner]) {
                  winner=j;
               }
            }
         }
      }
   }

   for(j=0; j<D_MMC; j++) {
      delete[](binaryDict[j]);
      binaryDict[j] = NULL;
   }

   return(predictionEstimate(correctCount, L-2, maxRunOfCorrects, 2, "MultiMMC", verbose));
}

// Section 6.3.9 - MultiMMC Prediction Estimate
double multi_mmc_test(byte *data, long len, int alph_size, const bool verbose){
	int winner, cur_winner;
	int entries[D_MMC];
	long i, d, N, C, run_len, max_run_len;
	long scoreboard[D_MMC] = {0};
	array<byte, D_MMC> x;
	bool found_x;

	if(alph_size == 2) return binaryMultiMMCPredictionEstimate(data, len, verbose);

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
