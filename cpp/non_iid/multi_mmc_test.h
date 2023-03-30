#pragma once
#include "../shared/utils.h"

#define D_MMC 16
#define MAX_ENTRIES 100000

static double binaryMultiMMCPredictionEstimate(const uint8_t *S, long L, const int verbose, const char *label)
{

   long scoreboard[D_MMC] = {0};
   long *binaryDict[D_MMC];
   long winner = 0;
   long curWinner;
   long curRunOfCorrects = 0;
   long maxRunOfCorrects = 0;
   long correctCount = 0;
   long j, d, i;
   uint32_t curPattern=0;
   long dictElems[D_MMC] = {0};

   assert(L>3);
   assert(D_MMC < 31); //D+1 < 32 to make the bit shifts well defined

   //Initialize the predictors
   for(j=0; j< D_MMC; j++) {
      //For a length m prefix, we need 2^m sets of length 2 arrays.
      //Here, j+1 is the length of the prefix, so we need 2^(j+1) prefixes, or 2*2^(j+1) = 2^(j+2) storage total.
      //Note: 2^(j+2) = 1<<(j+2).
      binaryDict[j] = new long[1U<<(j+2)];
      memset(binaryDict[j], 0, sizeof(long)*(1U<<(j+2)));
   }

   // initialize MMC counts
   for(d=0; d<D_MMC; d++) {
      curPattern = ((curPattern << 1) | (S[d]&1));

      //This is necessarily the first symbol of this length
      (BINARYDICTLOC(d+1, curPattern))[S[d+1]&1] = 1;
      dictElems[d] = 1;
   }

   
   //In C, arrays are 0 indexed.
   //i is the index of the new symbol to be predicted
   for(i=2; i<L; i++) {
      bool found_x = false;

      curWinner = winner;
      curPattern = 0;

      //d+1 is the number of symbols used by the predictor
      for(d=0; (d<D_MMC) && (d<=i-2); d++) {
         uint8_t curPrediction = 2;
         long curCount;
         long *binaryDictEntry;

         //Add S[i-d-1] to the start of curPattern
         curPattern = (curPattern | (((uint32_t)(S[i-d-1]&1))<<d));
         //curPattern should contain the d-tuple (S[i-d-1] ... S[i-1])

         binaryDictEntry = BINARYDICTLOC(d+1, curPattern);

         // check if the prefix x has been previously seen. If the prefix x has not occurred,
         // then do not make a prediction for current d and larger d's
         // as well, since it will not occur for them either. In other words,
         // prediction is NULL, so do not update the scoreboard.
         // Note that found_x is meaningless on the first round, but for that round d==0.
         // All future rounds use a meaningful found_x
         if((d == 0) || found_x) {
            //For the prediction, curPrediction is the max across all pairs (there are only 2 symbols here!)
            if((binaryDictEntry[0] > binaryDictEntry[1])) {
               curPrediction = 0;
               curCount = binaryDictEntry[0];
            } else {
               curPrediction = 1;
               curCount = binaryDictEntry[1];
            }

            if(curCount == 0) found_x = false;
            else found_x = true;
         }

         if(found_x) {
            // x is present as a prefix.
            // Check to see if the current prediction is correct.
            if(curPrediction == S[i]) {
               // prediction is correct, update scoreboard and (the next round's) winner
               scoreboard[d]++;
               if(scoreboard[d] >= scoreboard[winner]) winner = d;

               //If the best predictor was previously d, increment the relevant counters
               if(d == curWinner){
                  correctCount++;
                  curRunOfCorrects++;
                  if(curRunOfCorrects > maxRunOfCorrects) maxRunOfCorrects = curRunOfCorrects;
               }
            } else if(d == curWinner) {
               //This prediction was wrong;
               //If the best predictor was previously d, zero the run length counter
               curRunOfCorrects = 0;
            }

            //Now check to see in (x,y) needs to be counted or (x,y) added to the dictionary
            if(binaryDictEntry[S[i]&1] != 0) {
               //The (x,y) tuple has already been encountered.
               //Increment the existing entry
               binaryDictEntry[S[i]&1]++;
            } else if(dictElems[d] < MAX_ENTRIES) {
               //The x prefix has been encountered, but not (x,y)
               //We're allowed to make a new entry. Do so.
               binaryDictEntry[S[i]&1]=1;
               dictElems[d]++;
            }
         } else if(dictElems[d] < MAX_ENTRIES) {
            //We didn't find the x prefix, so (x,y) surely can't have occurred.
            //We're allowed to make a new entry. Do so.
            binaryDictEntry[S[i]&1]=1;
            dictElems[d]++;
         }
      }
   }

   for(j=0; j<D_MMC; j++) {
      delete[](binaryDict[j]);
      binaryDict[j] = NULL;
   }

   return(predictionEstimate(correctCount, L-2, maxRunOfCorrects, 2, "MultiMMC", verbose, label));
}

// Section 6.3.9 - MultiMMC Prediction Estimate
/* This implementation of the MultiMMC test is a based on NIST's really cleaver implementation,
 * which interleaves the predictions and updates. This makes optimization much easier.
 * It's opaque why it works correctly (in particular, the first few symbols are added in a different
 * order than in the reference implementation), but once the initialization is performed, the rest of the
 * operations are done in the correct order.
 * The general observations that explains why this approach works are that
 * 1) each prediction that could succeed (i.e., ignoring some of the early predictions that must fail due
 *    to lack of strings of the queried length) must occur only after all the correct (x,y) tuples for that
 *    length have been processed. One is free to reorder otherwise.
 * 2) If there is a distinct string of length n, then this induces corresponding unique strings of all
 *    lengths greater than n. We track all string lengths independently (thus conceptually, we
 *    could run out of a short-string length prior to a long string length, thus erroneously not add
 *    some long string to the dictionary after no longer looking for a string to the dictionary when
 *    we should have), this can't happen in practice because we add strings from shortest to longest.
 */
double multi_mmc_test(uint8_t *data, long len, int alph_size, const int verbose, const char *label){
	int winner, cur_winner;
	int entries[D_MMC];
	long i, d, N, C, run_len, max_run_len;
	long scoreboard[D_MMC] = {0};
	array<uint8_t, D_MMC> x;

	if(alph_size == 2) return binaryMultiMMCPredictionEstimate(data, len, verbose, label);

	array<map<array<uint8_t, D_MMC>, PostfixDictionary>, D_MMC> M;

	if(len < 3){	
		printf("\t*** Warning: not enough samples to run multiMMC test (need more than %d) ***\n", 3);
		return -1.0;
	}

	//Step 1
	N = len-2;

	//Step 3
	//scoreboard is initilized above.
	winner = 0;
	
	C = 0;
	run_len = 0;
	max_run_len = 0;

	// initialize MMC counts
	// this performs step 4.a and 4.b for the () case
	memset(x.data(), 0, D_MMC);
	for(d = 0; d < D_MMC; d++){
		if(d < N){
			memcpy(x.data(), data, d+1);
			(M[d][x]).incrementPostfix(data[d+1], true);
			entries[d] = 1;
		}
	}

	// perform predictions
	//i is the index of the new symbol to be predicted
	for (i = 2; i < len; i++){
		bool found_x = false;
		cur_winner = winner;
		memset(x.data(), 0, D_MMC);

		for(d = 0; (d < D_MMC) && (i-2 >= d); d++) {
			map<array<uint8_t, D_MMC>, PostfixDictionary>::iterator curp;
			// check if x has been previously seen as a prefix. If the prefix x has not occurred,
			// then do not make a prediction for current d and larger d's
			// as well, since it will not occur for them either. In other words,
			// prediction is NULL, so do not update the scoreboard.
			// Note that found_x is uninitialized on the first round, but for that round d==0.
			if((d == 0) || found_x) {
				//Get the prediction
				//predict S[i] by using the prior d+1 symbols and the current state
				//We need the d-tuple prior to S[i], that is (S[i-d-1], ..., S[i-1])

				//This populates the curp for the later increment

				memcpy(x.data(), data+i-d-1, d+1);
				curp = M[d].find(x);
				if(curp == M[d].end()) found_x = false;
				else found_x = true;
			}

			if(found_x){
				long predictCount;
				// x has occurred, find max (x,y) pair across all y's
				// Check to see if the current prediction is correct.
				if((curp->second).predict(predictCount) == data[i]){
					// prediction is correct, update scoreboard and winner
					if(++scoreboard[d] >= scoreboard[winner]) winner = d;
					if(d == cur_winner){
						C++;
						if(++run_len > max_run_len) max_run_len = run_len;
					}
				}
				else if(d == cur_winner) {
					//This prediction was wrong;
					//If the best predictor was previously d, zero the run length counter
					run_len = 0;
				}

				//Now check to see in (x,y) needs to be counted or (x,y) added to the dictionary
				if((curp->second).incrementPostfix(data[i], entries[d] < MAX_ENTRIES)) {
					//We had to make a new entry. Count this.
					entries[d]++;
				}
			} else if(entries[d] < MAX_ENTRIES) {
				//We didn't find the x prefix, so (x,y) surely can't have occurred.
				//We're allowed to make a new entry. Do so.
				//curp isn't populated here, because it wasn't found
				memcpy(x.data(), data+i-d-1, d+1);
				(M[d][x]).incrementPostfix(data[i], true);
				entries[d]++;
			}
		}
	}

	return(predictionEstimate(C, N, max_run_len, alph_size, "MultiMMC", verbose, label));
}
