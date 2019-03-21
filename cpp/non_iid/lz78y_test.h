#pragma once
#include "../shared/utils.h"

#define B 16
#define MAX_DICTIONARY_SIZE 65536

static double binaryLZ78YPredictionEstimate(const byte *S, long L, const bool verbose){
   byte curPrediction;
   byte prediction;
   long maxCount;
   long curPredictionCount;
   long curRunOfCorrects;
   long maxRunOfCorrects;
   long correctCount;
   long j, i;
   long dictElems;
   bool makeBranches;
   long *binaryDict[B];
   uint32_t curPattern=0;

   assert(L-B > 2);
   assert(B < 31); //B+1 < 32 to make the bit shifts well defined

   //Initialize the data structure tables
   for(j=0; j< B; j++) {
      //For a length m prefix, we need 2^m sets of length 2 arrays.
      //Here, j+1 is the length of the prefix, so we need 2^(j+1) prefixes, or 2*2^(j+1) = 2^(j+2) storage total.
      //Note: 2^(j+2) = 1<<(j+2).
      binaryDict[j] = new long[1U<<(j+2)];

      memset(binaryDict[j], 0, sizeof(long)*(1U<<(j+2)));
   }

   //But the first B bits into curPattern
   curPattern = compressedBitSymbols(S, B);

   maxRunOfCorrects = 0;
   curRunOfCorrects = 0;
   correctCount = 0;
   dictElems=0;

   //In C, arrays are 0 indexed.
   //i is the index of the bit to be predicted.
   //This is all rather confusing, but it helps to run the first and last few cycles, to verify that it works...
   for(i=B+1; i<L; i++) {
      //j is the length of the word to be used
      //3a
      for(j=B; j>0; j--) {
         //update the state to reflect last round's new value (add S[i-1] to the predictor) if our dictionary has room
         //We need the j-tuple prior to S[i-1], that is (S[i-j-1], ..., S[i-2])
         //Now update the state to reflect the new value.
         makeBranches = dictElems < MAX_DICTIONARY_SIZE;
         //This tuple is stored in curPattern. Take the lower j bits.
         if(incrementBinaryDict(binaryDict, j, curPattern, S[i-1], makeBranches, false) && makeBranches) {
            dictElems++;
         }
      }

      maxCount = 0;
      prediction = 0;

      //Add S[i-1] to the curPattern
      curPattern = ((curPattern << 1) | S[i-1]) & ((1U << B) - 1);

      //3b.
      for(j=B; j>0; j--) {
         //Get the predictions
         //predict S[i] by using the prior j bits and the current state
         //We need the j-tuple prior to S[i], that is (S[i-j], ..., S[i-1])
         curPredictionCount = predictBinaryDict(binaryDict, j, curPattern, &curPrediction);

         if(curPredictionCount > maxCount) {
            prediction = curPrediction;
            maxCount = curPredictionCount;
         }
      }

      if((maxCount!=0) && (S[i] == prediction)) {
         correctCount++;
         curRunOfCorrects++;
      } else {
         curRunOfCorrects = 0;
      }

      if(curRunOfCorrects > maxRunOfCorrects) {
         maxRunOfCorrects = curRunOfCorrects;
      }
   }

   for(j=0; j<B; j++) {
      delete[](binaryDict[j]);
      binaryDict[j] = NULL;
   }

   return(predictionEstimate(correctCount, L-B-1, maxRunOfCorrects, 2, "LZ78Y", verbose));
}

// Section 6.3.10 - LZ78Y Prediction Estimate
double LZ78Y_test(byte *data, long len, int alph_size, const bool verbose){
	int dict_size;
	long i, j, N, C, count, max_count, run_len, max_run_len;
	byte y, prediction;
	array<byte, B> prev;
	bool found_prev, have_prediction;


	if(alph_size==2) return binaryLZ78YPredictionEstimate(data, len, verbose);

	// j             prev          y          D[x,y]
	array<map<array<byte, B>, map<byte, long>>, B> D;

	if(len < B+2){	
		printf("\t*** Warning: not enough samples to run LZ78Y test (need more than %d) ***\n", B+2);
		return -1.0;
	}

	N = len-B-1;
	C = 0;
	run_len = 0;
	max_run_len = 0;

	// initialize dictionary counts
	dict_size = 0;
	memset(prev.data(), 0, B);
	for(j = 0; j < B; j++){
		memcpy(prev.data(), data+B-j-1, j+1);
		D[j][prev][data[B]] = 1;
		dict_size++;
	}

	// perform predictions
	for(i = B+1; i < len; i++){
		max_count = 0;
		have_prediction = false;
		memset(prev.data(), 0, B);
		for(j = 0; j < B; j++){
			// check if prev has been previously seen. If (prev,y) has not occurred, 
			// then do not make a prediction for current d and larger d's 
			// as well, since it will not occur for them either. In other words,
			// prediction is NULL, so do not update the scoreboard.
			if((j == 0) || found_prev){
				memcpy(prev.data(), data+i-j-1, j+1);
				if(D[j].find(prev) == D[j].end()) found_prev = false;
				else found_prev = true;
			}

			if(found_prev){
				y = max_map(D[j][prev]);
				count = D[j][prev][y];
				if(count >= max_count){
					max_count = count;
					prediction = y;
					have_prediction = true;
				}
				D[j][prev][data[i]]++;
			}
			else if(dict_size < MAX_DICTIONARY_SIZE){
				memcpy(prev.data(), data+i-j-1, j+1);
				D[j][prev][data[i]] = 1;
				dict_size++;
			}
		}
	
		// test	prediction of maximum (prev,y) pair
		if(have_prediction && (prediction == data[i])){
			C++;
			if(++run_len > max_run_len) max_run_len = run_len;
		}
		else run_len = 0;
	}

	return(predictionEstimate(C, N, max_run_len, alph_size, "LZ78Y", verbose));
}
