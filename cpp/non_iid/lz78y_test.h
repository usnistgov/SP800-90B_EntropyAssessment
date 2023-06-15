#pragma once
#include "../shared/utils.h"

#define B_len 16
#define MAX_DICTIONARY_SIZE 65536

static double binaryLZ78YPredictionEstimate(const uint8_t *S, long L, const int verbose, const char *label)
{
   long *binaryDict[B_len];
   long curRunOfCorrects=0;
   long maxRunOfCorrects=0;
   long correctCount=0;
   long i, j;
   uint32_t curPattern=0;
   long dictElems=0;

   assert(L>B_len);
   assert(L-B_len > 2);
   assert(B_len < 32); //B < 32 to make the bit shifts well defined

   //Initialize the data structure tables
   for(j=0; j< B_len; j++) {
      //For a length m prefix, we need 2^m sets of length 2 arrays.
      //Here, j+1 is the length of the prefix, so we need 2^(j+1) prefixes, or 2*2^(j+1) = 2^(j+2) storage total.
      //Note: 2^(j+2) = 1<<(j+2).
      binaryDict[j] = new long[1U<<(j+2)];

      memset(binaryDict[j], 0, sizeof(long)*(1U<<(j+2)));
   }

   // initialize B counts with {(S[15]), S[16]}, {(S[14], S[15]), S[16]}, ..., {(S[0]), S[1], ..., S[15]), S[16]},
   for(j=0; j<B_len; j++) {
      curPattern = curPattern | (((uint32_t)(S[B_len - j - 1]&1)) << j);

      //This is necessarily the first symbol of this length
      (BINARYDICTLOC(j+1, curPattern))[S[B_len]&0x1] = 1;
      dictElems++;
   }

   //In C, arrays are 0 indexed.
   //i is the index of the bit to be predicted.
   for(i=B_len+1; i<L; i++) {
      bool found_x;
      bool havePrediction = false;
      uint8_t roundPrediction=2;
      uint8_t curPrediction=2;
      long maxCount = 0;

      //But the first B bits into curPattern
      curPattern = compressedBitSymbols(S+i-B_len, B_len);

      //j is the length of the prefix to be used
      for(j=B_len; j>0; j--) {
         long curCount;
         long *binaryDictEntry;

         //curPattern starts off as long as possible. We then clear bits at the end
         //as we shorten curPattern
         curPattern = curPattern & ((1U<<j)-1);
         //curPattern should contain the j-tuple (S[i-j] ... S[i-1])

         binaryDictEntry = BINARYDICTLOC(j, curPattern);

          //check if x has been previously seen.
         //For the prediction, roundPrediction is the max across all pairs (there are only 2 symbols here!)
         if((binaryDictEntry[0] > binaryDictEntry[1])) {
            roundPrediction = 0;
            curCount = binaryDictEntry[0];
         } else {
            roundPrediction = 1;
            curCount = binaryDictEntry[1];
         }

         if(curCount == 0) {
            found_x = false;
         } else {
            found_x = true;
         }

         if(found_x) {
            // x is present in the dictionary as a prefix.
            if(curCount > maxCount) {
               maxCount = curCount;
               havePrediction = true;
               curPrediction = roundPrediction;
            }

            binaryDictEntry[S[i]&1]++;
         } else if(dictElems < MAX_DICTIONARY_SIZE) {
            //We didn't find the x prefix, so (x,y) surely can't have occurred.
            //We're allowed to make a new entry. Do so.
            binaryDictEntry[S[i]&1]=1;
            dictElems++;
         }
      }

      // Check to see if the current prediction is correct.
      if(havePrediction && (curPrediction == S[i])) {
            correctCount++;
            curRunOfCorrects++;
            if(curRunOfCorrects > maxRunOfCorrects) maxRunOfCorrects = curRunOfCorrects;
      } else {
            curRunOfCorrects = 0;
      }
   }

   for(j=0; j<B_len; j++) {
      delete[](binaryDict[j]);
      binaryDict[j] = NULL;
   }

   return(predictionEstimate(correctCount, L-B_len-1, maxRunOfCorrects, 2, "LZ78Y", verbose, label));
}

// Section 6.3.10 - LZ78Y Prediction Estimate
double LZ78Y_test(uint8_t *data, long len, int alph_size, const int verbose, const char *label) {
	int dict_size;
	long i, j, N, C, run_len, max_run_len;
	array<uint8_t, B_len> x;

	if(alph_size==2) return binaryLZ78YPredictionEstimate(data, len, verbose, label);

	array<map<array<uint8_t, B_len>, PostfixDictionary>, B_len> D;

	if(len < B_len+2){	
		printf("\t*** Warning: not enough samples to run LZ78Y test (need more than %d) ***\n", B_len+2);
		return -1.0;
	}

	N = len-B_len-1;
	C = 0;
	run_len = 0;
	max_run_len = 0;

	// initialize dictionary counts
	dict_size = 0;
	memset(x.data(), 0, B_len);
	// initialize LZ78Y counts with {(S[15]), S[16]}, {(S[14], S[15]), S[16]}, ..., {(S[0]), S[1], ..., S[15]), S[16]}
	for(j = 1; j <= B_len; j++){
		memcpy(x.data(), data+B_len-j, j);
		D[j-1][x].incrementPostfix(data[B_len], true);
		dict_size++;
	}

	// perform predictions
	for(i = B_len+1; i < len; i++) {
		bool found_x;
		bool have_prediction = false;
		uint8_t prediction = 0;
		long max_count = 0;

		for(j = B_len; j > 0; j--) {
			map<array<uint8_t, B_len>, PostfixDictionary>::iterator curp;

			// check if x has been previously seen. 
			//For the prediction, roundPrediction is the max across all pairs
			//The prefix string should contain the j-tuple (S[i-j] ... S[i-1])
			memset(x.data(), 0, B_len);
			memcpy(x.data(), data+i-j, j);
			curp = D[j-1].find(x);

			if(curp == D[j-1].end()) found_x = false;
			else found_x = true;

			if(found_x) {
				long count;
				uint8_t y;

				// x has occurred, find max (x,y) pair across all y's
				// Check to see if the current prediction is correct.
				y = (curp->second).predict(count);

				if(count > max_count){
					max_count = count;
					prediction = y;
					have_prediction = true;
				}
				//x exists as a prefix, so we always increment (and perhaps add a new postfix)
				(curp->second).incrementPostfix(data[i], true);
			} else if(dict_size < MAX_DICTIONARY_SIZE) {
				//We didn't find the x prefix, so (x,y) surely can't have occurred.
                                //We're allowed to make a new entry. Do so.
                                //curp isn't populated here, because it wasn't found
				D[j-1][x].incrementPostfix(data[i], true);
				dict_size++;
			}
		}
	
		// test	prediction of maximum (x,y) pair
		if(have_prediction && (prediction == data[i])){
			C++;
			if(++run_len > max_run_len) max_run_len = run_len;
		}
		else run_len = 0;
	}

	return(predictionEstimate(C, N, max_run_len, alph_size, "LZ78Y", verbose, label));
}
