//Version of the tool
#define VERSION "1.1.7"

#pragma once
#include <iostream>		// std::cout
#include <string>		// std::string
#include <map>			// std::map
#include <set>			// std::set
#include <string.h>		// strlen
#include <iomanip>		// setw / setfill
#include <stdio.h>
//#include <stdlib.h>
#include <cstdlib>
#include <vector>		// std::vector
#include <time.h>		// time
#include <algorithm>	// std::sort
#include <cmath>		// pow, log2
#include <array>		// std::array
#include <omp.h>		// openmp 4.0 with gcc 4.9
#include <bitset>
#include <mutex>		// std::mutex
#include <assert.h>
#include <cfloat>
#include <math.h>
#include "test_run_base.h"

#define SWAP(x, y) do { int s = x; x = y; y = s; } while(0)
#define INOPENINTERVAL(x, a, b) (((a)>(b))?(((x)>(b))&&((x)<(a))):(((x)>(a))&&((x)<(b))))
#define INCLOSEDINTERVAL(x, a, b) (((a)>(b))?(((x)>=(b))&&((x)<=(a))):(((x)>=(a))&&((x)<=(b))))

#define MIN_SIZE 1000000
#define PERMS 10000

//This is the smallest practical value (one can't do better with the double type)
#define RELEPSILON DBL_EPSILON
//This is clearly overkill, but it's difficult to do better without a view into the monotonic function
#define ABSEPSILON DBL_MIN
#define DBL_INFINITY __builtin_inf ()
#define ITERMAX 1076
#define ZALPHA 2.5758293035489008
#define ZALPHA_L 2.575829303548900384158L

//Make uint128_t a supported type (standard as of C23)
#ifdef __SIZEOF_INT128__
typedef unsigned __int128 uint128_t;
typedef unsigned __int128 uint_least128_t;
# define UINT128_MAX         ((uint128_t)-1)
# define UINT128_WIDTH       128
# define UINT_LEAST128_WIDTH 128
# define UINT_LEAST128_MAX   UINT128_MAX
# define UINT128_C(N)        ((uint_least128_t)+N ## WBU)
#endif

typedef struct data_t data_t;

struct data_t{
	int word_size; 		// bits per symbol
	int alph_size; 		// symbol alphabet size
	uint8_t maxsymbol; 	// the largest symbol present in the raw data stream
	uint8_t *rawsymbols; 	// raw data words
	uint8_t *symbols; 		// data words
	uint8_t *bsymbols; 	// data words as binary string
	long len; 		// number of words in data
	long blen; 		// number of bits in data
};



using namespace std;

//This generally performs a check for relative closeness, but (if that check would be nonsense)
//it can check for an absolute separation, using either the distance between the numbers, or
//the number of ULPs that separate the two numbers.
//See the following for details and discussion of this approach:
//https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
//https://floating-point-gui.de/errors/comparison/
//https://www.boost.org/doc/libs/1_62_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point/floating_points_comparison_theory.html
//Knuth AoCP vol II (section 4.2.2)
//Tested using modified test cases from https://floating-point-gui.de/errors/NearlyEqualsTest.java
bool relEpsilonEqual(double A, double B, double maxAbsFactor, double maxRelFactor, uint32_t maxULP)
{
   double diff;
   double absA, absB;
   uint64_t Aint;
   uint64_t Bint;

   assert(sizeof(uint64_t) == sizeof(double));
   assert(maxAbsFactor >= 0.0);
   assert(maxRelFactor >= 0.0);

   ///NaN is by definition not equal to anything (including itself)
   if(std::isnan(A) || std::isnan(B)) {
      return false;
   }

   //Deals with equal infinities, and the corner case where they are actually copies
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
   if(A==B) {
      return true;
   }
#pragma GCC diagnostic pop

   //If either is infinity, but they are not equal, then they aren't close.
   if(std::isinf(A) || std::isinf(B)) {
      return false;
   }

   absA = fabs(A);
   absB = fabs(B);
   //Make sure that A is the closest to 0.
   if(absA > absB) {
      double tmp;

      //Swap A and B
      tmp = B;
      B=A;
      A=tmp;

      //Swap absA and absB
      tmp = absB;
      absB = absA;
      absA = tmp;
   }

   //Capture the difference of the largest magnitude from the smallest magnitude
   diff=fabs(B-A);

   //Is absA, diff, or absB * maxRelFactor subnormal?
   //Did diff overflow?
   //if absA is subnormal (effectively 0) or 0, then relative difference isn't meaningful, as fabs(B-A)/B≈1 for all values of B
   //In the instance of overflows, the resulting relative comparison will be nonsense.
   if((absA < DBL_MIN) || (diff < DBL_MIN) || std::isinf(diff) || (absB * maxRelFactor < DBL_MIN)) {
      //Yes. Relative closeness is going to be nonsense
      return diff <= maxAbsFactor;
   } else {
      //No. Using relative closeness is probably the right thing to do.
      //Proceeding roughly as per Knuth AoCP vol II (section 4.2.2)
      if(diff <= absB * maxRelFactor) {
         //These are relatively close
         return true;
      } 
   }

   //Neither A or B is subnormal, and they aren't close in the conventional sense, 
   //but perhaps that's just due to IEEE representation. Check to see if the value is within maxULP ULPs.

   //We can't meaningfully compare non-zero values with 0.0 in this way,
   //but absA >= DBL_MIN if we're here, so neither value is 0.0.

   //if they aren't the same sign, then these can't be only a few ULPs away from each other
   if(signbit(A) != signbit(B)) {
      return false;
   }

   //Note, casting from one type to another is undefined behavior, but memcpy will necessarily work
   memcpy(&Aint, &absA, sizeof(double));
   memcpy(&Bint, &absB, sizeof(double));
   //This should be true by the construction of IEEE doubles
   assert(Bint > Aint);

   return (Bint - Aint <= maxULP);
}


void free_data(data_t *dp){
	if(dp->symbols != NULL) free(dp->symbols);
	if(dp->rawsymbols != NULL) free(dp->rawsymbols);
	if((dp->word_size > 1) && (dp->bsymbols != NULL)) free(dp->bsymbols);
} 


// Read in binary file to test
bool read_file_subset(const char *file_path, data_t *dp, unsigned long subsetIndex, unsigned long subsetSize, TestRunBase *testRun) {

	FILE *file; 
	int mask, j, max_symbols;
	long rc, i;
	long fileLen;

	file = fopen(file_path, "rb");
	if(!file){
                testRun->errorLevel = -1;
                testRun->errorMsg = "Error: could not open '%s'\n", file_path;
		printf("Error: could not open '%s'\n", file_path);
		return false;
	}

	rc = (long)fseek(file, 0, SEEK_END);
	if(rc < 0) {
    testRun->errorLevel = -1;
    testRun->errorMsg = "Error: fseek failed";
    printf("Error: fseek failed\n");
		fclose(file);
		return false;
	}

	fileLen = ftell(file);
	if(fileLen < 0){
    testRun->errorLevel = -1;
    testRun->errorMsg = "Error: ftell failed";
    printf("Error: ftell failed\n");
		fclose(file);
		return false;
	}

	rewind(file);

	if(subsetSize == 0) {
		dp->len = fileLen;
	} else {
		rc = (long)fseek(file, subsetIndex*subsetSize, SEEK_SET);
		if(rc < 0){
                        testRun->errorLevel = -1;
                        testRun->errorMsg = "Error: fseek failed";
			printf("Error: fseek failed\n");
			fclose(file);
			return false;
		}

		dp->len = min(fileLen - subsetIndex*subsetSize, subsetSize);
	}

	if(dp->len == 0){
    testRun->errorLevel = -1;
    testRun->errorMsg = "Error: '%s' is empty\n", file_path;
    printf("Error: '%s' is empty\n", file_path);
		fclose(file);
		return false;
	}

	dp->symbols = (uint8_t*)malloc(sizeof(uint8_t)*dp->len);
	dp->rawsymbols = (uint8_t*)malloc(sizeof(uint8_t)*dp->len);
	if((dp->symbols == NULL) || (dp->rawsymbols == NULL)){
    testRun->errorLevel = -1;
    testRun->errorMsg = "Error: failure to initialize memory for symbols";
    printf("Error: failure to initialize memory for symbols\n");
		fclose(file);
		if(dp->symbols != NULL) {
			free(dp->symbols);
			dp->symbols = NULL;
		}
		if(dp->rawsymbols != NULL) {
			free(dp->rawsymbols);
			dp->rawsymbols = NULL;
		}
		return false;
	}

	rc = fread(dp->symbols, sizeof(uint8_t), dp->len, file);
	if(rc != dp->len){
    testRun->errorLevel = -1;
    testRun->errorMsg = "Error: file read failure";
    printf("Error: file read failure\n");
		fclose(file);
		free(dp->symbols);
		dp->symbols = NULL;
		free(dp->rawsymbols);
		dp->rawsymbols = NULL;
		return false;
	}
	fclose(file);

	//Do we need to establish the word size?
	if(dp->word_size == 0) {
		uint8_t datamask = 0;
		uint8_t curbit = 0x80;

		for(i = 0; i < dp->len; i++) {
			datamask = datamask | dp->symbols[i];
		}

		for(i=8; (i>0) && ((datamask & curbit) == 0); i--) {
			curbit = curbit >> 1;
		}

		dp->word_size = i;
	} else {
		uint8_t datamask = 0;
		uint8_t curbit = 0x80;

		for(i = 0; i < dp->len; i++) {
			datamask = datamask | dp->symbols[i];
		}

		for(i=8; (i>0) && ((datamask & curbit) == 0); i--) {
			curbit = curbit >> 1;
		}

		if( i < dp->word_size ) {
			printf("Warning: Symbols appear to be narrower than described.\n");
                        testRun->errorMsg = "Warning: Symbols appear to be narrower than described.";
		} else if( i > dp->word_size ) {
                        testRun->errorLevel = -1;
                        testRun->errorMsg = "Error: Incorrect bit width specification: Data (" + std::to_string(i) + ") does not fit within described bit width: " + std::to_string(dp->word_size) + ".";
			printf("Incorrect bit width specification: Data (%ld) does not fit within described bit width: %d.\n",i,dp->word_size); 
                        free(dp->symbols);
			dp->symbols = NULL;
			free(dp->rawsymbols);
			dp->rawsymbols = NULL;
			return false;
		}
	}

	memcpy(dp->rawsymbols, dp->symbols, sizeof(uint8_t)* dp->len);
	dp->maxsymbol = 0;

	max_symbols = 1 << dp->word_size;
	int symbol_map_down_table[max_symbols];

	// create symbols (samples) and check if they need to be mapped down
	dp->alph_size = 0;
	memset(symbol_map_down_table, 0, max_symbols*sizeof(int));
	mask = max_symbols-1;
	for(i = 0; i < dp->len; i++){ 
		dp->symbols[i] &= mask;
		if(dp->symbols[i] > dp->maxsymbol) dp->maxsymbol = dp->symbols[i];
		if(symbol_map_down_table[dp->symbols[i]] == 0) symbol_map_down_table[dp->symbols[i]] = 1;
	}

	for(i = 0; i < max_symbols; i++){
		if(symbol_map_down_table[i] != 0) symbol_map_down_table[i] = (uint8_t)dp->alph_size++;
	}

	// create bsymbols (bitstring) using the non-mapped data
	dp->blen = dp->len * dp->word_size;
	if(dp->word_size == 1) dp->bsymbols = dp->symbols;
	else{
		dp->bsymbols = (uint8_t*)malloc(dp->blen);
		if(dp->bsymbols == NULL){
			printf("Error: failure to initialize memory for bsymbols\n");
			free(dp->symbols);
			dp->symbols = NULL;
			free(dp->rawsymbols);
			dp->rawsymbols = NULL;

			return false;
		}

		for(i = 0; i < dp->len; i++){
			for(j = 0; j < dp->word_size; j++){
				dp->bsymbols[i*dp->word_size+j] = (dp->symbols[i] >> (dp->word_size-1-j)) & 0x1;
			}
		}
	}

	// map down symbols if less than 2^bits_per_word unique symbols
	if(dp->alph_size < dp->maxsymbol + 1){
		for(i = 0; i < dp->len; i++) dp->symbols[i] = (uint8_t)symbol_map_down_table[dp->symbols[i]];
	} 

	return true;
}

bool read_file(const char *file_path, data_t *dp, TestRunBase *testRun){

	FILE *file; 
	int mask, j, max_symbols;
	long rc, i;

	file = fopen(file_path, "rb");
	if(!file){
                testRun->errorLevel = -1;
                testRun->errorMsg = "Error: could not open '%s'\n", file_path;
		printf("Error: could not open '%s'\n", file_path);
		return false;
	}

	rc = (long)fseek(file, 0, SEEK_END);
	if(rc < 0){
                testRun->errorLevel = -1;
                testRun->errorMsg = "Error: fseek failed";
		printf("Error: fseek failed\n");
		fclose(file);
		return false;
	}

	dp->len = ftell(file);
	if(dp->len < 0){
                testRun->errorLevel = -1;
                testRun->errorMsg = "Error: ftell failed";
		printf("Error: ftell failed\n");
		fclose(file);
		return false;
	}

	rewind(file);

	if(dp->len == 0){
                testRun->errorLevel = -1;
                testRun->errorMsg = "Error: '%s' is empty\n", file_path;
		printf("Error: '%s' is empty\n", file_path);
		fclose(file);
		return false;
	}

	dp->symbols = (uint8_t*)malloc(sizeof(uint8_t)*dp->len);
	dp->rawsymbols = (uint8_t*)malloc(sizeof(uint8_t)*dp->len);
        if((dp->symbols == NULL) || (dp->rawsymbols == NULL)){
                testRun->errorLevel = -1;
                testRun->errorMsg = "Error: failure to initialize memory for symbols";
                printf("Error: failure to initialize memory for symbols\n");
                fclose(file);
                if(dp->symbols != NULL) {
                        free(dp->symbols);
                        dp->symbols = NULL;
                }
                if(dp->rawsymbols != NULL) {
                        free(dp->rawsymbols);
                        dp->rawsymbols = NULL;
                }
                return false;
        }

	rc = fread(dp->symbols, sizeof(uint8_t), dp->len, file);
	if(rc != dp->len){
                testRun->errorLevel = -1;
                testRun->errorMsg = "Error: file read failure";       
		printf("Error: file read failure\n");
		fclose(file);
		free(dp->symbols);
		dp->symbols = NULL;
		free(dp->rawsymbols);
		dp->rawsymbols = NULL;
		return false;
	}
	fclose(file);

	//Do we need to establish the word size?
	if(dp->word_size == 0) {
		//Yes. Establish the word size using the highest order bit in use
		uint8_t datamask = 0;
		uint8_t curbit = 0x80;

		for(i = 0; i < dp->len; i++) {
			datamask = datamask | dp->symbols[i];
		}

		for(i=8; (i>0) && ((datamask & curbit) == 0); i--) {
			curbit = curbit >> 1;
		}

		dp->word_size = i;
       } else {
                uint8_t datamask = 0;
		uint8_t curbit = 0x80;

                for(i = 0; i < dp->len; i++) {
                        datamask = datamask | dp->symbols[i];
                }

                for(i=8; (i>0) && ((datamask & curbit) == 0); i--) {
                        curbit = curbit >> 1;
                }

                if( i < dp->word_size ) {
                        printf("Warning: Symbols appear to be narrower than described.\n");
                        testRun->errorMsg = "Warning: Symbols appear to be narrower than described.";                  
                } else if( i > dp->word_size ) {
                        testRun->errorLevel = -1;
                        testRun->errorMsg = "Error: Incorrect bit width specification: Data (" + std::to_string(i) + ") does not fit within described bit width: " + std::to_string(dp->word_size) + ".";
			printf("Incorrect bit width specification: Data (%ld) does not fit within described bit width: %d.\n",i,dp->word_size);
                        free(dp->symbols);
                        dp->symbols = NULL;
                        free(dp->rawsymbols);
                        dp->rawsymbols = NULL;
                        return false;
                }
        }


	memcpy(dp->rawsymbols, dp->symbols, sizeof(uint8_t)* dp->len);
	dp->maxsymbol = 0;

	max_symbols = 1 << dp->word_size;
	int symbol_map_down_table[max_symbols];

	// create symbols (samples) and check if they need to be mapped down
	dp->alph_size = 0;
	memset(symbol_map_down_table, 0, max_symbols*sizeof(int));
	mask = max_symbols-1;
	for(i = 0; i < dp->len; i++){ 
		dp->symbols[i] &= mask;
		if(dp->symbols[i] > dp->maxsymbol) dp->maxsymbol = dp->symbols[i];
		if(symbol_map_down_table[dp->symbols[i]] == 0) symbol_map_down_table[dp->symbols[i]] = 1;
	}

	for(i = 0; i < max_symbols; i++){
		if(symbol_map_down_table[i] != 0) symbol_map_down_table[i] = (uint8_t)dp->alph_size++;
	}

	// create bsymbols (bitstring) using the non-mapped data
	dp->blen = dp->len * dp->word_size;
	if(dp->word_size == 1) dp->bsymbols = dp->symbols;
	else{
		dp->bsymbols = (uint8_t*)malloc(dp->blen);
		if(dp->bsymbols == NULL){
                        testRun->errorLevel = -1;
                        testRun->errorMsg = "Error: failure to initialize memory for bsymbols";
			printf("Error: failure to initialize memory for bsymbols\n");
                        free(dp->symbols);
			dp->symbols = NULL;
			free(dp->rawsymbols);
			dp->rawsymbols = NULL;
			return false;
		}

		for(i = 0; i < dp->len; i++){
			for(j = 0; j < dp->word_size; j++){
				dp->bsymbols[i*dp->word_size+j] = (dp->symbols[i] >> (dp->word_size-1-j)) & 0x1;
			}
		}
	}

	// map down symbols if less than 2^bits_per_word unique symbols
	if(dp->alph_size < dp->maxsymbol + 1){
		for(i = 0; i < dp->len; i++) dp->symbols[i] = (uint8_t)symbol_map_down_table[dp->symbols[i]];
	} 

	return true;
}

/* This is xoshiro256** 1.0*/
/*This implementation is derived from David Blackman and Sebastiano Vigna, which they placed into
the public domain. See http://xoshiro.di.unimi.it/xoshiro256starstar.c
*/
static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

static inline uint64_t xoshiro256starstar(uint64_t *xoshiro256starstarState)
{
	const uint64_t result_starstar = rotl(xoshiro256starstarState[1] * 5, 7) * 9;
	const uint64_t t = xoshiro256starstarState[1] << 17;

	xoshiro256starstarState[2] ^= xoshiro256starstarState[0];
	xoshiro256starstarState[3] ^= xoshiro256starstarState[1];
	xoshiro256starstarState[1] ^= xoshiro256starstarState[2];
	xoshiro256starstarState[0] ^= xoshiro256starstarState[3];

	xoshiro256starstarState[2] ^= t;

	xoshiro256starstarState[3] = rotl(xoshiro256starstarState[3], 45);

   return result_starstar;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to xoshiro256starstar(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */
void xoshiro_jump(unsigned int jump_count, uint64_t *xoshiro256starstarState) {
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;

	for(unsigned int j=0; j < jump_count; j++) {
		for(unsigned int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
			for(unsigned int b = 0; b < 64; b++) {
				if (JUMP[i] & ((uint64_t)1) << b) {
					s0 ^= xoshiro256starstarState[0];
					s1 ^= xoshiro256starstarState[1];
					s2 ^= xoshiro256starstarState[2];
					s3 ^= xoshiro256starstarState[3];
				}
				xoshiro256starstar(xoshiro256starstarState);	
			}
			
		xoshiro256starstarState[0] = s0;
		xoshiro256starstarState[1] = s1;
		xoshiro256starstarState[2] = s2;
		xoshiro256starstarState[3] = s3;
	}
}

//This seeds using an external source
//We use /dev/urandom here. 
//We could alternately use the RdRand (or some other OS or HW source of pseudo-random numbers)
void seed(uint64_t *xoshiro256starstarState){
	FILE *infp;

	if((infp=fopen("/dev/urandom", "rb"))==NULL) {
		perror("Can't open random source. Reverting to a deterministic seed.");
		exit(-1);
	} 

	if(fread(xoshiro256starstarState, sizeof(uint64_t), 4, infp)!=4) {
		perror("Can't read random seed");
		exit(-1);
	}

	if(fclose(infp)!=0) {
		perror("Couldn't close random source");
		exit(-1);
	}
}

/*Return an integer in the range [0, high], without modular bias*/
/*This is a slight modification of Lemire's approach (as we want [0,s] rather than [0,s)*/
/*See "Fast Random Integer Generation in an Interval" by Lemire (2018) (https://arxiv.org/abs/1805.10941) */
 /* The relevant text explaining the central factor underlying this opaque approach is:
  * "Given an integer x ∈ [0, 2^L), we have that (x × s) ÷ 2^L ∈ [0, s). By multiplying by s, we take
  * integer values in the range [0, 2^L) and map them to multiples of s in [0, s × 2^L). By dividing by 2^L,
  * we map all multiples of s in [0, 2^L) to 0, all multiples of s in [2^L, 2 × 2^L) to one, and so forth. The
  * (i + 1)th interval is [i × 2^L, (i + 1) × 2^L). By Lemma 2.1, there are exactly floor(2^L/s) multiples of s in
  * intervals [i × 2^L + (2^L mod s), (i + 1) × 2^L) since s divides the size of the interval (2^L − (2^L mod s)).
  * Thus if we reject the multiples of s that appear in [i × 2^L, i × 2^L + (2^L mod s)), we get that all
  * intervals have exactly floor(2^L/s) multiples of s."
  *
  * This approach allows us to avoid _any_ modular reductions with high probability, and at worst case one
  * reduction. It's an opaque approach, but lovely.
  */
uint64_t randomRange64(uint64_t s, uint64_t *xoshiro256starstarState){
	uint64_t x;
	uint128_t m;
	uint64_t l;

	x = xoshiro256starstar(xoshiro256starstarState);

	if(UINT64_MAX == s) {
		return x;
	} else {
		s++; // We want an integer in the range [0,s], not [0,s)
		m = (uint128_t)x * (uint128_t)s;
		l = (uint64_t)m; //This is m mod 2^64

		if(l<s) {
			uint64_t t = ((uint64_t)(-s)) % s; //t = (2^64 - s) mod s (by definition of unsigned arithmetic in C)
			while(l < t) {
				x = xoshiro256starstar(xoshiro256starstarState);
				m = (uint128_t)x * (uint128_t)s;
				l = (uint64_t)m; //This is m mod 2^64
			}
		}

		return (uint64_t)(m >> 64U); //return floor(m/2^64)
	}
}

/*
 * This function produces a double that is uniformly distributed in the interval [0, 1).
 * Note that 2^53 is the largest integer that can be represented in a 64 bit IEEE 754 double, such that all 
 * smaller positive integers can also be represented. Shifting the initial random 64-bit value right by 11 
 * bits makes the result only in the lower 53 bits, so the resulting integer is in the range [0, 2^53 - 1].
 * 1.1102230246251565e-16 (0x1.0p-53) is 2^(-53). Multiplying by this value just effects the exponent of the 
 * resulting double, not the significand. We get a double uniformly distributed in the range [0, 1).  
 * The delta between adjacent values is 2^(-53).
 */
double randomUnit(uint64_t *xoshiro256starstarState) {
	return((xoshiro256starstar(xoshiro256starstarState) >> 11) * 1.1102230246251565e-16);
}

// Fisher-Yates Fast (in place) shuffle algorithm
void FYshuffle(uint8_t data[], uint8_t rawdata[], const int sample_size, uint64_t *xoshiro256starstarState) {
	long int r;
	static mutex shuffle_mutex;
	unique_lock<mutex> lock(shuffle_mutex);

	for (long int i = sample_size - 1; i > 0; --i) {
		r = (long int)randomRange64((uint64_t)i, xoshiro256starstarState);
		SWAP(data[r], data[i]);
		SWAP(rawdata[r], rawdata[i]);
	}
}

// Quick sum array  // TODO
long int sum(const uint8_t arr[], const int sample_size) {
	long int sum = 0;
	for (long int i = 0; i < sample_size; ++i) {
		sum += arr[i];
	}

	return sum;
}

// Quick sum std::array // TODO
template<size_t LENGTH>
int sum(const array<int, LENGTH> &arr) {
	int sum = 0;
	for (int i = 0; i < LENGTH; ++i) {
		sum += arr[i];
	}

	return sum;
}

// Quick sum vector
template<typename T>
T sum(const vector<T> &v) {
	T sum = 0;
	for (unsigned long int i = 0; i < v.size(); ++i) {
		sum += v[i];
	}

	return sum;
}

// Calculate baseline statistics
// Finds mean, median, and whether or not the data is binary
void calc_stats(const data_t *dp, double &rawmean, double &median) {

	// Calculate mean
	rawmean = sum(dp->rawsymbols, dp->len) / (double)dp->len;

	// Sort in a vector for median/min/max
	vector<uint8_t> v(dp->symbols, dp->symbols + dp->len);
	sort(v.begin(), v.end());

	long int half = dp->len / 2;
	if(dp->alph_size == 2) {
		//This isn't necessarily true, but we are supposed to set it this way.
		//See 5.1.5, 5.1.6.
		median = 0.5;
	} else {
		if((dp->len & 1) == 1) {
			//the length is odd
			median = v[half];
		} else {
			//the length is even
			median = (v[half] + v[half - 1]) / 2.0;
		}
	}
}


// Map initialization for integers
void map_init(map<uint8_t, int> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0;
	}
}

// Map initialization for doubles
void map_init(map<uint8_t, double> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0.0;
	}
}

// Map initialization for pair<uint8_t, uint8_t> to int
void map_init(map<pair<uint8_t, uint8_t>, int> &m) {
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			m[pair<uint8_t, uint8_t>(i, j)] = 0;
		}
	}
}

// Calculates proportions of each value as an index
void calc_proportions(const uint8_t data[], vector<double> &p, const int sample_size) {
	for (int i = 0; i < sample_size; i++) {
		p[data[i]] += (1.0 / sample_size);
	}
}

// Calculates proportions of each value as an index
void calc_counts(const uint8_t data[], vector<int> &c, const int sample_size) {
	for (int i = 0; i < sample_size; i++) {
		c[data[i]] ++;
	}
}

// Determines the standard deviation of a dataset
double std_dev(const vector<int> x, const double x_mean) {
	double sum = 0.0;

	for (unsigned int i = 0; i < x.size(); i++) {
		sum += pow(x[i] - x_mean, 2);
	}

	return sqrt(sum / x.size());
}

// Quick formula for n choose 2 (which can be simplified to [n^2 - n] / 2)
long int n_choose_2(const long int n) {
	return ((n*n) - n) / 2;
}

vector<uint8_t> substr(const uint8_t text[], const int pos, const int len, const int sample_size) {
	int substr_len = len;

	if (pos + len > sample_size) {
		substr_len = sample_size - pos;
	}

	vector<uint8_t> substring;

	for (int i = 0; i < substr_len; i++) {
		substring.push_back(text[pos + i]);
	}

	return substring;
}

// Fast substring with no bounds checking
array<uint8_t, 16> fast_substr(const uint8_t text[], const int pos, const int len) {
	array<uint8_t, 16> substring = { 0 };

	for (int i = 0; i < len; i++) {
		substring[i] = text[pos + i];
	}

	return substring;
}

template<typename T>
T max_vector(const vector<T> &vals) {
	T max = vals[0];
	for (unsigned int i = 0; i < vals.size(); i++) {
		if (vals[i] > max) {
			max = vals[i];
		}
	}

	return max;
}

template<typename T>
T max_arr(const T* vals, const unsigned int k){
	T max = vals[0];
	for (unsigned int i = 0; i < k; i++){
		if (vals[i] > max) {
			max = vals[i];
		}
	}

	return max;
}

double divide(const int a, const int b) {
	return ((double)a / (double)b);
}

double prediction_estimate_function(long double p, long r, long N){
	long double q, x, xlast=0.0L;
	assert(p > 0.0L); //In fact, it is >= 1/k
	assert(p < 1.0L);

	q = 1.0L-p;
	x = 1.0L;

	//We know that 
	// * x is in the interval [1,1/p] which is a subset of [1,k]
	// * the sequence is monotonic up (so x >= xlast).
	//As such we don't need much fancyness for looking for "equality"
	for(int i = 0; (i <= 65) && ((x - xlast) > (LDBL_EPSILON*x)); i++) {
		xlast = x;
		x = 1.0L + q*powl(p, r)*powl(x, r+1.0L);
		//We expect this convergence to be monotonic up.
		assert(x >= xlast);
		//We expect x<=1/p
		assert(p*x <= 1.0L);
	}

	return((double)(logl(1.0L-p*x) - logl((r+1.0L-r*x)*q) - (N+1.0L)*logl(x)));
}

double calc_p_local(long max_run_len, long N, double ldomain){
	int j;
	double p, log_alpha;
	double lastP, pVal;
	double lvalue, hvalue;
	double hbound, lbound;
	double hdomain;

	// binary search for p_local
	log_alpha = log(0.99);
	
	hdomain = 1.0;

	lbound = ldomain;
	hbound = hdomain;

	lvalue = DBL_INFINITY;
	hvalue = -DBL_INFINITY;

	//Note that the bounds are in [0,1], so overflows aren't an issue
	//But underflows are.
	p = (lbound + hbound) / 2.0;
	pVal = prediction_estimate_function(p, max_run_len+1, N);

	//We don't need the initial pVal invariant, as our initial bounds are infinite.
	//We don't need the initial bounds, as they are set to the domain bounds
	for(j=0; j<ITERMAX; j++) {
		//Have we reached "equality"?
		if(relEpsilonEqual(pVal, log_alpha, ABSEPSILON, RELEPSILON, 4)) break;

		//Now update based on the found pVal
		if(log_alpha < pVal) {
			lbound = p;
			lvalue = pVal;
		} else {
			hbound = p;
			hvalue = pVal;
		}

		//We now verify that ldomain <= lbound < p < hbound <= hdomain
		//and that target in [ lvalue, hvalue ]
		if(lbound >= hbound) {
			p = fmin(fmax(lbound, hbound),hdomain);
			break;
		}

		//invariant. If this isn't true, then we can't evaluate here.
		if(!(INCLOSEDINTERVAL(lbound, ldomain, hdomain) && INCLOSEDINTERVAL(hbound,  ldomain, hdomain))) {
			p = hdomain;
			break;
		}

		//invariant. If this isn't true, then seeking the value within this interval doesn't make sense.
		if(!INCLOSEDINTERVAL(log_alpha, lvalue, hvalue)) {
			p = hdomain;
			break;
		}

		//Update p
		lastP = p;
		p = (lbound + hbound) / 2.0;

		//invariant. If this isn't true, then further calculation isn't really meaningful.
		if(!INOPENINTERVAL(p,  lbound, hbound)) {
			p = hbound;
			break;
		}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
		//Look for a cycle
		if(lastP == p) {
			p = hbound;
			break;
		}
#pragma GCC diagnostic pop

		pVal = prediction_estimate_function(p, max_run_len+1, N);

		//invariant: If this isn't true, then this isn't loosely monotonic
		if(!INCLOSEDINTERVAL(pVal, lvalue, hvalue)) {
			p = hbound;
			break;
		}
	}//for loop

	return p;
}

double predictionEstimate(long C, long N, long max_run_len, long k, const char *testname, const int verbose, const char *label) {
	double curMax;
	double p_local=-1.0;
	double entEst;
	double p_global;
	double p_globalPrime;

	curMax = 1.0 / ((double)k);

	p_global = (double)C/(double)N;

	if(p_global > 0) p_globalPrime = min(1.0, p_global + ZALPHA*sqrt((p_global*(1.0-p_global))/((double)N-1.0)));
	else p_globalPrime = 1.0 - pow(0.01, 1.0/(double)N);

	curMax = fmax(curMax, p_globalPrime);
	if((curMax < 1.0) && (prediction_estimate_function(curMax, max_run_len+1, N) > log(0.99))) {
		p_local = calc_p_local(max_run_len, N, curMax);
		curMax = fmax(curMax, p_local);
	}

	entEst = -log2(curMax);

	if(verbose == 2) {
		if(p_local > 0.0) printf("%s %s Prediction Estimate: N = %ld, Pglobal' = %.17g (C = %ld) Plocal = %.17g (r = %ld)\n", label, testname, N, p_globalPrime, C, p_local, max_run_len+1);
		else printf("%s %s Prediction Estimate: N = %ld, Pglobal' = %.17g (C = %ld) Plocal can't affect result (r = %ld)\n", label, testname, N, p_globalPrime, C, max_run_len+1);
	} else if(verbose == 3) {
		printf("%s %s Prediction Estimate: C = %ld\n", label, testname, C);
		printf("%s %s Prediction Estimate: r = %ld\n", label, testname, max_run_len + 1);
		printf("%s %s Prediction Estimate: N = %ld\n", label, testname, N);
		printf("%s %s Prediction Estimate: P_global = %.17g\n", label, testname, p_global);
		printf("%s %s Prediction Estimate: P_global' = %.17g\n", label, testname, p_globalPrime);

		if(p_local > 0.0) printf("%s %s Prediction Estimate: P_local = %.17g\n", label, testname, p_local);
		else printf("%s %s Prediction Estimate: P_local can't change the result.\n", label, testname);

		printf("%s %s Prediction Estimate: min entropy = %.17g\n", label, testname, entEst);
	}
	return entEst;
}

//The idea here is that we've given an array of pointers (binaryDict). 
//We are trying to produce the address of the length-2 array associated with the length-d prefix "b".
// array The dth index is d-1, so we first find the start of the address space (binaryDict[(d)-1])
//We take the least significant d bits from "b": this is the expression "(b) & ((1U << (d)) - 1)"
//We then multiply this by 2 (as each pattern is associated with a length-2 array) by left shifting by 1.
#define BINARYDICTLOC(d, b) (binaryDict[(d)-1] + (((b) & ((1U << (d)) - 1))<<1))

static uint32_t compressedBitSymbols(const uint8_t *S, long length)
{
   uint32_t retPattern;
   long j;

   assert(length<=32);

   retPattern = 0;

   for(j=0; j<length; j++) {
      assert(S[j] <= 1);
      retPattern = (retPattern << 1) | S[j];
   }

   return retPattern;
}

static void printVersion(string name) {
    cout << name << " " << VERSION << "\n\n";
    cout << "Disclaimer: ";
    cout << "NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.";
    cout << "\n\n";
    cout << "NIST-developed software is expressly provided \"AS IS.\" NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.";
    cout << "\n\n";
    cout << "You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.";
    cout << "\n\n";
}

static string recreateCommandLine(int argc, char* argv[]) {
    string commandLine = "";
    for(int i = 0; i < argc; ++i) {
        commandLine.append(argv[i]);
        if((i + 1) != argc) {
            commandLine.append(" ");
        }
    }
    return commandLine;
}

class PostfixDictionary {
	map<uint8_t, long> postfixes;
	long curBest;
	uint8_t curPrediction;
public:
	PostfixDictionary() { curBest = 0; curPrediction = 0;}
	uint8_t predict(long &count) {assert(curBest > 0); count = curBest; return curPrediction;}
	bool incrementPostfix(uint8_t in, bool makeNew) {
		map<uint8_t, long>::iterator curp = postfixes.find(in);
		long curCount;
		bool newEntry=false;

		if(curp != postfixes.end()) {
			//The entry is already there. We always increment in this case.
			curCount = ++(curp->second);
		} else if(makeNew) {
			//The entry is not here, but we are allowed to create a new entry
			newEntry = true;
			curCount = postfixes[in] = 1;
		} else {
			//The entry is not here, we are not allowed to create a new entry
			return false;
		}

		//Only instances where curCount is set and an increment was performed get here
		if((curCount > curBest) || ((curCount == curBest) && (in > curPrediction))) { 
			curPrediction = in; 
			curBest = curCount; 
		} 

		return newEntry;
	}
};
