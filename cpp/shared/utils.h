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

typedef unsigned char byte;

typedef struct data_t data_t;

struct data_t{
	int word_size; 		// bits per symbol
	int alph_size; 		// symbol alphabet size
	byte maxsymbol; 	// the largest symbol present in the raw data stream
	byte *rawsymbols; 	// raw data words
	byte *symbols; 		// data words
	byte *bsymbols; 	// data words as binary string
    long len;           // number of symbols in data
	long blen; 		    // number of bits in data
    bool packed;        // whether the data is bit-packed; this is important to retain
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
    /* For packed data, rawsymbols is the same structure as symbols, so don't free twice */
	if(!dp->packed && dp->symbols != NULL) { free(dp->symbols); dp->symbols = NULL; }
	if(dp->rawsymbols != NULL) { free(dp->rawsymbols); dp->rawsymbols = NULL; }
    /* For bit-based data, the binary data is the same as the symbols data, so don't free twice */
	if((dp->word_size > 1) && (dp->bsymbols != NULL)) { free(dp->bsymbols); dp->bsymbols = NULL; }
} 


bool read_file_subset(const char *file_path, data_t *dp, bool packed, unsigned long subsetIndex, unsigned long subsetSize){
	FILE *file = NULL; 
	int mask = 0, j = 0, max_symbols = 0;
	long rc = 0, i = 0;
    byte datamask = 0;
    byte curbit = 0x80;
    int *symbol_map_down_table = NULL;
    int shift_start = 0;
    int tmp_symbol = 0;
    bool status = false;
    byte *inputdata = NULL;
    int inputlen = 0;

	file = fopen(file_path, "rb");
	if(!file){
		printf("Error: could not open '%s'\n", file_path);
        goto error_die;
	}

	rc = (long)fseek(file, 0, SEEK_END);
	if(rc < 0){
		printf("Error: fseek failed\n");
        goto error_die;
	}

    inputlen = ftell(file);
    if(inputlen < 0){
        printf("Error: ftell failed\n");
        goto error_die;
    }

	rewind(file);

    /* Disable bit-packing for bit sizes that have harder semantics to deal with.
     * Technically the loading code below will work for all bit sizes, but for things
     * like ranges and extents and transposition, bit sizes like 3 are difficult to
     * deal with at the end of the range.
     */
    dp->packed = packed;
    if(dp->packed)  {
        /* If dp->word_size is 0, then the parameter wasn't provided which isn't allowed with
         * bit-packing.
         */
        switch(dp->word_size)  {
            case 0: printf("The word size must be provided explicitly when bit-packed samples are used.\n");
                    goto error_die;
            case 3:
            case 5:
            case 6:
            case 7: printf("Bit-packed samples with sample sizes of size %d are not yet supported for all cases.\n", dp->word_size);
                    goto error_die;
        }
    }

    if(subsetSize > 0 && dp->packed)  {
        /* Subsets are essentially a request to locate a substring of samples.
         * The samples are of size word_size; the size of the sample is the number of samples.
         * The index is the ith substring of this length.
         * This means the input has to have at LEAST subsetIndex*subsetSize samples.
         * Indexing starts at 0.
         * A subset size must be > 0 to have any meaning.
         * Units for the index is in "number of substrings of the given size".
         * Units for the size is in "samples".
         * If we want the 3rd substring of 12 samples, where a sample is 3 bits, then
         * we need an input file that has AT LEAST ceil(3*12*3 bits) (14 bytes).
         */
        /* Then we get into a situation where the ith subset could start in the middle
         * of a byte which definitely is an issue when trying to read from the raw
         * data file.  So we don't support these semantics (yet).
         */
        printf("Subset semantics with bit-packed samples is not yet supported.\n");
        goto error_die;
    }

    if(subsetSize > 0) {
        rc = (long)fseek(file, subsetIndex*subsetSize, SEEK_SET);
        if(rc < 0){
            printf("Error: fseek failed\n");
            goto error_die;
        }

        inputlen = min(inputlen - subsetIndex*subsetSize, subsetSize);
    }

	if(inputlen == 0){
		printf("Error: '%s' is empty\n", file_path);
        goto error_die;
	}

	inputdata = (byte*)malloc(sizeof(byte)*inputlen);
    if(!inputdata){
        printf("Error: failure to initialize memory for %d input data\n", inputlen);
        goto error_die;
    }

	rc = fread(inputdata, sizeof(byte), inputlen, file);
	if(rc != inputlen){
		printf("Error: file read failure\n");
        goto error_die;
	}
	fclose(file);
    file = NULL;


    for(i = 0; i < inputlen; i++) {
        datamask = datamask | inputdata[i];
    }

    for(i=8; (i>0) && ((datamask & curbit) == 0); i--) {
        curbit = curbit >> 1;
    }

	//Do we need to establish the word size?
	if(dp->word_size == 0) {
		//Yes. Establish the word size using the highest order bit in use
		dp->word_size = i;
    } else {
        /* Else just do sanity checking and report back to user */
        if( i < dp->word_size ) {
            printf("Warning: Symbols appear to be narrower than described.\n");
        } else { 
            if( !packed && i > dp->word_size ) {
                /* This error only occurs when the data is not packed. */
                printf("Warning: Incorrect bit width specification: Data does not fit within described bit width.\n");
                goto error_die;
            }
            /* When using bit-packed data, we can run into situations where we might not be able to determine what
             * the last sample is supposed to be if it only takes up a fraction of a byte.  For instance, if there
             * are 3-bit samples and the last byte has only one sample in it, then what are the last 5 bits supposed
             * to be interpreted as?  For this reason, we check to make sure that -- at the very least -- the last
             * byte could potentially be composed completely.  We do this based on the input length.  It isn't
             * totally accurate and the worst case scenario occurs when when the input is bit-oriented and there
             * is (for example) only one significant sample bit in the very last byte set (meaning the remaining 7
             * bits have no interpretative value).  A better solution might be to simply trim the input to the 
             * largest set of non-ambiguous interpretation.
             * Instead, we will just warn of the possibility when it is obvious (eg. 3 bit sample in an 8-bit byte).
             */
            if( packed && (inputlen * 8 % dp->word_size) != 0)  {
                printf("Warning: Final byte has ambiguous number of samples. Consider truncating the file to %d bytes.\n", inputlen / dp->word_size * dp->word_size);
            }
        }
    }

    /* After reading the raw symbols, we can now determine how large the symbols will be
     * that we want to process.  If word_size == 8, then it's just a bit-for-bit copy of
     * the raw symbols. Else, we need to mask off the data words OR do bit ops to get the
     * packed data, depending on how the user wanted to process the file.
     * This is the fundamental difference between "rawsymbols" and "symbols".
     * If the data is NOT bit-packed, 'rawsymbols' is before any bit-masking and 'symbols'
     * is after bit-masking.  This difference is not apparent with packed data and
     * 'rawsymbols' will be the same as 'symbols' in this case.  This is only important
     * with IID processing.
     */
    if(packed)  {
        /* Truncation division makes sense here because we can't round up to bits that
         * don't exist.
         */
        dp->len = (inputlen * 8) / dp->word_size;
    } else  {
        /* If we are not doing bit-packing, then we have significant bits in each
         * byte (just depends on how many we care about).  Thus, the number of 
         * samples is the same as the byte size.
         */
        dp->len = inputlen;
    }


    dp->rawsymbols = (byte*)malloc(sizeof(byte)*dp->len);
    if(!dp->rawsymbols)  {
        printf("Error: failure to initialize memory for %d raw symbols\n", dp->len);
        goto error_die;
    }

    /* If not packed, the rawsymbols and symbols have different semantics.
     * If data is packed, then these structures are the same.  No need to duplicate memory.
     */
    if(dp->packed)
        dp->symbols = dp->rawsymbols;
    else  {
        /* And raw symbols will be the input data. */
        memcpy(dp->rawsymbols, inputdata, sizeof(byte)*dp->len);
        dp->symbols = (byte*)malloc(sizeof(byte)*dp->len);
    }

    if(!dp->symbols)  {
        printf("Error: failure to initialize memory for %d symbols\n", dp->len);
        goto error_die;
    }

	dp->maxsymbol = 0;

	max_symbols = 1 << dp->word_size;
	symbol_map_down_table = (int *)malloc(sizeof(int)* max_symbols);
    if(!symbol_map_down_table)  {
        printf("Error! Unable to allocate memory for symbol map down table.\n");
        goto error_die;
    }

	// create symbols (samples) and check if they need to be mapped down
	dp->alph_size = 0;
	memset(symbol_map_down_table, 0, max_symbols*sizeof(int));
	mask = max_symbols-1;
    shift_start = (8-dp->word_size);
    tmp_symbol = 0;

    /* We need to process all of the input data */
	for(i = 0, j = 0; i < inputlen; i++){ 
        /* If we are dealing with packed bits, then we need to do bit-shifting, possibly
         * between two adjacent samples. 
         */
        if(packed)  {
            dp->symbols[j] |= tmp_symbol;   /* Possibly continuing from before */
            do {
                /* Sanity check */
                if (j > dp->len)  {
                    printf("Error: unpacking bits resulted in symbol overflow\n");
                    goto error_die;
                }

                dp->symbols[j] = ((inputdata[i] >> shift_start) & mask);
                shift_start -= dp->word_size;

                /* Need to ensure we get all map down symbols captured in the bit-packed scenario. */
    		    if(dp->symbols[j] > dp->maxsymbol) dp->maxsymbol = dp->symbols[j];
    	    	if(symbol_map_down_table[dp->symbols[j]] == 0) symbol_map_down_table[dp->symbols[j]] = 1;

                j++;
            } while(shift_start >= 0);

            /* At this point, we've exited the bit shifting loop because we ran out of bits in the current
             * sample. This means shift_start < 0.  It is possible that we didn't get *enough* bits to
             * form a viable sample.  So we need to use the residual bits as a temporary symbol and
             * augment with the remaining bits on the next iteration.
             * The way we know if we ran out of bits while in the middle of a sample extraction is if
             * |shift_start| < wordsize.
             */
            if((-1*shift_start) < dp->word_size)  {
                /* Ran out of bits before a full symbol was realized. So concat with next sample on next iteration. */
                tmp_symbol = inputdata[i] & ((1 << (shift_start + dp->word_size)) - 1);  /* Make a temp symbol based on remaining bits */
                shift_start = 8 - (shift_start * -1);  /* new start location of shifting depends on how many residual bits we need */
            }
            else  {
                /* Did not run out of bits when composing. Therefore, tmp_symbol should be zeroed out and shift_start
                 * starts at the beginning again. 
                 */
                tmp_symbol = 0;
                shift_start = 8 - dp->word_size;
            }
        }   /* End of bit-packing scenario */
        else  {
		    dp->symbols[i] = dp->rawsymbols[i] & mask;
		    if(dp->symbols[i] > dp->maxsymbol) dp->maxsymbol = dp->symbols[i];
    		if(symbol_map_down_table[dp->symbols[i]] == 0) symbol_map_down_table[dp->symbols[i]] = 1;
        }
	}

	for(i = 0; i < max_symbols; i++){
		if(symbol_map_down_table[i] != 0) symbol_map_down_table[i] = (byte)dp->alph_size++;
	}

	// create bsymbols (bitstring) using the non-mapped (non-raw) data
	dp->blen = dp->len * dp->word_size;
	if(dp->word_size == 1) dp->bsymbols = dp->symbols;
	else{
		dp->bsymbols = (byte*)malloc(dp->blen);
		if(dp->bsymbols == NULL){
			printf("Error: failure to initialize memory for bsymbols\n");
            goto error_die;
		}

		for(i = 0; i < dp->len; i++){
			for(j = 0; j < dp->word_size; j++){
				dp->bsymbols[i*dp->word_size+j] = (dp->symbols[i] >> (dp->word_size-1-j)) & 0x1;
			}
		}
	}

	// map down symbols if less than 2^bits_per_word unique symbols
	if(dp->alph_size < dp->maxsymbol + 1){
		for(i = 0; i < dp->len; i++) dp->symbols[i] = (byte)symbol_map_down_table[dp->symbols[i]];
	} 

    /* In both success and failure, need to clear the symbol_map_down_table memory */
    free(symbol_map_down_table); 
    symbol_map_down_table = NULL;

    status = true;
    goto exit_function;

error_die:
    free_data(dp);
    status = false;

    /* Common exit point */
exit_function:
    if(symbol_map_down_table) free(symbol_map_down_table); 
    if(file) fclose(file); 
    if(inputdata) free(inputdata);

    return status;
}

/* Function alias */
bool read_file(const char *file_path, data_t *dp, bool packed)  {
    return read_file_subset(file_path, dp, packed, 0, 0);
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

    if (getenv("__SP80090B_TESTING__") == NULL)  {
    	if(fread(xoshiro256starstarState, sizeof(uint64_t), 4, infp)!=4) {
	    	perror("Can't read random seed");
		    exit(-1);
	    }
    } else  {
        /* This static seed obtained by the following: 
         * dd if=/dev/urandom bs=$[64/8] count=4 | hexdump -v -e '16/1 "0x%x," "\n"'
         */
        byte staticTestSeed[sizeof(uint64_t)*4] = {
            0xdd,0xe9,0x96,0x88,0x0,0x39,0xc,0xd0,0xa4,0xf4,0x64,0x30,0xd0,0x8d,0x2c,0x71,
            0xf7,0xf,0x41,0x5a,0x7e,0xd8,0x10,0x3b,0xf7,0x11,0x1,0x54,0x5b,0x87,0x12,0x99
        };
        memcpy(xoshiro256starstarState, staticTestSeed, sizeof(uint64_t) * 4);
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
	__uint128_t m;
	uint64_t l;

	x = xoshiro256starstar(xoshiro256starstarState);

	if(UINT64_MAX == s) {
		return x;
	} else {
		s++; // We want an integer in the range [0,s], not [0,s)
		m = (__uint128_t)x * (__uint128_t)s;
		l = (uint64_t)m; //This is m mod 2^64

		if(l<s) {
			uint64_t t = ((uint64_t)(-s)) % s; //t = (2^64 - s) mod s (by definition of unsigned arithmetic in C)
			while(l < t) {
				x = xoshiro256starstar(xoshiro256starstarState);
				m = (__uint128_t)x * (__uint128_t)s;
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
void FYshuffle(byte data[], byte rawdata[], const int sample_size, uint64_t *xoshiro256starstarState) {
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
long int sum(const byte arr[], const int sample_size) {
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
	vector<byte> v(dp->symbols, dp->symbols + dp->len);
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
void map_init(map<byte, int> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0;
	}
}

// Map initialization for doubles
void map_init(map<byte, double> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0.0;
	}
}

// Map initialization for pair<byte, byte> to int
void map_init(map<pair<byte, byte>, int> &m) {
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			m[pair<byte, byte>(i, j)] = 0;
		}
	}
}

// Calculates proportions of each value as an index
void calc_proportions(const byte data[], vector<double> &p, const int sample_size) {
	for (int i = 0; i < sample_size; i++) {
		p[data[i]] += (1.0 / sample_size);
	}
}

// Calculates proportions of each value as an index
void calc_counts(const byte data[], vector<int> &c, const int sample_size) {
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

vector<byte> substr(const byte text[], const int pos, const int len, const int sample_size) {
	int substr_len = len;

	if (pos + len > sample_size) {
		substr_len = sample_size - pos;
	}

	vector<byte> substring;

	for (int i = 0; i < substr_len; i++) {
		substring.push_back(text[pos + i]);
	}

	return substring;
}

// Fast substring with no bounds checking
array<byte, 16> fast_substr(const byte text[], const int pos, const int len) {
	array<byte, 16> substring = { 0 };

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

	if(verbose == 1) {
		if(p_local > 0.0) printf("%s %s Prediction Estimate: N = %ld, Pglobal' = %.17g (C = %ld) Plocal = %.17g (r = %ld)\n", label, testname, N, p_globalPrime, C, p_local, max_run_len+1);
		else printf("%s %s Prediction Estimate: N = %ld, Pglobal' = %.17g (C = %ld) Plocal can't affect result (r = %ld)\n", label, testname, N, p_globalPrime, C, max_run_len+1);
	} else if(verbose == 2) {
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

static uint32_t compressedBitSymbols(const byte *S, long length)
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

class PostfixDictionary {
	map<byte, long> postfixes;
	long curBest;
	byte curPrediction;
public:
	PostfixDictionary() { curBest = 0; curPrediction = 0;}
	byte predict(long &count) {assert(curBest > 0); count = curBest; return curPrediction;}
	bool incrementPostfix(byte in, bool makeNew) {
		map<byte, long>::iterator curp = postfixes.find(in);
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
