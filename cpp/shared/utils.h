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

#define RELEPSILON 4.0*DBL_EPSILON
#define ABSEPSILON 4.0*DBL_EPSILON
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
	long len; 		// number of words in data
	long blen; 		// number of bits in data
};

using namespace std;

//This generally performs a check for relative closeness, but (if that check would be nonsense)
//it can check for an absolute separation, using either the distance between the numbers, or
//the number of ULPs that separate the two numbers.
//See the following for details and discussion of this approach:
//See https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
//https://floating-point-gui.de/errors/comparison/
bool relEpsilonEqual(double A, double B, double maxAbsFactor, double maxRelFactor, uint32_t maxULP)
{
   double diff;
   double absA, absB;
   uint64_t Aint;
   uint64_t Bint;

   assert(sizeof(uint64_t) == sizeof(double));

   ///NaN is by definition not equal to anything (including itself)
   if(std::isnan(A) || std::isnan(B)) {
      return false;
   }

   //We now know that absB>absA
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

      tmp = B;
      B=A;
      A=tmp;

      tmp = absA;
      absA = absB;
      absB = tmp;
   }

   diff=fabs(B-A);

   //Is absA or diff subnormal?
   if((absA <= DBL_MIN) || (diff < DBL_MIN)) {
      //Yes. relitive closeness is going to be nonsense
      return diff < maxAbsFactor;
   } else if(diff <= absB * maxRelFactor) {
      return true;
   }
   //They aren't close in the normal sense, but perhaps that's just due to IEEE representation.
   //Check to see if the value is within maxULP ULPs

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
   //Can't meaningfully compare non-zero values with 0.0 in this way.
   if(A==0.0) {
      return false;
   }
#pragma GCC diagnostic pop

   //if they aren't the same sign, they can't be equal
   if(signbit(A) !=  signbit(B)) {
      return false;
   }

   memcpy(&Aint, &absA, sizeof(double));
   memcpy(&Bint, &absB, sizeof(double));
   assert(Bint > Aint);

   return (Bint - Aint <= maxULP);
}


void free_data(data_t *dp){
	if(dp->symbols != NULL) free(dp->symbols);
	if(dp->rawsymbols != NULL) free(dp->rawsymbols);
	if((dp->word_size > 1) && (dp->bsymbols != NULL)) free(dp->bsymbols);
} 

// Read in binary file to test
bool read_file(const char *file_path, data_t *dp){
	FILE *file; 
	int mask, j, max_symbols = 1 << dp->word_size;
	int symbol_map_down_table[max_symbols];
	long rc, i;

	file = fopen(file_path, "rb");
	if(!file){
		printf("Error: could not open '%s'\n", file_path);
		return false;
	}

	rc = (long)fseek(file, 0, SEEK_END);
	if(rc < 0){
		printf("Error: fseek failed\n");
		fclose(file);
		return false;
	}

	dp->len = ftell(file);
	if(dp->len < 0){
		printf("Error: ftell failed\n");
		fclose(file);
		return false;
	}

	rewind(file);

	if(dp->len == 0){
		printf("Error: '%s' is empty\n", file_path);
		fclose(file);
		return false;
	}

	dp->symbols = (byte*)malloc(sizeof(byte)*dp->len);
	dp->rawsymbols = (byte*)malloc(sizeof(byte)*dp->len);
	if(dp->symbols == NULL){
		printf("Error: failure to initialize memory for symbols\n");
		fclose(file);
		return false;
	}

	rc = fread(dp->symbols, sizeof(byte), dp->len, file);
	if(rc != dp->len){
		printf("Error: file read failure\n");
		fclose(file);
		free(dp->symbols);
		dp->symbols = NULL;
		return false;
	}
	fclose(file);

	memcpy(dp->rawsymbols, dp->symbols, sizeof(byte)* dp->len);
	dp->maxsymbol = 0;

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
		if(symbol_map_down_table[i] != 0) symbol_map_down_table[i] = (byte)dp->alph_size++;
	}

	// create bsymbols (bitstring) using the non-mapped data
	dp->blen = dp->len * dp->word_size;
	if(dp->word_size == 1) dp->bsymbols = dp->symbols;
	else{
		dp->bsymbols = (byte*)malloc(dp->blen);
		if(dp->bsymbols == NULL){
			printf("Error: failure to initialize memory for bsymbols\n");
			free(dp->symbols);
			dp->symbols = NULL;
			return false;
		}

		for(i = 0; i < dp->len; i++){
			for(j = 0; j < dp->word_size; j++){
				dp->bsymbols[i*dp->word_size+j] = (dp->symbols[i] >> (dp->word_size-1-j)) & 0x1;
			}
		}
	}

	// map down symbols if less than 2^bits_per_word unique symbols
	if(dp->alph_size < max_symbols){
		for(i = 0; i < dp->len; i++) dp->symbols[i] = (byte)symbol_map_down_table[dp->symbols[i]];
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
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */
void xoshiro_jump(unsigned int jump_count, uint64_t *xoshiro256starstarState) {
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;

	for(int j=0; j < jump_count; j++) {
		for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
			for(int b = 0; b < 64; b++) {
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

// Return the key that leads to the maximum value
byte max_map(const map<byte, long>& m){
	long max_cnt = 0;
	byte max_key;
	map<byte, long>::const_iterator itr;

	for(itr = m.begin(); itr != m.end(); ++itr){
		if((itr->second > max_cnt) || ((itr->second == max_cnt) && (itr->first > max_key))){
			max_key = itr->first;
			max_cnt = itr->second;
		}
	}

	return max_key;
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

double calc_p_global(long C, long N){
	double p = C/(double)N;

	if(p > 0) p = min(1.0, p + ZALPHA*sqrt((p*(1.0-p))/(N-1.0)));
	else p = 1 - pow(0.01, 1.0/(double)N);
	return p;
}

double calc_p_local(long max_run_len, long N){
	int i, j;
	double p, q, r, log_alpha;
	long double x;
	double lastP, pVal;
	double lvalue, hvalue;
	double hbound, lbound;
	double hdomain, ldomain;

	// binary search for p_local
	r = (double)max_run_len+1;
	log_alpha = log(0.99);
	
	ldomain = 0.0;
	hdomain = 1.0;

	lbound = ldomain;
	hbound = hdomain;

	lvalue = DBL_INFINITY;
	hvalue = -DBL_INFINITY;

	//Note that the bounds are in [0,1], so overflows aren't an issue
	//But underflows are.
	p = (lbound + hbound) / 2.0;

	q = 1.0-p;
	x = 1.0;
	for(i = 0; i < 10; i++) x = 1.0 + q*powl(p, r)*powl(x, r+1.0);
	pVal = (double)(logl(1.0-p*x) - logl((r+1.0-r*x)*q) - (N+1.0)*logl(x));

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

		q = 1.0-p;
		x = 1.0;
		for(i = 0; i < 10; i++) x = 1.0 + q*powl(p, r)*powl(x, r+1.0);
		pVal = log(1.0-p*x) - log((r+1.0-r*x)*q) - (N+1.0)*log(x);

		//invariant: If this isn't true, then this isn't loosly monotonic
		if(!INCLOSEDINTERVAL(pVal, lvalue, hvalue)) {
			p = hbound;
			break;
		}
	}//for loop

	return p;
}
