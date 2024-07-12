#pragma once

#include <utility>
#include "../shared/utils.h"
#include <cmath>
#include <cstdint>
#include <assert.h>

/*
* ---------------------------------------------
* 		  HELPER FUNCTIONS / VARIABLES
* ---------------------------------------------
*/


/*
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1995, 2000 by Stephen L. Moshier
   *
   * This software is derived from the Cephes Math Library and is
   * incorporated herein by permission of the author.
   *
   * Copyright (c) 1984, 1987, 1989, 2000 by Stephen L. Moshier
   * All rights reserved.
   * 
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *     * Redistributions of source code must retain the above copyright
   *       notice, this list of conditions and the following disclaimer.
   *     * Redistributions in binary form must reproduce the above copyright
   *       notice, this list of conditions and the following disclaimer in the
   *       documentation and/or other materials provided with the distribution.
   *     * Neither the name of the organization nor the
   *       names of its contributors may be used to endorse or promote products
   *       derived from this software without specific prior written permission.
   * 
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
   * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   */

/*The author allowed for its use under the BSD license*/
//https://raw.githubusercontent.com/deepmind/torch-cephes/master/LICENSE.txt
//https://lists.debian.org/debian-legal/2004/12/msg00295.html

static double MACHEP = 1.11022302462515654042E-16;     // 2**-53
static double MAXLOG = 7.09782712893383996732224E2; // log(MAXNUM)
static double MAXNUM = 1.7976931348623158E308;         // 2**1024*(1-MACHEP)
static double PI     = 3.14159265358979323846;         // pi, duh!

static double big = 4.503599627370496e15;
static double biginv =  2.22044604925031308085e-16;

static int sgngam = 0;

/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */

static double A[] = {
   8.11614167470508450300E-4,
   -5.95061904284301438324E-4,
   7.93650340457716943945E-4,
   -2.77777777730099687205E-3,
   8.33333333333331927722E-2
};
static double B[] = {
   -1.37825152569120859100E3,
   -3.88016315134637840924E4,
   -3.31612992738871184744E5,
   -1.16237097492762307383E6,
   -1.72173700820839662146E6,
   -8.53555664245765465627E5
};
static double C[] = {
   /* 1.00000000000000000000E0, */
   -3.51815701436523470549E2,
   -1.70642106651881159223E4,
   -2.20528590553854454839E5,
   -1.13933444367982507207E6,
   -2.53252307177582951285E6,
   -2.01889141433532773231E6
};

#define MAXLGM 2.556348e305

double cephes_igamc(double a, double x);

static double cephes_polevl(double x, double *coef, int N)
{
   double   ans;
   int      i;
   double   *p;

   p = coef;
   ans = *p++;
   i = N;

   do {
      ans = ans * x  +  *p++;
   } while ( --i );

   return ans;
}

static double cephes_p1evl(double x, double *coef, int N)
{
   double   ans;
   double   *p;
   int      i;

   p = coef;
   ans = x + *p++;
   i = N-1;

   do {
      ans = ans * x  + *p++;
   } while ( --i );

   return ans;
}

/* Logarithm of gamma function */
static double cephes_lgam(double x)
{
   double   p, q, u, w, z;
   int      i;

   sgngam = 1;

   if ( x < -34.0 ) {
      q = -x;
      w = cephes_lgam(q); /* note this modifies sgngam! */
      p = floor(q);

      if ( relEpsilonEqual(p, q, DBL_EPSILON, DBL_EPSILON, 4) ) {
         goto loverf;
      }

      i = (int)p; //Note, p is the output of floor.

      if ( (i & 1) == 0 ) {
         sgngam = -1;
      } else {
         sgngam = 1;
      }

      z = q - p;

      if ( z > 0.5 ) {
         p += 1.0;
         z = p - q;
      }

      z = q * sin( PI * z );

      if (  relEpsilonEqual(z, 0.0, DBL_EPSILON, DBL_EPSILON, 4) ) {
         goto loverf;
      }

      /*      z = log(PI) - log( z ) - w;*/
      z = log(PI) - log( z ) - w;
      return z;
   }

   if ( x < 13.0 ) {
      z = 1.0;
      p = 0.0;
      u = x;

      while ( u >= 3.0 ) {
         p -= 1.0;
         u = x + p;
         z *= u;
      }

      while ( u < 2.0 ) {
         if ( relEpsilonEqual(u, 0.0, DBL_EPSILON, DBL_EPSILON, 4) ) {
            goto loverf;
         }

         z /= u;
         p += 1.0;
         u = x + p;
      }

      if ( z < 0.0 ) {
         sgngam = -1;
         z = -z;
      } else {
         sgngam = 1;
      }

      if ( relEpsilonEqual(u, 2.0, DBL_EPSILON, DBL_EPSILON, 4) ) {
         return( log(z) );
      }

      p -= 2.0;
      x = x + p;
      p = x * cephes_polevl( x, B, 5 ) / cephes_p1evl( x, (double *)C, 6);

      return log(z) + p;
   }

   if ( x > MAXLGM ) {
loverf:
      fprintf(stderr, "lgam: OVERFLOW\n");

      return sgngam * MAXNUM;
   }

   q = ( x - 0.5 ) * log(x) - x + log( sqrt( 2*PI ) );

   if ( x > 1.0e8 ) {
      return q;
   }

   p = 1.0/(x*x);

   if ( x >= 1000.0 )
      q += ((   7.9365079365079365079365e-4 * p
                - 2.7777777777777777777778e-3) *p
            + 0.0833333333333333333333) / x;
   else {
      q += cephes_polevl( p, A, 4 ) / x;
   }

   return q;
}

static double cephes_igam(double a, double x)
{
   double ans, ax, c, r;

   if ( (x <= 0) || ( a <= 0) ) {
      return 0.0;
   }

   if ( (x > 1.0) && (x > a ) ) {
      return 1.e0 - cephes_igamc(a,x);
   }

   /* Compute  x**a * exp(-x) / gamma(a)  */
   ax = a * log(x) - x - cephes_lgam(a);

   if ( ax < -MAXLOG ) {
      fprintf(stderr, "igam: UNDERFLOW\n");
      return 0.0;
   }

   ax = exp(ax);

   /* power series */
   r = a;
   c = 1.0;
   ans = 1.0;

   do {
      r += 1.0;
      c *= x/r;
      ans += c;
   } while ( c/ans > MACHEP );

   return ans * ax/a;
}

double cephes_igamc(double a, double x)
{
   double ans, ax, c, yc, r, t, y, z;
   double pk, pkm1, pkm2, qk, qkm1, qkm2;

   if ( (x <= 0) || ( a <= 0) ) {
      return( 1.0 );
   }

   if ( (x < 1.0) || (x < a) ) {
      return( 1.e0 - cephes_igam(a,x) );
   }

   ax = a * log(x) - x - cephes_lgam(a);

   if ( ax < -MAXLOG ) {
      fprintf(stderr, "igamc: UNDERFLOW\n");
      return 0.0;
   }

   ax = exp(ax);

   /* continued fraction */
   y = 1.0 - a;
   z = x + y + 1.0;
   c = 0.0;
   pkm2 = 1.0;
   qkm2 = x;
   pkm1 = x + 1.0;
   qkm1 = z * x;
   ans = pkm1/qkm1;

   do {
      c += 1.0;
      y += 1.0;
      z += 2.0;
      yc = y * c;
      pk = pkm1 * z  -  pkm2 * yc;
      qk = qkm1 * z  -  qkm2 * yc;

      if ( ! relEpsilonEqual(qk, 0.0, DBL_EPSILON, DBL_EPSILON, 4) ) {
         r = pk/qk;
         t = fabs( (ans - r)/r );
         ans = r;
      } else {
         t = 1.0;
      }

      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      if ( fabs(pk) > big ) {
         pkm2 *= biginv;
         pkm1 *= biginv;
         qkm2 *= biginv;
         qkm1 *= biginv;
      }
   } while ( t > MACHEP );

   return ans*ax;
}


//This document is using Pearson's chi-squared test
//https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test
//The underlying distribution for this test is the chi-squared distribution.
//https://en.wikipedia.org/wiki/Chi-squared_distribution
//The CDF for the chi square distribution with test statistic x and k degrees of freedom is
//P(k/2, x/2), where P is the regularized gamma function.
//For the p-value, we thus need 1-P(k/2, x/2), which is a different regularized gamma function, Q(k/2, x/2)
//https://en.wikipedia.org/wiki/Incomplete_gamma_function#Regularized_Gamma_functions_and_Poisson_random_variables
//Using cephes, we use the igamc:
/* y = igamc( a, x );
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                            inf.
 *                              -
 *                     1       | |  -t  a-1
 *               =   -----     |   e   t   dt.
 *                    -      | |
 *                   | (a)    -
 *                             x
 */
//In Wikipedia terms, this is Gamma(a,x) / Gamma(a) = Q(a,x).
//Thus, the p-value associated with the test statistica T in a Pearson's chi-square test with k degrees of freedom is
// igamc( k/2, x/2 )
double chi_square_pvalue(double x, double k){
	return cephes_igamc(k/2.0, x/2.0);
}

/*
* ---------------------------------------------
* 	  HELPERS FOR CHI_SQUARE_INDEPENDENCE
* ---------------------------------------------
*/
struct tupleTranslateEntry {
        uint16_t  tuple; //a tuple is rendered into a single value
        double expectation;
        int bin;
};

void independence_calc_expectations(const vector<double> &p, vector<struct tupleTranslateEntry> &e, const int sample_size){
	uint16_t index;
	assert(p.size() <= UINT8_MAX + 1);
	for(unsigned long i = 0; i < p.size(); i++){
		for(unsigned long j = 0; j < p.size(); j++){
			index = (uint16_t) ((i*p.size()) + j);
			e[index].tuple = index;
			e[index].expectation = p[i] * p[j] * floor(sample_size * 0.5);
			e[index].bin = -1;
		}
	}
}

void allocate_bins(vector<struct tupleTranslateEntry> &e, vector<double> &bin_exp){
	int current_bin = 0;
	double current_expectation = 0.0;

	for(unsigned int i = 0; i < e.size(); i++){
		if(current_expectation >= 5.0) {
			bin_exp.push_back(current_expectation);
			current_bin ++;
			current_expectation = 0.0;
		}

		e[i].bin = current_bin;
		current_expectation += e[i].expectation;
	}

	//If the current_bin is 0, we can't combine anything
	if((current_bin != 0) && (current_expectation < 5.0)) {
		//Combine the last two bins
		for(unsigned int i =  e.size() - 1; e[i].bin == current_bin; i--) {
			e[i].bin = current_bin - 1;
		}
		bin_exp[current_bin-1] += current_expectation;
	} else {
		bin_exp.push_back(current_expectation);
	}
}

void independence_calc_observed(const uint8_t data[], const vector<struct tupleTranslateEntry> &e, vector<int> &o, const int sample_size, const int alphabet_size){
	for(int j = 0; j < sample_size-1; j+=2){
		uint16_t index = (uint16_t)((data[j] * alphabet_size) + data[j+1]);
		o[e[index].bin]++;
	}
}

double calc_T(const vector<double> &bin_expectations, const vector<int> &o){
	double T = 0.0;

	assert(bin_expectations.size() == o.size());
	
	//fprintf(stderr, "nbins: %zu\n", bin_expectations.size());

	for (unsigned int i = 0; i < bin_expectations.size(); i++){
		//fprintf(stderr, "bin index %u: binCount = %u, binExp = %.17g\n", i, o[i], bin_expectations[i]);
		T += pow((o[i] - bin_expectations[i]), 2) / bin_expectations[i];
	}

	return T;
}

void goodness_of_fit_calc_observed(const uint8_t data[], const vector<struct tupleTranslateEntry> &e, vector<int> &o, const int sample_size){
	for(int j = 0; j < sample_size; j++){
		o[e[data[j]].bin]++;
	}
}

/*
* ---------------------------------------------
* 		  			 TESTS
* ---------------------------------------------
*/

void binary_chi_square_independence(const uint8_t data[], double &score, int &df, const int sample_size){

	// Compute proportion of 0s and 1s
	double p0 = 0.0, p1 = 0.0;
	unsigned int tuple_count;

	for(int i = 0; i < sample_size; i++){
		p1 += data[i];
	}

	p1 /= sample_size;
	p0 = 1.0 - p1;

	// Compute m
	double min_p = min(p0, p1);
	int m = 11;
	int threshhold = 5;
	while(m > 1){
		if (pow(min_p, m) * (sample_size / m) >= threshhold){
			break;
		}else{
			m--;
		}
	}

	//fprintf(stderr, "chi_square m: %u\n", m);
	tuple_count = 1 << m;

	if (m < 2){
		score = 0.0;
		df = 0;
		return;
	}

	// Test is only run if m >= 2
	double T = 0;

	// Count occurances of m-bit tuples by converting to decimal and using as index in a vector
	vector<int> occ(tuple_count, 0);
	int block_count = sample_size / m;

	for(int i = 0; i < block_count; i++){

		int symbol = 0;
		for(int j = 0; j < m; j++){
			symbol = (symbol << 1) | data[i*m + j];
		}

		occ[symbol]++;
	}

	for(unsigned int i = 0; i < occ.size(); i++){

		// GCC only, counts the number of 1s in an integer
		int w = __builtin_popcount(i);

		double e = pow(p1, w) * pow(p0, m - w) * block_count;

		T += pow(occ[i] - e, 2) / e;
	}

	score = T;
	df = pow(2, m) - 2;
}

bool expectationOrder(const struct tupleTranslateEntry &a, const struct tupleTranslateEntry &b){
	if(a.expectation != b.expectation) return a.expectation < b.expectation;
	else return a.tuple < b.tuple;
}

bool tupleOrder(const struct tupleTranslateEntry &a, const struct tupleTranslateEntry &b){
	return a.tuple < b.tuple;
}

void chi_square_independence(const uint8_t data[], double &score, int &df,  const int sample_size, const int alphabet_size){
	// Proportion of each element to the entire set
	vector<double> p(alphabet_size, 0.0);
	calc_proportions(data, p, sample_size);

	// Calculate the expected number of occurrences for each possible pair of symbols
	vector<struct tupleTranslateEntry> e(alphabet_size*alphabet_size);
	independence_calc_expectations(p, e, sample_size);
	//Sort by expectation, from smallest to largest. Secondary sort on tuple value, from smallest to largest
	sort(e.begin(), e.end(), expectationOrder);

	// Allocate sorted expected values into bins and accumulate corresponding expected values of entire bins
	vector<double> bin_expectations;
	allocate_bins(e, bin_expectations);

	//Sort by tuple, from smallest to largest, so we can use this as a lookup table
	sort(e.begin(), e.end(), tupleOrder);
	
	// Calculate the observed frequency of each pair of symbols
	vector<int> o(bin_expectations.size(), 0);
	independence_calc_observed(data, e, o, sample_size, alphabet_size);

	// Calcualte T 
	score = calc_T(bin_expectations, o);

	// Return score and degrees of freedom
	df = bin_expectations.size() - alphabet_size;
}

void binary_goodness_of_fit(const uint8_t data[], double &score, int &df, const int sample_size){

	// Find proportion of 1s to the whole data set
	int sublength = sample_size / 10;
	int ones = 0;

	for(int i = 0; i < sample_size; i++){
		ones += data[i];
	}

	double p = divide(ones, sample_size);
	double T = 0;

	// Compute expected 0s and 1s in each sub-sequence
	double e0 = (1.0 - p) * sublength;
	double e1 = p * sublength;

	for(int i = 0; i < 10; i++){

		// Count actual 0s and 1s in each sub-sequence
		int o0 = 0, o1 = 0;

		for(int j = 0; j < sublength; j++){
			o1 += data[i*sublength + j];
		}

		o0 = sublength - o1;

		// Compute T
		T += (pow(o0 - e0, 2) / e0) + (pow(o1 - e1, 2) / e1);
	}

	score = T;
	df = 9;
}

void goodness_of_fit(const uint8_t data[], double &score, int &df, const int sample_size, const int alphabet_size){
        vector<double> p(alphabet_size, 0.0);
        calc_proportions(data, p, sample_size);

        // Calculate the expected number of occurrences for each possible pair of symbols
        // Calculate the expected number of occurrences for each possible pair of symbols
        vector<struct tupleTranslateEntry> e(alphabet_size);
	for(long int j=0; j < alphabet_size; j++) {
		e[j].tuple = j;
		e[j].expectation = p[j] * floor((double) sample_size / 10.0);
		e[j].bin = -1;
	}

        //Sort by expectation, from smallest to largest. Secondary sort on tuple value, from smallest to largest
        sort(e.begin(), e.end(), expectationOrder);

        // Allocate sorted expected values into bins and accumulate corresponding expected values of entire bins
        vector<double> bin_expectations;
        allocate_bins(e, bin_expectations);

        //Sort by tuple, from smallest to largest, so we can use this as a lookup table
        sort(e.begin(), e.end(), tupleOrder);

	// Calculate the observed frequency of each symbol in each subset
	int block_size = sample_size/10;
	double T = 0.0;
	vector<int> o(bin_expectations.size());

	for(int j=0; j<10; j++) {
		for(unsigned int i=0; i<o.size(); i++) o[i] = 0;
		goodness_of_fit_calc_observed(data+j*block_size, e, o, block_size);
		T += calc_T(bin_expectations, o);
	}

	// Return score and degrees of freedom
	score = T;
	df = 9*(bin_expectations.size()-1);
}

bool chi_square_tests(const uint8_t data[], const int sample_size, const int alphabet_size, const int verbose){

	double score = 0.0;
	double pvalue;
	int df = 0;
	bool result = true;

	// Chi Square independence test
	if(alphabet_size == 2){
		binary_chi_square_independence(data, score, df, sample_size);
	}else{
		chi_square_independence(data, score, df, sample_size, alphabet_size);
	}

	pvalue = chi_square_pvalue(score, df);

	// Print results
	if(verbose == 2) {
		printf("Chi square independence\n");
		printf("\tscore = %f\n", score);
		printf("\tdegrees of freedom = %d\n", df);
		printf("\tp-value = %f\n\n", pvalue);
	} else if(verbose >= 3) {
		printf("Chi square independence: T = %.17g\n", score);
		printf("Chi square independence: df = %d\n", df);
		printf("Chi square independence: P-value = %.17g\n", pvalue);
	}

	// Check result to return if test failed
	if(pvalue < 0.001){
		result = false;
	}

	// Reset score and df
	score = 0.0;
	df = 0;

	// Chi Square goodness of fit test
	if(alphabet_size == 2){
		binary_goodness_of_fit(data, score, df, sample_size);
	}else{
		goodness_of_fit(data, score, df, sample_size, alphabet_size);
	}

	pvalue = chi_square_pvalue(score, df);

	// Print results
	if(verbose == 2) {
		printf("Chi square goodness of fit\n");
		printf("\tscore = %f\n", score);
		printf("\tdegrees of freedom = %d\n", df);
		printf("\tp-value = %f\n\n", pvalue);
	} else if(verbose >= 3) {
		printf("Chi square goodness of fit: T = %.17g\n", score);
		printf("Chi square goodness of fit: df = %d\n", df);
		printf("Chi square goodness of fit: P-value = %.17g\n", pvalue);
	}

	// Check result to return if test failed
	if(pvalue < 0.001){
		result = false;
	}

	return result;
}
