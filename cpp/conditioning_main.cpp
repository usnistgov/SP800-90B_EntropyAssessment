#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <climits>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <getopt.h>
#include <mpfr.h>
#include <errno.h>
#include <fenv.h>

[[ noreturn ]] static void print_usage() {
	printf("Usage is: ea_conditioning [-v] <n_in> <n_out> <nw> <h_in>\n");
	printf("\tor \n\tea_conditioning -n <n_in> <n_out> <nw> <h_in> <h'>\n\n");
	printf("\t <n_in>: input number of bits to conditioning function.\n");
	printf("\t <n_out>: output number of bits from conditioning function.\n");
	printf("\t <nw>: narrowest internal width of conditioning function.\n");
	printf("\t <h_in>: input entropy to conditioning function.\n");
	printf("\t <-v|-n>: '-v' for vetted conditioning function, '-n' for non-vetted conditioning function. Vetted conditioning is the default.\n");
	printf("\t <h'>: entropy estimate per bit of conditioned sequential dataset (only for '-n' option).\n");
	printf("\n");
	printf("\t This program computes the entropy of the output of a conditioning function 'h_out' (Section 3.1.5).\n"); 
	printf("\t If the conditioning function is vetted, then\n\n"); 
	printf("\t\t h_out = Output_Entropy(n_in, n_out, nw, h_in)\n\n");
	printf("\t where 'Output_Entropy' is specified in Section 3.1.5.1.2. If the conditioning function is non-vetted then\n\n");
	printf("\t\t h_out = min(Output_Entropy(n_in, n_out, nw, h_in), 0.999*n_out, h'*n_out)\n\n");
	printf("\t as stated in Section 3.1.5.2.\n");
	printf("\n");
	exit(-1);
}

static long double inputLongDoubleOption(const char *input, long double low, long double high, const char *label) {
	char *nextoptchar;
	long double indouble;

	assert(input != NULL);
	assert(!isnan(low));
	assert(!isnan(high));

	if(label == NULL) label="parameter";

	indouble = strtold(input, &nextoptchar);
	assert(nextoptchar != NULL);

	if((nextoptchar == input) || (*nextoptchar != '\0')) {
		printf("Non-numeric characters in %s: '%c'\n", label, *nextoptchar);
		print_usage();
	}

	if((errno == ERANGE) || !isfinite(indouble)) {
		printf("Provided value for %s is out of range of a long double or isn't a finite value\n", label);
		print_usage();
	}

	if(indouble < low) {
		printf("%s must be greater than or equal to %.22Lg.\n", label, low);
		print_usage();
	}

	if(indouble > high) {
		printf("%s must be less than or equal to %.22Lg.\n", label, high);
		print_usage();
	}

	return indouble;
}

static unsigned int inputUnsignedOption(const char *input, unsigned int low, unsigned int high, const char *label) {
	char *nextoptchar;
	unsigned long inint;

	assert(input != NULL);

	if(label == NULL) label="parameter";

	inint = strtoul(input, &nextoptchar, 0);
	assert(nextoptchar != NULL);

	if((nextoptchar == input) || (*nextoptchar != '\0')) {
		printf("Non-integer characters in %s: '%c'\n", label, *nextoptchar);
		print_usage();
	}

	if(inint < (unsigned long int)low) {
		printf("%s must be greater than or equal to %u.\n", label, low);
		print_usage();
	}

	if(inint > (unsigned long int)high) {
		printf("%s must be less than or equal to %u.\n", label, high);
		print_usage();
	}

	return (unsigned int)inint;
}

//Check to see how close the provided value is to its maximal value
//1 - epsilon = value/max
// epsilon = 1 - value/max so -log2(epsilon) = -log2(1-value/max)
//To be conservative, round so that -log2(epsilon) epsilon is as small as possible
//(that is epsilon should be as large as possible)
static long double calculateEpsilon(mpfr_t calcValue, mpfr_t maxValue,  mpfr_prec_t precision) {
	mpfr_t ratio, output, ap_log2;

	mpfr_inits2(precision, ratio, output, ap_log2, NULL);

	//We're going to need an arbitrary precision version of log(2)
	mpfr_set_ui(ap_log2, 2U, MPFR_RNDZ);
	mpfr_log(ap_log2, ap_log2, MPFR_RNDU);

	mpfr_log_ui(ap_log2, 2UL, MPFR_RNDU);

	//Calculate the ratio value/max
	mpfr_set(ratio, calcValue, MPFR_RNDU);
	mpfr_neg(ratio, ratio, MPFR_RNDD);
	mpfr_div(ratio, ratio, maxValue, MPFR_RNDZ);

	//Calculate log(1-value/max)
	mpfr_log1p(output, ratio, MPFR_RNDZ);

	//Calculate log_2(1-value/max)
	mpfr_div(output, output, ap_log2, MPFR_RNDZ);

	//Calculate -log_2(1-value/max)
	mpfr_neg(output, output, MPFR_RNDZ);

	//return this value
	return mpfr_get_ld(output, MPFR_RNDZ);
}

int main(int argc, char* argv[]) {
	bool vetted;
	long double h_p = -1.0L;
	long double h_in, h_out; 
	long double outputEntropy = -1.0L;
	unsigned int n_in, n_out, nw, n;
	mpfr_prec_t precision;
	int opt;
	bool adaquatePrecision = false;
	unsigned int maxval;
	long double noutEpsilonExp = -1.0L;
	long double hinEpsilonExp = -1.0L;
	long double nwEpsilonExp = -1.0L;

	mpfr_t ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out;

	vetted = true;

	//Setting this rounding method helps prevent us from overestimating the input parameters
	fesetround(FE_TOWARDZERO);

        while ((opt = getopt(argc, argv, "vn")) != -1) {
                switch(opt) {
                        case 'v':
                                vetted = true;
                                break;
                        case 'n':
                                vetted = false;
                                break;
                        default:
                                print_usage();
                }
        }

        argc -= optind;
        argv += optind;

	// Parse args
	if((vetted && (argc != 4)) || (!vetted && (argc != 5))){
		printf("Incorrect usage.\n");
		print_usage();
	} else{
                // get n_in     
		n_in = inputUnsignedOption(argv[0], 1, UINT_MAX, "n_in");

	    	// get n_out     
                n_out = inputUnsignedOption(argv[1], 1, UINT_MAX, "n_out");

                // get nw     
                nw = inputUnsignedOption(argv[2], 1, UINT_MAX, "nw");

                // get h_in; note that h_in <= n_in
		h_in = inputLongDoubleOption(argv[3], 0.0L, (long double) n_in, "h_in");

		if(!vetted){
			h_p = inputLongDoubleOption(argv[4], 0.0L, 1.0L, "h_p");
			if(h_p <= 0) print_usage();
		} 
	}


	if(h_in <= 0.0L) print_usage();

	//Step 2 is invariant, and not subject to precision problems.
	if(nw > n_in) nw = n_in; //By 90B Appendix E
	n = std::min(n_out, nw);

	//Print out the inputs
	printf("n_in = %u\n", n_in);
	printf("n_out = %u\n", n_out);
	printf("nw = %u\n", nw);
	printf("h_in = %.22Lg\n", h_in);
	if(!vetted) printf("h' = %.22Lg\n", h_p);

	// Establish the maximum precision that ought to be necessary
	// If something goes wrong, we can increase this precision automatically.
	maxval = 64; // Always be large enough to faithfully represent h_in.
	maxval = (maxval>n_in)?maxval:n_in;
	maxval = (maxval>n_out)?maxval:n_out;
	maxval = (maxval>nw)?maxval:nw;
	precision = 2*maxval;

	//Initialize all the arbitrary precision values
	mpfr_inits2(precision, ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);

	adaquatePrecision = false;

	while(!adaquatePrecision) {
		//General goal: want to round to cause psi and omega to be as large as possible (to provide a conservative estimate)
		adaquatePrecision = true;
		//Initialize arbitrary precision versions of h_in
		//We want to make sure not to lose precision here.
		if(mpfr_set_ld(ap_h_in, h_in, MPFR_RNDZ) != 0) {
			adaquatePrecision=false;
		}

		// compute Output Entropy (Section 3.1.5.1.2)
		// Step 1.
		// Want to round so that both P_low and P_high are as large as possible.
		// P_high
		mpfr_neg(ap_entexp, ap_h_in, MPFR_RNDZ);
		mpfr_ui_pow(ap_p_high, 2UL, ap_entexp, MPFR_RNDU);

		// p_high must be in the interval (0,1)
		if(mpfr_cmp_ui(ap_p_high, 0UL)<=0) {
			adaquatePrecision=false;
		}
		if(mpfr_cmp_ui(ap_p_high, 1UL)>=0) {
			adaquatePrecision=false;
		}

		//P_low
		mpfr_ui_sub(ap_p_low, 1UL, ap_p_high, MPFR_RNDU); //p_low = 1-p_high

		//This is an integer value, and should be exact
		if(mpfr_ui_pow_ui (ap_denom, 2UL, n_in, MPFR_RNDZ)!=0) { //ap_denom = 2^(n_in)
			adaquatePrecision=false;
		}

		mpfr_sub_ui(ap_denom, ap_denom, 1UL, MPFR_RNDZ); // ap_denom = 2^(n_in) - 1
		mpfr_div (ap_p_low, ap_p_low, ap_denom, MPFR_RNDU); // p_low = (1-p_high)/(2^(n_in)-1)
		// p_low must be in the interval (0,1)
		if(mpfr_cmp_ui(ap_p_low, 0UL)<=0) {
			adaquatePrecision=false;
		}
		if(mpfr_cmp_ui(ap_p_low, 1UL)>=0) {
			adaquatePrecision=false;
		}

		//Prior to moving on, calculate a reused power term
		//This is an integer value, and should be exact
		if(mpfr_ui_pow_ui(ap_power_term, 2UL, n_in - n, MPFR_RNDU)!=0) { //power_term = 2^(n_in - n)
			adaquatePrecision=false;
		}

		//Step 3: Calculate Psi
		mpfr_mul(ap_psi, ap_power_term, ap_p_low, MPFR_RNDU);
		mpfr_add(ap_psi, ap_psi, ap_p_high,  MPFR_RNDU);
		// h_in > 0 so Psi > P_high. If this isn't so, then we're doing the calculation at too low of a precision.
		if(mpfr_cmp(ap_p_high, ap_psi) >= 0) {
			adaquatePrecision=false;
		}

		assert(mpfr_cmp_ui(ap_psi, 0UL)>=0);

		if(mpfr_cmp_ui(ap_psi, 0UL)==0) {
			adaquatePrecision=false;
		}

		//Is psi > 1?
		if(mpfr_cmp_ui(ap_psi, 1UL)>0) {
			//This is quite unlikely for most parameters, but it is possible for some allowed values.
			//Set this value to the largest meaningful value.
			mpfr_set_ui(ap_psi, 1UL, MPFR_RNDZ);
		}

		//We're going to need an arbitrary precision version of log(2)
		if(mpfr_set_ui(ap_omega, 2U, MPFR_RNDZ) != 0) { //omega = 2
			adaquatePrecision=false;
		}
		mpfr_log(ap_omega, ap_omega, MPFR_RNDU); // omega = log(2)

		//Step 4: Calculate U (goes into the ap_omega variable)
		mpfr_mul(ap_omega, ap_omega, ap_power_term, MPFR_RNDU); //omega = log(2) 2^(n_in - n)
		mpfr_mul_ui(ap_omega, ap_omega, 2UL*n, MPFR_RNDU); //omega = log(2) 2^(n_in - n) * 2 * n
		mpfr_sqrt(ap_omega, ap_omega, MPFR_RNDU); // omega = Sqrt(log(2) 2^(n_in - n) * 2 * n)
		mpfr_add(ap_omega, ap_omega, ap_power_term,  MPFR_RNDU); // omega = 2^(n_in-n) + Sqrt(log(2) 2^(n_in - n) * 2 * n)

		//Step 5: Calculate omega
		mpfr_mul(ap_omega, ap_omega, ap_p_low, MPFR_RNDU); // omega = (2^(n_in-n) + Sqrt(log(2) 2^(n_in - n) * 2 * n)) * p_low
		//Omega is expected to be non-negative
		assert(mpfr_cmp_ui(ap_omega, 0UL)>=0);

		if(mpfr_cmp_ui(ap_omega, 0UL)==0) {
			//Omega is expected to be non-zero for all parameters
			adaquatePrecision=false;
		}

		//Is omega > 1?
		if(mpfr_cmp_ui(ap_omega, 1UL)>0) {
			//This is quite unlikely for most parameters, but it is possible for some allowed values.
			//Set this value to the largest meaningful value.
			mpfr_set_ui(ap_omega, 1UL, MPFR_RNDZ);
		}

		//Step 6: Compare the values
		//we want to round so that the log is (in absolute value) as small as possible.
		if(mpfr_cmp(ap_omega, ap_psi) > 0) {
			//omega > psi
			mpfr_log2(ap_outputEntropy, ap_omega, MPFR_RNDZ);
		} else {
			//omega <= psi
			mpfr_log2(ap_outputEntropy, ap_psi, MPFR_RNDZ);
		}

		//finalize outputEntropy
		mpfr_neg(ap_outputEntropy, ap_outputEntropy, MPFR_RNDZ);

		//Could outputEntropy be valid?
		//We know that n_out > ap_outputEntropy for all finite inputs...
		if(mpfr_cmp_ui(ap_outputEntropy, n_out) >= 0) {
			adaquatePrecision = false;
		}

		//We know that h_in > ap_outputEntropy for all finite inputs...
		if(mpfr_cmp(ap_outputEntropy, ap_h_in) >= 0) {
			adaquatePrecision = false;
		}

		if(adaquatePrecision) {
			//Check to see if meets the definition of "full entropy".
			//iff -log2(epsilon) = - log(1 - (h_out)/n_out)/log(2) > 64 (2012 90B draft) or > 32 (2021 90C draft)
			//To be conservative, round so that -log2(epsilon) epsilon is as small as possible
			//(that is epsilon should be as large as possible)
			mpfr_set_ui(ap_n_out, n_out, MPFR_RNDZ);
			noutEpsilonExp = calculateEpsilon(ap_outputEntropy, ap_n_out, precision);

			//We may also be interested in other ways that this output may have been limited.
			hinEpsilonExp = calculateEpsilon(ap_outputEntropy, ap_h_in, precision);

			mpfr_set_ui(ap_nw, nw, MPFR_RNDZ);
			nwEpsilonExp = calculateEpsilon(ap_outputEntropy, ap_nw, precision);

			//If we get here, then adequate precision was used
			//Extract a value for display.
			//Note, this may round up, but we'll deal with this later.
			outputEntropy = mpfr_get_ld(ap_outputEntropy, MPFR_RNDN);
		} else {
			//We either inappropriately rounded previously,
			//or outputEntropy is greater than or equal to n_out
			//or outputEntropy is greater than or equal to h_in
			//In any case, there was a precision problem.
			//Double the precision of everything and try again
			precision = 2*precision;
			fprintf(stderr, "Increasing precision to %ld bits and trying again.\n", precision);
			mpfr_set_prec(ap_h_in, precision);
			mpfr_set_prec(ap_entexp, precision);
			mpfr_set_prec(ap_p_high, precision);
			mpfr_set_prec(ap_p_low, precision);
			mpfr_set_prec(ap_denom, precision);
			mpfr_set_prec(ap_power_term, precision);
			mpfr_set_prec(ap_psi, precision);
			mpfr_set_prec(ap_omega, precision);
			mpfr_set_prec(ap_outputEntropy, precision);
			mpfr_set_prec(ap_nw, precision);
			mpfr_set_prec(ap_n_out, precision);
		}
	}

	//Check some basic bounds.
	assert(outputEntropy <= (long double)n_out);
	assert(outputEntropy <= h_in);
	assert(outputEntropy <= (long double)nw);
	assert(outputEntropy >= 0.0L);

	//We're done with the calculation. Now print results.
	printf("Output_Entropy(*) = %.22Lg", outputEntropy);
	if(outputEntropy == (long double)n_out) {
		//outputEntropy rounded to full entropy, so the difference between this and full entropy is less than 1/2 ULP.
		printf("; Close to n_out (epsilon = 2^(-%.22Lg))", noutEpsilonExp);
	}
	if(outputEntropy == h_in) {
		//outputEntropy rounded to the input entropy, so the difference between this and the input entropy is less than 1/2 ULP.
		printf("; Close to h_in (epsilon = 2^(-%.22Lg))", hinEpsilonExp);
	}
	if(outputEntropy == (long double)nw) {
		//outputEntropy rounded to the nw, so the difference between this and nw is less than 1/2 ULP.
		printf("; Close to nw (epsilon = 2^(-%.22Lg))", nwEpsilonExp);
	}
	printf("\n");

	if(vetted) {
		printf("(Vetted) h_out = %.22Lg\n", outputEntropy);

		if(outputEntropy > 0.999L*((long double)n_out)) {
			//h_out = (1 - epsilon) * n_out
			printf("epsilon = 2^(-%.22Lg)", noutEpsilonExp);
			//Should this qualify as "full entropy" under the 2012 draft of SP 800-90B?
			//Should this qualify as "full entropy" under the 2021 SP 800-90C draft?
			if(noutEpsilonExp >= 64.0L) {
				printf(": SP 800-90B 2012 Draft and SP 800-90C 2021 Draft Full Entropy");
			} else if(noutEpsilonExp >= 32.0L) {
				printf(": SP 800-90C 2021 Full Entropy");
			}
			printf("\n");
		}
	} else {
		long double bound90B = 0.999L*((long double)n_out);
		long double statBound = h_p*((long double)n_out);

		//Note, we can't assess as full entropy in this case.
		printf("0.999 * n_out = %.22Lg\n", bound90B);
		printf("h' * n_out = %.22Lg\n", statBound);
		h_out = std::min(outputEntropy, std::min(bound90B, statBound));
		printf("(Non-vetted) h_out = %.22Lg\n", h_out);
	}

	return 0;
}
