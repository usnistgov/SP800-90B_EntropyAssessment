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

int main(int argc, char* argv[]) {
	bool vetted;
	long double h_p = -1.0L;
	long double h_in, h_out; 
	long double outputEntropy = -1.0L;
	unsigned int n_in, n_out, nw, n;
	mpfr_prec_t precision;
	int opt;
	bool adaquatePrecision = false;
	bool closeToFullEntropy = false;
	bool fullEntropy = false;
	unsigned int maxval;
	long double epsilonExp = -1.0L;

	mpfr_t ap_h_in, ap_p_high, ap_p_low, ap_denom, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_log2epsilon, ap_log2;

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

                // get h_in     
		h_in = inputLongDoubleOption(argv[3], 0.0L, (long double) n_in, "h_in");

		if(!vetted){
			h_p = inputLongDoubleOption(argv[4], 0.0L, 1.0L, "h_p");
		} 
	}


	//Step 2 is invariant, and not subject to precision problems.
	if(nw > n_in) nw = n_in; //By 90B Appendix E
	n = std::min(n_out, nw);

	//Print out the inputs
	printf("n_in: %u\n", n_in);
	printf("n_out: %u\n", n_out);
	printf("nw: %u\n", nw);
	printf("h_in: %.22Lg\n", h_in);
	if(!vetted) printf("h': %.22Lg\n", h_p);

	// Establish the maximum precision that ought to be necessary
	// If something goes wrong, we can increase this precision automatically.
	maxval = 64; // Always be large enough to faithfully represent h_in.
	maxval = (maxval>n_in)?maxval:n_in;
	maxval = (maxval>n_out)?maxval:n_out;
	maxval = (maxval>nw)?maxval:nw;
	precision = 2*maxval;

	//Initialize all the arbitrary precision values
	mpfr_inits2(precision, ap_h_in, ap_p_high, ap_p_low, ap_denom, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_log2epsilon, ap_log2, NULL);

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
		mpfr_neg(ap_h_in, ap_h_in, MPFR_RNDZ);
		mpfr_ui_pow(ap_p_high, 2UL, ap_h_in, MPFR_RNDU);

		//P_low
		mpfr_ui_sub(ap_p_low, 1UL, ap_p_high, MPFR_RNDU);

		//This is an integer value, and should be exact
		if(mpfr_ui_pow_ui (ap_denom, 2UL, n_in, MPFR_RNDZ)!=0) {
			adaquatePrecision=false;
		}

		mpfr_sub_ui(ap_denom, ap_denom, 1UL, MPFR_RNDZ);
		mpfr_div (ap_p_low, ap_p_low, ap_denom, MPFR_RNDU);

		//Prior to moving on, calculate a reused power term
		//This is an integer value, and should be exact
		if(mpfr_ui_pow_ui(ap_power_term, 2UL, n_in - n, MPFR_RNDU)!=0) {
			adaquatePrecision=false;
		}

		//Step 3: Calculate Psi
		mpfr_mul(ap_psi, ap_power_term, ap_p_low, MPFR_RNDU);
		mpfr_add(ap_psi, ap_psi, ap_p_high,  MPFR_RNDU);

		//We're going to need an arbitrary precision version of log(2)
		if(mpfr_set_ui(ap_log2, 2U, MPFR_RNDZ) != 0) {
			adaquatePrecision=false;
		}
		mpfr_log(ap_log2, ap_log2, MPFR_RNDU);

		//Step 4: Calculate U (goes into the ap_omega variable)
		mpfr_set(ap_omega, ap_log2, MPFR_RNDU);
		mpfr_mul(ap_omega, ap_omega, ap_power_term, MPFR_RNDU);
		mpfr_mul_ui(ap_omega, ap_omega, 2UL*n, MPFR_RNDU);
		mpfr_sqrt(ap_omega, ap_omega, MPFR_RNDU);
		mpfr_add(ap_omega, ap_omega, ap_power_term,  MPFR_RNDU);

		//Step 5: Calculate omega
		mpfr_mul(ap_omega, ap_omega, ap_p_low, MPFR_RNDU);

		//Step 6: Compare the values
		//we want to round so that the log is (in absolute value) as small as possible.
		if(mpfr_cmp(ap_omega, ap_psi) > 0) {
			//omega > psi
			mpfr_log2(ap_outputEntropy, ap_omega, MPFR_RNDZ);
		} else {
			//omega <= psi
			mpfr_log2(ap_outputEntropy, ap_psi, MPFR_RNDZ);
		}

		//Keep -outputEntropy for calculation of log2epsilon
		mpfr_set(ap_log2epsilon, ap_outputEntropy, MPFR_RNDZ);

		//finalize outputEntropy
		mpfr_neg(ap_outputEntropy, ap_outputEntropy, MPFR_RNDZ);

		//Could outputEntropy be valid?
		//We know that n_out > ap_outputEntropy for all finite inputs...
		if(mpfr_cmp_ui(ap_outputEntropy, n_out) >= 0) {
			adaquatePrecision = false;
		}

		if(adaquatePrecision) {
			//If we get here, then adequate precision was used
			//Extract a value for display.
			//Note, this may round up, but we'll deal with this later.
			outputEntropy = mpfr_get_ld(ap_outputEntropy, MPFR_RNDN);

			//Check to see if meets the definition of "full entropy".
			//iff -log2(epsilon) = - log(1 - (h_out)/n_out)/log(2) > 64
			//To be conservative, round so that -log2(epsilon) epsilon is as small as possible
			//(that is epsilon should be as large as possible)
			if(outputEntropy > 0.999L * (long double)n_out) {
				closeToFullEntropy = true;
				//Calculate -log2(epsilon)
				mpfr_div_ui(ap_log2epsilon, ap_log2epsilon, n_out, MPFR_RNDZ);
				mpfr_log1p(ap_log2epsilon, ap_log2epsilon, MPFR_RNDZ);
				//Convert the result to log base 2
				mpfr_div(ap_log2epsilon, ap_log2epsilon, ap_log2, MPFR_RNDZ);
				//Make the result positive
				mpfr_neg(ap_log2epsilon, ap_log2epsilon, MPFR_RNDZ);

				//We now have a (arbitrary precision) version of the exponent.
				//Convert it back to a long double...
				epsilonExp = mpfr_get_ld(ap_log2epsilon, MPFR_RNDZ);

				//Should this qualify as "full entropy" under the 2012 draft of SP800-90B?
				//Compare using the arbitrary precision version of the exponent.
				if(mpfr_cmp_ui(ap_log2epsilon, 64) >= 0) {
					fullEntropy = true;
				}
			}
		} else {
			//We either inappropriately rounded previously, or
			//outputEntropy is greater than or equal to n_out
			//In either case, there was a precision problem.
			//Double the precision of everything and try again
			precision = 2*precision;
			fprintf(stderr, "Increasing precision to %ld bits and trying again.\n", precision);
			mpfr_set_prec(ap_h_in, precision);
			mpfr_set_prec(ap_p_high, precision);
			mpfr_set_prec(ap_p_low, precision);
			mpfr_set_prec(ap_denom, precision);
			mpfr_set_prec(ap_power_term, precision);
			mpfr_set_prec(ap_psi, precision);
			mpfr_set_prec(ap_omega, precision);
			mpfr_set_prec(ap_outputEntropy, precision);
			mpfr_set_prec(ap_log2epsilon, precision);
			mpfr_set_prec(ap_log2, precision);
		}
	}

	assert((outputEntropy <= (long double)n_out) && (outputEntropy >= 0.0L));

	//We're done with the calculation. Now print results.
	if(vetted) {
		printf("(Vetted) ");
		if(closeToFullEntropy) {
			if(outputEntropy == (long double)n_out) {
				//outputEntropy rounded to full entropy, so the difference between this and full entropy is less than 1/2 ULP.
				printf("h_out: Close to %u\n", n_out);
			} else {
				printf("h_out: %.22Lg\n", outputEntropy);
			}

			//h_out = (1 - epsilon) * n_out
			printf("epsilon: 2^(-%.22Lg)", epsilonExp);
			if(fullEntropy) {
				printf(": SP800-90B 2012 Full Entropy\n");
			} else {
				printf("\n");
			}
		} else {
			printf("h_out: %.22Lg\n", outputEntropy);
		}
	} else {
		//Note, we can't assess as full entropy in this case.
		printf("(Non-vetted) ");
		h_out = std::min(outputEntropy, std::min(0.999L*((long double)n_out), h_p*(long double)n_out));
		printf("h_out: %.22Lg\n", h_out);
	}

	return 0;
}
