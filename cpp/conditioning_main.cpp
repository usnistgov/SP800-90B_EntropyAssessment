#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <algorithm>

#include <getopt.h>
#include <mpfr.h>

[[ noreturn ]] void print_usage() {
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

int main(int argc, char* argv[]) {
	bool vetted;
	long double h_p = -1.0;
	long double h_in, h_out; 
	long double outputEntropy;
	unsigned int n_in, n_out, nw, n;
	mpfr_prec_t precision;
	int opt;
	bool adaquatePrecision = false;
	bool closeToFullEntropy = false;
	bool fullEntropy = false;
	unsigned int maxval;
	long double epsilonExp = -1.0;
	mpfr_t ap_h_in, ap_p_high, ap_p_low, ap_denom, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_log2epsilon;

	vetted = true;

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
	}
	else{
                // get n_in     
                n_in = strtoul(argv[0], NULL, 0);
		if(n_in < 1) {
                        printf("n_in must be greater than 0.\n");
                        print_usage();
		}

	    	// get n_out     
                n_out = strtoul(argv[1], NULL, 0);

                // get n_out     
                nw = strtoul(argv[2], NULL, 0);

                // get h_in     
                h_in = strtold(argv[3], NULL);
                if((h_in <= 0) || (h_in > (long double)n_in)){
                        printf("h_in must be positive and at most n_in.\n");
                        print_usage();
                }

		if(!vetted){
			h_p = strtold(argv[4], NULL);
			if((h_p <= 0) || (h_p > 1.0)){
				printf("h' must be the interval (0 and 1].\n");
				print_usage();
			}
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
	maxval = (n_in>n_out)?n_in:n_out;
	maxval = (maxval>nw)?maxval:nw;
	precision = 2*maxval;

	//Initialize all the arbitrary precision values
	mpfr_inits2(precision, ap_h_in, ap_p_high, ap_p_low, ap_denom, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_log2epsilon, NULL);

	adaquatePrecision = false;

	while(!adaquatePrecision) {
		adaquatePrecision = true;
		//Initialize arbitrary precision versions of h_in and needed constants.
		mpfr_set_ld(ap_h_in, h_in, MPFR_RNDN);

		// compute Output Entropy (Section 3.1.5.1.2)
		// Step 1.
		// P_high
		mpfr_neg(ap_h_in, ap_h_in, MPFR_RNDN);
		mpfr_ui_pow(ap_p_high, 2UL, ap_h_in, MPFR_RNDN);

		//P_low
		mpfr_ui_sub(ap_p_low, 1UL, ap_p_high, MPFR_RNDN);
		mpfr_ui_pow_ui (ap_denom, 2UL, n_in, MPFR_RNDN);
		mpfr_sub_ui(ap_denom, ap_denom, 1UL, MPFR_RNDN);
		mpfr_div (ap_p_low, ap_p_low, ap_denom, MPFR_RNDN);

		//Prior to moving on, calculate a reused power term
		mpfr_ui_pow_ui(ap_power_term, 2UL, n_in - n, MPFR_RNDN);

		//Step 3: Calculate Psi
		mpfr_mul(ap_psi, ap_power_term, ap_p_low, MPFR_RNDN);
		mpfr_add(ap_psi, ap_psi, ap_p_high,  MPFR_RNDN);

		//Step 4: Calculate U (goes into the ap_omega variable)
		mpfr_log_ui(ap_omega, 2UL, MPFR_RNDN);
		mpfr_mul(ap_omega, ap_omega, ap_power_term, MPFR_RNDN);
		mpfr_mul_ui(ap_omega, ap_omega, 2UL*n, MPFR_RNDN);
		mpfr_sqrt(ap_omega, ap_omega, MPFR_RNDN);
		mpfr_add(ap_omega, ap_omega, ap_power_term,  MPFR_RNDN);

		//Step 5: Calculate omega
		mpfr_mul(ap_omega, ap_omega, ap_p_low, MPFR_RNDN);

		//Step 6: Compare the values
		if(mpfr_cmp(ap_omega, ap_psi) > 0) {
			//omega > psi
			mpfr_log2(ap_outputEntropy, ap_omega, MPFR_RNDN);
		} else {
			//omega <= psi
			mpfr_log2(ap_outputEntropy, ap_psi, MPFR_RNDN);
		}

		//Keep -outputEntropy for calculation of log2epsilon
		mpfr_set(ap_log2epsilon, ap_outputEntropy, MPFR_RNDN);

		//finalize outputEntropy
		mpfr_neg(ap_outputEntropy, ap_outputEntropy, MPFR_RNDN);


		if(mpfr_cmp_ui(ap_outputEntropy, n_out) >= 0) {
			//Double the precision of everything and try again
			adaquatePrecision = false;
			precision = 2*precision;
			fprintf(stderr, "Indeterminate result. Increasing precision to %ld bits and trying again.\n", precision);
			mpfr_set_prec(ap_h_in, precision);
			mpfr_set_prec(ap_p_high, precision);
			mpfr_set_prec(ap_p_low, precision);
			mpfr_set_prec(ap_denom, precision);
			mpfr_set_prec(ap_power_term, precision);
			mpfr_set_prec(ap_psi, precision);
			mpfr_set_prec(ap_omega, precision);
			mpfr_set_prec(ap_outputEntropy, precision);
			mpfr_set_prec(ap_log2epsilon, precision);
			continue;
		}

		outputEntropy = mpfr_get_ld(ap_outputEntropy, MPFR_RNDN);

		if(outputEntropy > 0.999 * (long double)n_out) {
			closeToFullEntropy = true;
			//Calculate -log2(epsilon)
			mpfr_div_ui(ap_log2epsilon, ap_log2epsilon, n_out, MPFR_RNDN);
			mpfr_log1p(ap_log2epsilon, ap_log2epsilon, MPFR_RNDN);
			//Reuse ap_denom to hold log(2)
			mpfr_log_ui(ap_denom, 2UL, MPFR_RNDN);
			mpfr_div(ap_log2epsilon, ap_log2epsilon, ap_denom, MPFR_RNDN);
			//Make it positive
			mpfr_neg(ap_log2epsilon, ap_log2epsilon, MPFR_RNDN);

			//Now convert back to a long double
			epsilonExp = mpfr_get_ld(ap_log2epsilon, MPFR_RNDN);

			//Should this qualify as "full entropy" under the 2012 draft of SP800-90B?
			if(mpfr_cmp_ui(ap_log2epsilon, 64) >= 0) {
				fullEntropy = true;
			}
		}
	}

	//We're done with the calculation. Now print results.
	if(vetted) {
		printf("\n(Vetted) ");
		if(closeToFullEntropy) {
			if(outputEntropy == (long double)n_out) {
				printf("h_out: Close to %.22Lg\n", outputEntropy);
			} else {
				printf("h_out: %.22Lg\n", outputEntropy);
			}
			printf("epsilon = 2^(-%.22Lg)", epsilonExp);
			if(fullEntropy) {
				printf(": SP800-90B 2012 Full Entropy\n");
			} else {
				printf("\n");
			}
		} else {
			printf("h_out: %.22Lg\n", outputEntropy);
		}
	} else {
		printf("\n(Non-vetted) ");
		h_out = std::min(outputEntropy, std::min(0.999L*((long double)n_out), h_p*(long double)n_out));
		printf("h_out: %.22Lg\n", h_out);
	}

	return 0;
}
