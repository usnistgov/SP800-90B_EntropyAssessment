/* VERSION information is kept in utils.h. Please update when a new version is released */

#include "non_iid/collision_test.h"
#include "non_iid/lz78y_test.h"
#include "non_iid/multi_mmc_test.h"
#include "non_iid/lag_test.h"
#include "non_iid/multi_mcw_test.h"
#include "non_iid/compression_test.h"
#include "non_iid/markov_test.h"
#include "non_iid/non_iid_test_run.h"
#include "iid/iid_test_run.h"
#include "shared/TestRunUtils.h"
#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include <string.h>
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
#include <iostream>
#include <fstream>

[[ noreturn ]] void print_usage() {
    printf("Usage is: ea_conditioning -v [-q] <n_in> <n_out> <nw> <h_in> [-o filename.json]\n");
    printf("\tor \n\tea_conditioning -n <n_in> <n_out> <nw> <h_in> [h' | -i filename] [-o filename.json]\n\n");
    printf("\t <n_in>: input number of bits to conditioning function.\n");
    printf("\t <n_out>: output number of bits from conditioning function.\n");
    printf("\t <nw>: narrowest internal width of conditioning function.\n");
    printf("\t <h_in>: input entropy to conditioning function.\n");
    printf("\t <-v|-n>: '-v' for vetted conditioning function, '-n' for non-vetted conditioning function. Vetted conditioning is the default.\n");
    printf("\t <h'>: entropy estimate per bit of conditioned sequential dataset (only for '-n' option).\n");
    printf("\t -q: Quiet mode, less output to screen.\n");
    printf("\t -i: Input file name, to run an entropy assessment on a non-vetted conditioned data file and use that value as h'.\n");
    printf("\n");
    printf("\t This program computes the entropy of the output of a conditioning function 'h_out' (Section 3.1.5).\n");
    printf("\t If the conditioning function is vetted, then\n\n");
    printf("\t\t h_out = Output_Entropy(n_in, n_out, nw, h_in)\n\n");
    printf("\t where 'Output_Entropy' is specified in Section 3.1.5.1.2. If the conditioning function is non-vetted then\n\n");
    printf("\t\t h_out = min(Output_Entropy(n_in, n_out, nw, h_in), 0.999*n_out, h'*n_out)\n\n");
    printf("\t as stated in Section 3.1.5.2.\n");
    printf("\n");
    printf("\t -o: Set Output Type to JSON\n");
    printf("\n");
    printf("\t\t Changes the output format to JSON and sets the file location for the output file.\n");
    printf("\n");
    printf("\t --version: Prints tool version information");
    printf("\n");
    exit(-1);
}

static long double inputLongDoubleOption(const char *input, long double low, long double high, const char *label) {
    char *nextoptchar;
    long double indouble;

    assert(input != NULL);
    assert(!isnan(low));
    assert(!isnan(high));

    if (label == NULL) label = "parameter";

    indouble = strtold(input, &nextoptchar);
    assert(nextoptchar != NULL);

    if ((nextoptchar == input) || (*nextoptchar != '\0')) {
        printf("Non-numeric characters in %s: '%c'\n", label, *nextoptchar);
        print_usage();
    }

    if ((errno == ERANGE) || !isfinite(indouble)) {
        printf("Provided value for %s is out of range of a long double or isn't a finite value\n", label);
        print_usage();
    }

    if (indouble < low) {
        printf("%s must be greater than or equal to %.22Lg.\n", label, low);
        print_usage();
    }

    if (indouble > high) {
        printf("%s must be less than or equal to %.22Lg.\n", label, high);
        print_usage();
    }

    return indouble;
}

static unsigned int inputUnsignedOption(const char *input, unsigned int low, unsigned int high, const char *label) {
    char *nextoptchar;
    unsigned long inint;

    assert(input != NULL);

    if (label == NULL) label = "parameter";

    inint = strtoul(input, &nextoptchar, 0);
    assert(nextoptchar != NULL);

    if ((nextoptchar == input) || (*nextoptchar != '\0')) {
        printf("Non-integer characters in %s: '%c'\n", label, *nextoptchar);
        print_usage();
    }

    if (inint < (unsigned long int) low) {
        printf("%s must be greater than or equal to %u.\n", label, low);
        print_usage();
    }

    if (inint > (unsigned long int) high) {
        printf("%s must be less than or equal to %u.\n", label, high);
        print_usage();
    }

    return (unsigned int) inint;
}

/*
 * Check to see how close the provided value is to its maximal value
 * 1 - epsilon = value/max  =>  epsilon = 1 - value/max  =>  -log2(epsilon) = -log2(1 - value/max)
 * To be conservative, round so that -log2(epsilon) epsilon is as small as possible
 * (that is epsilon should be as large as possible)
 */
static long double calculateEpsilon(mpfr_t calcValue, mpfr_t maxValue, mpfr_prec_t precision) {
    mpfr_t ratio, output, ap_log2;
    mpfr_inits2(precision, ratio, output, ap_log2, NULL);
    long double value;

    // We're going to need an arbitrary precision version of log(2)
    mpfr_set_ui(ap_log2, 2U, MPFR_RNDZ);
    mpfr_log(ap_log2, ap_log2, MPFR_RNDU);

    // Calculate the ratio value/max
    mpfr_set(ratio, calcValue, MPFR_RNDU);
    mpfr_neg(ratio, ratio, MPFR_RNDD);
    mpfr_div(ratio, ratio, maxValue, MPFR_RNDZ);

    // Calculate log(1 - value/max)
    mpfr_log1p(output, ratio, MPFR_RNDZ);

    // Calculate log_2(1 - value/max)
    mpfr_div(output, output, ap_log2, MPFR_RNDZ);

    // Calculate -log_2(1 - value/max)
    mpfr_neg(output, output, MPFR_RNDZ);

    // return this value
    value = mpfr_get_ld(output, MPFR_RNDZ);
    mpfr_clears(ratio, output, ap_log2, NULL);

    return value;
}

// General goal: want to round to cause psi and omega to be as large as possible (to provide a conservative estimate)
// If any estimate is not appropriate, increase the precision and start again

static long double computeEntropyWithPrecision(mpfr_prec_t precision, long double h_in, unsigned int n_in, unsigned int n, unsigned int n_out, unsigned int nw, long double &noutEpsilonExp, long double &hinEpsilonExp, long double &nwEpsilonExp) {
    long double value;

    // TODO quietmode?
    printf("Attempting to compute entropy with %ld bits of precision.\n", precision);

    // Initialize all the arbitrary precision values
    mpfr_t ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out;
    mpfr_inits2(precision, ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);

    // Initialize arbitrary precision versions of h_in
    // We want to make sure not to lose precision here.
    if (mpfr_set_ld(ap_h_in, h_in, MPFR_RNDZ) != 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // Compute Output Entropy (Section 3.1.5.1.2)
    // Step 1.
    // Want to round so that both P_low and P_high are as large as possible.
    // p_high
    mpfr_neg(ap_entexp, ap_h_in, MPFR_RNDZ);
    mpfr_ui_pow(ap_p_high, 2UL, ap_entexp, MPFR_RNDU);

    // p_high must be in the interval (0,1)
    if (mpfr_cmp_ui(ap_p_high, 0UL) <= 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    if (mpfr_cmp_ui(ap_p_high, 1UL) >= 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // p_low = 1 - p_high
    mpfr_ui_sub(ap_p_low, 1UL, ap_p_high, MPFR_RNDU);

    // This is an integer value, and should be exact
    // ap_inputSpaceSize = 2^(n_in)
    if (mpfr_ui_pow_ui(ap_inputSpaceSize, 2UL, n_in, MPFR_RNDZ) != 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // ap_denom = 2^(n_in) - 1
    mpfr_sub_ui(ap_denom, ap_inputSpaceSize, 1UL, MPFR_RNDZ);

    // Is the difference correct?
    mpfr_sub(ap_diff, ap_inputSpaceSize, ap_denom, MPFR_RNDZ);
    if (mpfr_cmp_ui(ap_diff, 1UL) != 0) {
        // Evidently not. Increase the precision.
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // p_low = (1-p_high)/(2^(n_in)-1)
    mpfr_div(ap_p_low, ap_p_low, ap_denom, MPFR_RNDU);

    // p_low must be in the interval (0,1)
    if (mpfr_cmp_ui(ap_p_low, 0UL) <= 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    if (mpfr_cmp_ui(ap_p_low, 1UL) >= 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // Prior to moving on, calculate a reused power term
    // This is an integer value, and should be exact
    // power_term = 2^(n_in - n)
    if (mpfr_ui_pow_ui(ap_power_term, 2UL, n_in - n, MPFR_RNDU) != 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // Step 3: Calculate Psi
    // ap_psi = 2^(n_in - n) * p_low
    mpfr_mul(ap_psi, ap_power_term, ap_p_low, MPFR_RNDU);

    // ap_psi = 2^(n_in - n) * p_low + p_high
    mpfr_add(ap_psi, ap_psi, ap_p_high, MPFR_RNDU);

    // h_in > 0 so Psi > P_high. If this isn't so, then we're doing the calculation at too low of a precision.
    if (mpfr_cmp(ap_p_high, ap_psi) >= 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // Psi > 0 is expected
    assert(mpfr_cmp_ui(ap_psi, 0UL) >= 0);

    // If we have equality, then we didn't use adaquate precision.
    if (mpfr_cmp_ui(ap_psi, 0UL) == 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // Is psi > 1?
    if (mpfr_cmp_ui(ap_psi, 1UL) > 0) {
        // This is quite unlikely for most parameters, but it is possible for some allowed values.
        // Set this value to the largest meaningful value.
        mpfr_set_ui(ap_psi, 1UL, MPFR_RNDZ);
    }

    // We're going to need an arbitrary precision version of log(2)
    // omega = 2
    if (mpfr_set_ui(ap_omega, 2U, MPFR_RNDZ) != 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // omega = log(2)
    mpfr_log(ap_omega, ap_omega, MPFR_RNDU);

    // Step 4: Calculate U (goes into the ap_omega variable)
    mpfr_mul(ap_omega, ap_omega, ap_power_term, MPFR_RNDU); //omega = log(2) 2^(n_in - n)
    mpfr_mul_ui(ap_omega, ap_omega, 2UL * n, MPFR_RNDU); //omega = log(2) 2^(n_in - n) * 2 * n
    mpfr_sqrt(ap_omega, ap_omega, MPFR_RNDU); // omega = Sqrt(log(2) 2^(n_in - n) * 2 * n)
    mpfr_add(ap_omega, ap_omega, ap_power_term, MPFR_RNDU); // omega = 2^(n_in-n) + Sqrt(log(2) 2^(n_in - n) * 2 * n)

    // Step 5: Calculate omega
    mpfr_mul(ap_omega, ap_omega, ap_p_low, MPFR_RNDU); // omega = (2^(n_in-n) + Sqrt(log(2) 2^(n_in - n) * 2 * n)) * p_low

    // Omega is expected to be non-negative
    assert(mpfr_cmp_ui(ap_omega, 0UL) >= 0);

    if (mpfr_cmp_ui(ap_omega, 0UL) == 0) {
        // Omega is expected to be non-zero for all parameters
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // Is omega > 1?
    if (mpfr_cmp_ui(ap_omega, 1UL) > 0) {
        // This is quite unlikely for most parameters, but it is possible for some allowed values.
        // Set this value to the largest meaningful value.
        mpfr_set_ui(ap_omega, 1UL, MPFR_RNDZ);
    }

    // Step 6: Compare the values
    // We want to round so that the log is (in absolute value) as small as possible.
    if (mpfr_cmp(ap_omega, ap_psi) > 0) {
        // omega > psi
        mpfr_log2(ap_outputEntropy, ap_omega, MPFR_RNDZ);
    } else {
        // omega <= psi
        mpfr_log2(ap_outputEntropy, ap_psi, MPFR_RNDZ);
    }

    // Finalize outputEntropy
    mpfr_neg(ap_outputEntropy, ap_outputEntropy, MPFR_RNDZ);

    // Could outputEntropy be valid?
    // We know that n_out > ap_outputEntropy for all finite inputs...
    if (mpfr_cmp_ui(ap_outputEntropy, n_out) >= 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    //We know that h_in > ap_outputEntropy for all finite inputs...
    if (mpfr_cmp(ap_outputEntropy, ap_h_in) >= 0) {
        mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);
        return computeEntropyWithPrecision(precision * 2, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);
    }

    // Check to see if meets the definition of "full entropy".
    // iff -log2(epsilon) = - log(1 - (h_out)/n_out)/log(2) > 64 (2012 90B draft) or > 32 (2021 90C draft)
    // To be conservative, round so that -log2(epsilon) epsilon is as small as possible
    // (that is epsilon should be as large as possible)
    mpfr_set_ui(ap_n_out, n_out, MPFR_RNDZ);
    noutEpsilonExp = calculateEpsilon(ap_outputEntropy, ap_n_out, precision);

    // We may also be interested in other ways that this output may have been limited.
    hinEpsilonExp = calculateEpsilon(ap_outputEntropy, ap_h_in, precision);

    mpfr_set_ui(ap_nw, nw, MPFR_RNDZ);
    nwEpsilonExp = calculateEpsilon(ap_outputEntropy, ap_nw, precision);

    // If we get here, then adequate precision was used
    // Extract a value for display.
    // Note, this may round up, but we'll deal with this later.
    value =  mpfr_get_ld(ap_outputEntropy, MPFR_RNDN);
    mpfr_clears(ap_h_in, ap_entexp, ap_p_high, ap_p_low, ap_denom, ap_inputSpaceSize, ap_diff, ap_power_term, ap_psi, ap_omega, ap_outputEntropy, ap_nw, ap_n_out, NULL);

    return value;
}

//This function performs statistical testing on the bitwise input data.
//It returns h', which is a value in the range [0,1].
static long double computeEntropyOfConditionedData(string inputfilename, bool iid, TestRunBase *testRun) {

    data_t data;
    double h_bitstring = 1.0;
    const int verbose = 0;


    //Have the code establish the symbol width.
    data.word_size = 0;

    //Read in the complete file.
    if(!read_file_subset(inputfilename.c_str(), &data, 0, 0, testRun)) {
      fprintf(stderr, "Can't read and trandlate supplied data.\n");
      exit(-1);
    }

    if (iid) {
        // IID path
        //All of these run the bitstring version of the test (as per SP 800-90B Section 3.1.5.2 Paragraph 2)
        // Section 6.3.1 - Estimate entropy with Most Common Value
        h_bitstring = min(h_bitstring, most_common(data.bsymbols, data.blen, 2, verbose, "Bitstring"));
    } else {
        // NON-IID path
        double ret_min_entropy;
        double bin_t_tuple_res = -1.0, bin_lrs_res = -1.0;

        //All of these run the bitstring version of the test (as per SP 800-90B Section 3.1.5.2 Paragraph 2)
        // Section 6.3.1 - Estimate entropy with Most Common Value
        ret_min_entropy = most_common(data.bsymbols, data.blen, 2, verbose, "Bitstring");
        h_bitstring = min(ret_min_entropy, h_bitstring);

        // Section 6.3.2 - Estimate entropy with Collision Test
        ret_min_entropy = collision_test(data.bsymbols, data.blen, verbose, "Bitstring");
        h_bitstring = min(ret_min_entropy, h_bitstring);

        // Section 6.3.3 - Estimate entropy with Markov Test
        ret_min_entropy = markov_test(data.bsymbols, data.blen, verbose, "Bitstring");
        h_bitstring = min(ret_min_entropy, h_bitstring);

        // Section 6.3.4 - Estimate entropy with Compression Test
        ret_min_entropy = compression_test(data.bsymbols, data.blen, verbose, "Bitstring");
        if (ret_min_entropy >= 0) {
            h_bitstring = min(ret_min_entropy, h_bitstring);
        }

        //This call performs both the t-Tuple Test and the LRS Test
        SAalgs(data.bsymbols, data.blen, 2, bin_t_tuple_res, bin_lrs_res, verbose, "Bitstring");

        // Section 6.3.5 - Estimate entropy with t-Tuple Test
        if (bin_t_tuple_res >= 0.0) {
            h_bitstring = min(bin_t_tuple_res, h_bitstring);
        }

        // Section 6.3.6 - Estimate entropy with LRS Test
        if (bin_lrs_res >= 0) {
            h_bitstring = min(bin_lrs_res, h_bitstring);
        }

        // Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test
        ret_min_entropy = multi_mcw_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
        if (ret_min_entropy >= 0) {
            h_bitstring = min(ret_min_entropy, h_bitstring);
        }

        // Section 6.3.8 - Estimate entropy with Lag Prediction Test
        ret_min_entropy = lag_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
        if (ret_min_entropy >= 0) {
            h_bitstring = min(ret_min_entropy, h_bitstring);
        }

        // Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)
        ret_min_entropy = multi_mmc_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
        if (ret_min_entropy >= 0) {
            h_bitstring = min(ret_min_entropy, h_bitstring);
        }

        // Section 6.3.10 - Estimate entropy with LZ78Y Test
        ret_min_entropy = LZ78Y_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
        if (ret_min_entropy >= 0) {
            h_bitstring = min(ret_min_entropy, h_bitstring);
        }
    }

    free_data(&data);

    return h_bitstring;
}

int main(int argc, char* argv[]) {
    bool vetted, quietMode = false, iid = false;
    long double h_p = -1.0L;
    long double h_in, h_out;
    unsigned int n_in, n_out, nw, n;
    mpfr_prec_t precision;
    int opt;
    unsigned int maxval;

    long double noutEpsilonExp = -1.0L;
    long double hinEpsilonExp = -1.0L;
    long double nwEpsilonExp = -1.0L;
    long double outputEntropy = -1.0L;

    // Setting this rounding method helps prevent us from overestimating the input parameters
    fesetround(FE_TOWARDZERO);

    vetted = true;

    bool jsonOutput = false;
    string timestamp = getCurrentTimestamp();
    string outputfilename;
    string inputfilename;
    char *file_path;
    string commandline = recreateCommandLine(argc, argv);
    
    for (int i = 0; i < argc; i++) {
        std::string Str = std::string(argv[i]);
        if ("--version" == Str) {
            printVersion("conditioning");
            exit(0);
        }
    }

    while ((opt = getopt(argc, argv, "vnqo:i:c:")) != -1) {
        switch (opt) {
            case 'v':
                vetted = true;
                break;
            case 'n':
                vetted = false;
                break;
            case 'q':
                quietMode = true;
                break;
            case 'o':
                jsonOutput = true;
                outputfilename = optarg;
                break;
            case 'i':
                inputfilename = optarg;
                file_path = optarg;
                break;
            case 'c':
                iid = (strcmp(optarg, "iid") == 0);
                break;
            default:
                print_usage();
        }
    }

    argc -= optind;
    argv += optind;

    // Set up results
    // TODO individual structure for conditioning, do not reuse nonIID structure
    NonIidTestRun testRunNonIid;
    testRunNonIid.type = "Conditioning";
    testRunNonIid.timestamp = timestamp;
    testRunNonIid.commandline = commandline;
    
    if(!inputfilename.empty()) {
        testRunNonIid.filename = inputfilename;
    
        // Record hash of input file
        char hash[2*SHA256_DIGEST_LENGTH+1];

        sha256_file(file_path, hash);
        testRunNonIid.sha256 = hash;
    }
    
    // Parse args
    if (vetted && argc != 4) {
        printf("Incorrect usage.\n");
        print_usage();
    } else if (!vetted && (argc != 4 && argc != 5)) {
        printf("Incorrect usage.\n");
        printf("argc: %d\n", argc);
        print_usage();
    } else {

        // get n_in     
        n_in = inputUnsignedOption(argv[0], 1, UINT_MAX, "n_in");

        // get n_out     
        n_out = inputUnsignedOption(argv[1], 1, UINT_MAX, "n_out");

        // get nw     
        nw = inputUnsignedOption(argv[2], 1, UINT_MAX, "nw");

        // get h_in; note that h_in <= n_in
        h_in = inputLongDoubleOption(argv[3], 0.0L, (long double) n_in, "h_in");
        if (h_in <= 0.0L) {
            if (jsonOutput) {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "Error with input: generating h_in.";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
            print_usage();
        }

        if (!vetted) {
            if (argc == 4) {
                // If h_p is not provided but a file is provided instead, compute the entropy of that
                h_p = computeEntropyOfConditionedData(inputfilename, iid, &testRunNonIid);
            } else {
                // If h_p is provided via command line, use that value
                h_p = inputLongDoubleOption(argv[4], 0.0L, 1.0L, "h_p");
                if (h_p <= 0) {
                    if (jsonOutput) {
                        testRunNonIid.errorLevel = -1;
                        testRunNonIid.errorMsg = "Error with input: generating h_p.";
                        ofstream output;
                        output.open(outputfilename);
                        output << testRunNonIid.GetAsJson();
                        output.close();
                    }
                    print_usage();
                }
            }
        }
    }

    // Step 2 is invariant, and not subject to precision problems.
    nw = std::min(nw, n_in); // By 90B Appendix E
    n = std::min(n_out, nw);

    // Print out the inputs
    if (!quietMode) {
        printf("n_in = %u\n", n_in);
        printf("n_out = %u\n", n_out);
        printf("nw = %u\n", nw);
        printf("h_in = %.22Lg\n", h_in);
        if (!vetted) printf("h' = %.22Lg\n", h_p);
    }

    // Establish the maximum precision that ought to be necessary
    // If something goes wrong, we can increase this precision automatically.
    maxval = 53; // Always be large enough to faithfully represent h_in.
    maxval = (maxval > n_in) ? maxval : n_in;
    maxval = (maxval > n_out) ? maxval : n_out;
    maxval = (maxval > nw) ? maxval : nw;
    precision = 2 * maxval;

    // Check to see if this environment is going to support the needed exponent range
    assert(mpfr_get_emax() > maxval);
    assert(mpfr_get_emin() < -maxval);

    // Compute entropy
    outputEntropy = computeEntropyWithPrecision(precision, h_in, n_in, n, n_out, nw, noutEpsilonExp, hinEpsilonExp, nwEpsilonExp);

    // Check some basic bounds.
    assert(outputEntropy <= (long double) n_out);
    assert(outputEntropy <= h_in);
    assert(outputEntropy <= (long double) nw);
    assert(outputEntropy >= 0.0L);

    // We're done with the calculation. Now print results.
    if (!quietMode) {

        printf("Output_Entropy(*) = %.22Lg", outputEntropy);

        if (outputEntropy == (long double) n_out) {
            // outputEntropy rounded to full entropy, so the difference between this and full entropy is less than 1/2 ULP.
            printf("; Close to n_out (epsilon = 2^(-%.22Lg))", noutEpsilonExp);
        }
        
        if (outputEntropy == h_in) {
            // outputEntropy rounded to the input entropy, so the difference between this and the input entropy is less than 1/2 ULP.
            printf("; Close to h_in (epsilon = 2^(-%.22Lg))", hinEpsilonExp);
        }

        if (outputEntropy == (long double) nw) {
            // outputEntropy rounded to the nw, so the difference between this and nw is less than 1/2 ULP.
            printf("; Close to nw (epsilon = 2^(-%.22Lg))", nwEpsilonExp);
        }

        printf("\n");
    }

    if (vetted) {

        if (!quietMode)
            printf("(Vetted) h_out = %.22Lg\n", outputEntropy);

        h_out = outputEntropy;

        if (outputEntropy > 0.999L * ((long double) n_out)) {

            //h_out = (1 - epsilon) * n_out
            if (!quietMode) {
                printf("epsilon = 2^(-%.22Lg)", noutEpsilonExp);

                //Should this qualify as "full entropy" under FIPS 140-3 IG D.K Resolution 19
		if (h_in >= n_out + 64.0) {
                    printf(": FIPS 140-3 IG D.K Resolution 19 Full Entropy if the conditioning component security strength is >= %u", n_out);
		}
                printf("\n");
            }
        }
    } else {
        long double bound90B = 0.999L * ((long double) n_out);
        long double statBound = h_p * ((long double) n_out);

        //Note, we can't assess as full entropy in this case.
        if (!quietMode) {
            printf("0.999 * n_out = %.22Lg\n", bound90B);
            printf("h' * n_out = %.22Lg\n", statBound);
        }

        h_out = std::min(outputEntropy, std::min(bound90B, statBound));

        if (!quietMode)
            printf("(Non-vetted) h_out = %.22Lg\n", h_out);
    }

    NonIidTestCase tcOverallnonIid;

    tcOverallnonIid.testCaseNumber = "Overall";
    //tcOverallnonIid.vetted = vetted;
    tcOverallnonIid.n_in = n_in;
    tcOverallnonIid.n_out = n_out;
    tcOverallnonIid.nw = nw;
    tcOverallnonIid.h_in = h_in;
    tcOverallnonIid.h_out = h_out;
    tcOverallnonIid.h_p = h_p;

    testRunNonIid.testCases.push_back(tcOverallnonIid);
    testRunNonIid.errorLevel = 0;

    if (jsonOutput) {
        ofstream output;
        output.open(outputfilename);
        output << testRunNonIid.GetAsJson();
        output.close();
    }

    return 0;
}
