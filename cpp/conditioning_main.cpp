#include "other/conditioning_test_run.h"
#include "non_iid/collision_test.h"
#include "non_iid/lz78y_test.h"
#include "non_iid/multi_mmc_test.h"
#include "non_iid/lag_test.h"
#include "non_iid/multi_mcw_test.h"
#include "non_iid/compression_test.h"
#include "non_iid/markov_test.h"
#include "iid/iid_test_run.h"
#include "shared/TestRunUtils.h"
#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include <string.h>
#include <stdio.h>
#include <cstdlib>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <limits.h>


#include <getopt.h>

#include <iostream>
#include <fstream>

[[ noreturn ]] void print_usage() {
    printf("Usage is: ea_conditioning [-v] [-q] <n_in> <n_out> <nw> <h_in> [-c 'iid' | 'non_iid'] [-i inputfile_.json] [-o output_file.json]\n");
    printf("\tor \n\tea_conditioning -n <n_in> <n_out> <nw> <h_in> <h'> [-i inputfile_.json] [-w word size] [-o output_file.json]\n\n");
    printf("\t <n_in>: input number of bits to conditioning function.\n");
    printf("\t <n_out>: output number of bits from conditioning function.\n");
    printf("\t <nw>: narrowest internal width of conditioning function.\n");
    printf("\t <h_in>: input entropy to conditioning function.\n");
    printf("\t <-v|-n>: '-v' for vetted conditioning function, '-n' for non-vetted conditioning function. Vetted conditioning is the default.\n");
    printf("\t <h'>: entropy estimate per bit of conditioned sequential dataset (only for '-n' option).\n");
    printf("\t -q: Quiet mode, less output to screen.\n");
    printf("\t -i inputfile.json: select inputfile\n");
    printf("\t -c ['iid'|'non-iid']: select iid or non-iid\n");
    printf("\n");
    printf("\t -o: Set Output Type to JSON\n");
    printf("\n");
    printf("\t\t Changes the output format to JSON and sets the file location for the output file.\n");
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
    bool quietMode = false;
    double h_p = -1.0;
    double h_in, h_out, n_in, n_out, nw, n, p_high, p_low, psi, omega, output_entropy, power_term;
    char opt;

    vetted = true;

    bool jsonOutput = false;
    string timestamp = getCurrentTimestamp();
    string outputfilename;
    string inputfilename = "";
    bool iid = false;

    // Parse args
    // if ((argc < 4)) {
    //     printf("Incorrect usage.\n");
    //     print_usage();
    // }

    while ((opt = getopt(argc, argv, "c:vnqi:o:")) != -1) {
   
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
                break;
            case 'c':
                if(!optarg) {
                    print_usage();
                    break;
                }
                if (strcmp(optarg, "iid") == 0)
                    iid = true;
                else
                    iid = false;
                break;
            default:
                print_usage();
        }
    }

    argc -= optind;
    argv += optind;

    ConditioningTestRun testRun;
    testRun.type = "Conditioning";
    testRun.timestamp = timestamp;
    testRun.IID = -1;

    // get n_in 
    n_in = atof(argv[0]);
    if ((n_in <= 0) || (floor(n_in) != n_in)) {
        
        printf("n_in must be a positive integer.\n");
        if (jsonOutput) {
            testRun.errorLevel = -1;
            testRun.errorMsg = "n_in must be a positive integer";
            ofstream output;
            output.open(outputfilename);
            output << testRun.GetAsJson();
            output.close();
        }
        print_usage();
    }

    // get n_out
    n_out = atof(argv[1]);
    if ((n_out <= 0) || (floor(n_out) != n_out)) {
        printf("n_out must be a positive integer.\n");

        if (jsonOutput) {
            testRun.errorLevel = -1;
            testRun.errorMsg = "n_out must be a positive integer";
            ofstream output;
            output.open(outputfilename);
            output << testRun.GetAsJson();
            output.close();
        }

        print_usage();
    }

    // get n_w
    nw = atof(argv[2]);
    if ((nw <= 0) || (floor(nw) != nw)) {
        printf("n_w must be a positive integer.\n");
        if (jsonOutput) {
            testRun.errorLevel = -1;
            testRun.errorMsg = "n_w must be a positive integer";
            ofstream output;
            output.open(outputfilename);
            output << testRun.GetAsJson();
            output.close();
        }
        print_usage();
    }

    // get h_in 
    h_in = atof(argv[3]);
    if ((h_in <= 0) || (h_in > n_in)) {
        printf("h_in must be positive and at most n_in.\n");
        if (jsonOutput) {
            testRun.errorLevel = -1;
            testRun.errorMsg = "h_in must be a positive and at most n_in.";
            ofstream output;
            output.open(outputfilename);
            output << testRun.GetAsJson();
            output.close();
        }
        print_usage();
    }

    if (!vetted) {
        if (argc == 5) {
            // Reading h_p from commandline
            h_p = atof(argv[4]);
            if ((h_p < 0) || (h_p > 1.0)) {
                printf("h' must be between 0 and 1 inclusive.\n");
                if (jsonOutput) {
                    testRun.errorLevel = -1;
                    testRun.errorMsg = "h' must be between 0 and 1 inclusive.";
                    ofstream output;
                    output.open(outputfilename);
                    output << testRun.GetAsJson();
                    output.close();
                }
                print_usage();
            }
        } else {
            // We calculate h_p based on input file

            char hash[65];
            sha256_file(&inputfilename[0], hash);
            testRun.sha256 = hash;
            testRun.filename = inputfilename.c_str();

            data_t data;
            data.word_size = 0;

            double H_bitstring = 1.0;
            bool all_bits = true; // this is the default setting in ea_iid & ea_non_iid
            bool initial_entropy = false; // mimics the -c flag
            bool verbose = false;
            unsigned long subsetIndex = ULONG_MAX;
            unsigned long subsetSize = 0;

            read_file_subset(inputfilename.c_str(), &data, subsetIndex, subsetSize);
            double h_assessed = data.word_size;
            double H_original = data.word_size;
            int sample_size = data.len;
            int alphabet_size = data.alph_size;
            
            if (iid) {
                // IID path
                if (!all_bits && (data.blen > MIN_SIZE)) data.blen = MIN_SIZE;
                if (initial_entropy) {
                    H_original = most_common(data.symbols, sample_size, alphabet_size, verbose, "Literal");
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    H_bitstring = most_common(data.bsymbols, data.blen, 2, verbose, "Bitstring");
                }

                if ((data.alph_size > 2)) {
                    h_assessed = min(h_assessed, H_bitstring * data.word_size);
                }
                h_assessed = min(h_assessed, H_bitstring * data.word_size);


                H_original = most_common(data.symbols, sample_size, alphabet_size, false, "Literal");
                h_assessed = min(h_assessed, H_original);

                if (h_assessed != 0) {
                    h_assessed = h_assessed / data.word_size;
                }
                h_p = h_assessed;
            } else {
                // NON-IID path
                h_assessed = data.word_size;
                double ret_min_entropy = 0.0;
                double H_original = data.word_size;
                double bin_t_tuple_res = -1.0, bin_lrs_res = -1.0;
                double t_tuple_res = -1.0, lrs_res = -1.0;

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = most_common(data.bsymbols, data.blen, 2, verbose, "Bitstring");
                    H_bitstring = min(ret_min_entropy, H_bitstring);
                }

                if (initial_entropy) {
                    ret_min_entropy = most_common(data.symbols, data.len, data.alph_size, verbose, "Literal");
                    H_original = min(ret_min_entropy, H_original);
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = collision_test(data.bsymbols, data.blen, verbose, "Bitstring");

                    H_bitstring = min(ret_min_entropy, H_bitstring);
                }

                if (initial_entropy && (data.alph_size == 2)) {
                    ret_min_entropy = collision_test(data.symbols, data.len, verbose, "Literal");
                    H_original = min(ret_min_entropy, H_original);
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = markov_test(data.bsymbols, data.blen, verbose, "Bitstring");

                    H_bitstring = min(ret_min_entropy, H_bitstring);
                }

                if (initial_entropy && (data.alph_size == 2)) {
                    ret_min_entropy = markov_test(data.symbols, data.len, verbose, "Literal");
                    H_original = min(ret_min_entropy, H_original);
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = compression_test(data.bsymbols, data.blen, verbose, "Bitstring");
                    if (ret_min_entropy >= 0) {
                        H_bitstring = min(ret_min_entropy, H_bitstring);
                    }
                }

                if (initial_entropy && (data.alph_size == 2)) {
                    ret_min_entropy = compression_test(data.symbols, data.len, verbose, "Literal");
                    H_original = min(ret_min_entropy, H_original);
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    SAalgs(data.bsymbols, data.blen, 2, bin_t_tuple_res, bin_lrs_res, verbose, "Bitstring");
                    if (bin_t_tuple_res >= 0.0) {

                        H_bitstring = min(bin_t_tuple_res, H_bitstring);
                    }
                }

                if (initial_entropy) {
                    SAalgs(data.symbols, data.len, data.alph_size, t_tuple_res, lrs_res, verbose, "Literal");
                    if (t_tuple_res >= 0.0) {

                        H_original = min(t_tuple_res, H_original);
                    }
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    H_bitstring = min(bin_lrs_res, H_bitstring);
                }

                if (initial_entropy) {
                    H_original = min(lrs_res, H_original);
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = multi_mcw_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
                    if (ret_min_entropy >= 0) {

                        H_bitstring = min(ret_min_entropy, H_bitstring);
                    }
                }

                if (initial_entropy) {
                    ret_min_entropy = multi_mcw_test(data.symbols, data.len, data.alph_size, verbose, "Literal");
                    if (ret_min_entropy >= 0) {

                        H_original = min(ret_min_entropy, H_original);
                    }
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = lag_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
                    if (ret_min_entropy >= 0) {
                        if (verbose == 1) printf("\tLag Prediction Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
                        H_bitstring = min(ret_min_entropy, H_bitstring);
                    }
                }

                if (initial_entropy) {
                    ret_min_entropy = lag_test(data.symbols, data.len, data.alph_size, verbose, "Literal");
                    if (ret_min_entropy >= 0) {
                        if (verbose == 1) printf("\tLag Prediction Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
                        H_original = min(ret_min_entropy, H_original);
                    }
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = multi_mmc_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
                    if (ret_min_entropy >= 0) {

                        H_bitstring = min(ret_min_entropy, H_bitstring);
                    }
                }

                if (initial_entropy) {
                    ret_min_entropy = multi_mmc_test(data.symbols, data.len, data.alph_size, verbose, "Literal");
                    if (ret_min_entropy >= 0) {

                        H_original = min(ret_min_entropy, H_original);
                    }
                }

                if (((data.alph_size > 2) || !initial_entropy)) {
                    ret_min_entropy = LZ78Y_test(data.bsymbols, data.blen, 2, verbose, "Bitstring");
                    if (ret_min_entropy >= 0) {
                        if (verbose == 1) printf("\tLZ78Y Prediction Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
                        H_bitstring = min(ret_min_entropy, H_bitstring);
                    }
                }

                if (initial_entropy) {
                    ret_min_entropy = LZ78Y_test(data.symbols, data.len, data.alph_size, verbose, "Literal");
                    if (ret_min_entropy >= 0) {
                        if (verbose == 1) printf("\tLZ78Y Prediction Test Estimate = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
                        H_original = min(ret_min_entropy, H_original);
                    }
                }

                h_assessed = data.word_size;

                if ((data.alph_size > 2) || !initial_entropy) {
                    h_assessed = min(h_assessed, H_bitstring * data.word_size);
                    printf("H_bitstring = %.17g\n", H_bitstring);
                }

                if (initial_entropy) {
                    h_assessed = min(h_assessed, H_original);
                    printf("H_original: %.17g\n", H_original);
                }

                printf("Assessed min entropy: %g\n", h_assessed);
                if (h_assessed != 0) {
                    h_assessed = h_assessed / data.word_size;
                }
                h_p = h_assessed;
            }
            
        }
    }

    if (nw > n_in) nw = n_in;

    if (!quietMode) {
        printf("n_in: %f\n", n_in);
        printf("n_out: %f\n", n_out);
        printf("nw: %f\n", nw);
        printf("h_in: %f\n", h_in);
        if (!vetted) printf("h': %f\n", h_p);
    }

    // compute Output Entropy (Section 3.1.5.1.2)
    p_high = pow(2.0, -h_in);
    p_low = (1.0 - p_high) / (pow(2.0, n_in) - 1);
    n = std::min(n_out, nw);
    power_term = pow(2.0, n_in - n);
    psi = power_term * p_low + p_high;
    omega = (power_term + sqrt(2 * n * power_term * log(2))) * p_low;
    output_entropy = -log2(std::max(psi, omega));

    if (output_entropy > n_out) output_entropy = n_out;
    if (output_entropy < 0) output_entropy = 0;

    if (!quietMode) {
        if (vetted) {
            printf("\n(Vetted) ");
            h_out = output_entropy;
        } else {
            printf("\n(Non-vetted) ");
            h_out = std::min(output_entropy, std::min(0.999 * n_out, h_p * n_out));
        }

        printf("h_out: %f\n", h_out);
    }

    ConditioningTestCase tcOverall;

    tcOverall.testCaseNumber = "Overall";
    tcOverall.vetted = vetted;
    tcOverall.n_in = n_in;
    tcOverall.n_out = n_out;
    tcOverall.nw = nw;
    tcOverall.h_in = h_in;
    tcOverall.h_out = h_out;
    tcOverall.h_p = h_p;

    testRun.testCases.push_back(tcOverall);
    testRun.errorLevel = 0;

    if (jsonOutput) {
        ofstream output;
        output.open(outputfilename);
        output << testRun.GetAsJson();
        output.close();
    }

    return 0;
}
