/* VERSION information is kept in utils.h. Please update when a new version is released */

#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "non_iid/non_iid_test_run.h"
#include "iid/iid_test_run.h"
#include "shared/TestRunUtils.h"
#include "non_iid/collision_test.h"
#include "non_iid/lz78y_test.h"
#include "non_iid/multi_mmc_test.h"
#include "non_iid/lag_test.h"
#include "non_iid/multi_mcw_test.h"
#include "non_iid/compression_test.h"
#include "non_iid/markov_test.h"
#include "iid/chi_square_tests.h"
#include "iid/permutation_tests.h"

#include <cstdint>
#include <getopt.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <openssl/sha.h>

//Each test has a targeted chance of roughly 0.000005, and we need to witness at least 5 failures, so this should be no less than 1000000
#define DEFAULT_SIMULATION_ROUNDS 5000000UL

[[ noreturn ]] void print_usage() {
    printf("Usage is: ea_restart [-i|-n] [-v] [-q] [-s <simulation count>] <file_name> [bits_per_symbol] <H_I>\n\n");
    printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (samples),\n");
    printf("\t and in the \"row dataset\" format described in SP800-90B Section 3.1.4.1.\n");
    printf("\t [bits_per_symbol]: Must be between 1-8, inclusive.\n");
    printf("\t <H_I>: Initial entropy estimate.\n");
    printf("\t [-i|-n]: '-i' for IID data, '-n' for non-IID data. Non-IID is the default.\n");
    printf("\t -s <simulation count>: Establish cutoff using <simulation count> rounds.\n");
    printf("\t -v: Optional verbosity flag for more output.\n");
    printf("\t -q: Quiet mode, less output to screen.\n");
    printf("\n");
    printf("\t Restart samples are assumed to be packed into 8-bit values, where the rightmost 'bits_per_symbol'\n");
    printf("\t bits constitute the sample.\n");
    printf("\n");
    printf("\t This program performs restart testing as described in Restart Tests (Section 3.1.4). The data\n");
    printf("\t consists of 1000 restarts, each with 1000 samples. The data is converted to rows and columns\n");
    printf("\t as described Section 3.1.4.1. The sanity check (Section 3.1.4.3) and the validation test\n");
    printf("\t (Section 3.1.4.2) are performed on this data.\n");
    printf("\n");
    printf("\t If the restart data passes the sanity check and validation test, this program returns\n");
    printf("\t min(H_r, H_c, H_I), which is either the validated entropy assessment or used to derive\n");
    printf("\t 'h_in' if conditioning is used (Section 3.1.5).\n");
    printf("\n");
    printf("\t --version: Prints tool version information");
    printf("\n");
    exit(-1);
}

// Here, we simulate a "worst case" for the restart sanity test. This is "worst case" in the sense that the adopted distribution
// results in the largest acceptable collision bound for a given assessed entropy level, so if a data sample fails this
// test, it is likely to indicate an underlying problem.
//
// This "worst case" uses the "inverted near-uniform" family (see Hagerty-Draper "Entropy Bounds and Statistical Tests" for
// a full definition of this distribution and justification for its use here).
//
// This distribution has as many maximal probability symbols as possible (each occurring with probability p), and possibly one
// additional symbol that contains all the residual probability.
//
// If the probability for the most likely symbol is p, then there are floor(1/p) most likely symbols,
// each occurring with probability p and possibly one additional symbol that has all the remaining (1 - p floor(1/p)) chance.
// In this code, we generate a random unit value in the range [0, 1), and we need to map this to one of the ceil(1/p) possible
// output symbols.
//
// Note that the function x -> floor(x/p) yields
// [0p,1p) -> 0
// [1p, 2p) -> 1
// [2p, 3p) -> 2
// ...
// [(floor(1/p)-1)p, floor(1/p)p) -> floor(1/p)-1
// [ floor(1/p)p, 1 ) -> floor(1/p)
//
// As such, each of the first floor(1/p) symbols (0 through floor(1/p)-1) have probability p of occurring, and
// the symbol floor(1/p) has probability 1-floor(1/p)p of occurring, as desired.
//
// Note that if floor(1/p) = ceil(1/p) = 1/p, then there is no "residual" symbol, only 1/p most likely symbols.
//
// The array is 0-indexed, so we can use this map to establish the index directly.
uint16_t simulateCount(int k_effective, double p, uint64_t *xoshiro256starstarState) {
    uint16_t counts[256] = {0};
    uint16_t max_count = 0;

    for (int j = 0; j < 1000; j++) {
        // Note that (int)floor(randomUnit(xoshiro256starstarState) / p) is the index map discussed in the above comments.
        counts[(int)floor(randomUnit(xoshiro256starstarState) / p)]++;
    }

    // We could have tracked this during the above loop, but that would yield 1000 comparisons,
    // rather than k_effective (<= 256) comparisons, as here.
    for (int j = 0; j < k_effective; j++) {
        if (max_count < counts[j]) max_count = counts[j];
    }

    return max_count;
}

//This returns the bound (cutoff) for the test. Counts equal to this value should pass.
//Larger values should fail.

int simulateBound(double alpha, int k, double H_I, unsigned long int simulation_rounds) {
    uint64_t xoshiro256starstarMainSeed[4];
    uint16_t *results;
    long int returnIndex;
    double p;
    int k_effective;
    int returnValue;

    assert((k > 1) && (k <= 256));

    // A few constraints:
    // This array may be very large (many gigabytes) so can't go onto the stack
    // Many of the C++ STL-derived types are not thread safe. Our mode of access is
    // quite straight forward, but there are clearly issues in some cases.
    // In C, calloc is possibly faster, but this is probably the best we can do
    // using somewhat idiomatic C++.
    results = new uint16_t[simulation_rounds];
    memset(results, 0, sizeof(uint16_t)*simulation_rounds);

    //The probability of the most likely symbol (MLS) only needs to be calculated once...
    p = pow(2.0, -H_I);

    //if floor(1/p) = ceil(1/p) = 1/p, then there are exactly that many symbols (e.g., p=1/2, then there are 2 symbols expected).
    //If ceil(1/p) > 1/p, then ceil(1/p) = floor(1/p)+1, that is there are the floor(1/p) most likely symbols, and then the extra
    //symbol that fills the rest of the space (with probability < p).
    k_effective = ceil(1.0 / p);
    assert(k_effective <= k);

    seed(xoshiro256starstarMainSeed);

#pragma omp parallel
    {
        uint64_t xoshiro256starstarSeed[4];

        memcpy(xoshiro256starstarSeed, xoshiro256starstarMainSeed, sizeof (xoshiro256starstarMainSeed));
        //Cause the RNG to jump omp_get_thread_num() * 2^128 calls
        xoshiro_jump(omp_get_thread_num(), xoshiro256starstarSeed);

#pragma omp for
        for (unsigned long int i = 0; i < simulation_rounds; i++) {
            results[i] = simulateCount(k_effective, p, xoshiro256starstarSeed);
        }
    }

    sort(results, results+simulation_rounds);
    assert((results[0] >= (1000 / k)) && (results[0] <= 1000));
    assert((results[simulation_rounds - 1] >= (1000 / k)) && (results[simulation_rounds - 1] <= 1000));

    returnIndex = ((size_t) floor((1.0 - alpha) * ((double) simulation_rounds))) - 1;
    assert((returnIndex >= 0) && (returnIndex < simulation_rounds));

    returnValue = (int)results[returnIndex];

    delete results;

    return returnValue;
}

int main(int argc, char* argv[]) {
    bool iid;
    int verbose = 1; //verbose 0 is for JSON output, 1 is the normal mode, 2 is the NIST tool verbose mode, and 3 is for extra verbose output
    bool quietMode = false;
    char *file_path;
    int r = 1000, c = 1000;
    int counts[256];
    unsigned long int simulation_rounds = DEFAULT_SIMULATION_ROUNDS;
    int X_cutoff;
    int i, j;
    int X_i, X_r, X_c, X_max;
    double H_I, H_r, H_c, alpha, ret_min_entropy;
	double rawmean, median;
    uint8_t *rdata, *cdata;
    unsigned long int inul;
    data_t data;
    int opt;

    iid = false;
    data.word_size = 0;

    bool jsonOutput = false;
    string timestamp = getCurrentTimestamp();
    string outputfilename = timestamp + ".json";
    string commandline = recreateCommandLine(argc, argv);

    for (int i = 0; i < argc; i++) {
        std::string Str = std::string(argv[i]);
        if ("--version" == Str) {
            printVersion("restart");
            exit(0);
        }
    }

    while ((opt = getopt(argc, argv, "invqo:s:")) != -1) {
        switch (opt) {
            case 'i':
                iid = true;
                break;
            case 'n':
                iid = false;
                break;
            case 'v':
                verbose++;
                break;
            case 'q':
                quietMode = true;
                break;
            case 'o':
                jsonOutput = true;
                outputfilename = optarg;
                break;
            case 's':
                inul = strtoul(optarg, NULL, 10);
		if((inul == 0) || (inul == ULONG_MAX) || (inul < simulation_rounds)) {
                    print_usage();
                } else {
                    simulation_rounds = inul;
                }
                break;
            default:
                print_usage();
        }
    }

    argc -= optind;
    argv += optind;

    // Parse args
    if ((argc != 3) && (argc != 2)) {
        printf("Incorrect usage.\n");
        print_usage();
    }

    // get filename
    file_path = argv[0];
    argv++;
    argc--;

    if (quietMode) verbose = 0;

    char hash[2*SHA256_DIGEST_LENGTH+1];
    sha256_file(file_path, hash);

    IidTestRun testRunIid;
    testRunIid.type = "Restart";
    testRunIid.timestamp = timestamp;
    testRunIid.filename = file_path;
    testRunIid.commandline = commandline;
    testRunIid.sha256 = hash;

    NonIidTestRun testRunNonIid;
    testRunNonIid.type = "Restart";
    testRunNonIid.timestamp = timestamp;
    testRunNonIid.sha256 = hash;
    testRunNonIid.filename = file_path;
    testRunNonIid.commandline = commandline;

    if (argc == 2) {
        // get bits per word
        data.word_size = atoi(argv[0]);
        if (data.word_size < 1 || data.word_size > 8) {
            printf("Invalid bits per symbol.\n");
            if (jsonOutput) {
                if (iid) {
                    testRunIid.errorLevel = -1;
                    testRunIid.errorMsg = "Invalid bits per symbol.";
                    ofstream output;
                    output.open(outputfilename);
                    output << testRunIid.GetAsJson();
                    output.close();
                } else {
                    testRunNonIid.errorLevel = -1;
                    testRunNonIid.errorMsg = "Invalid bits per symbol.";
                    ofstream output;
                    output.open(outputfilename);
                    output << testRunNonIid.GetAsJson();
                    output.close();
                }
            }

            print_usage();
        }
        argv++;
        argc--;
    }

    // get H_I	
    H_I = atof(argv[0]);
    if (H_I < 0) {
        printf("H_I %f must be nonnegative.\n", H_I);
        if (jsonOutput) {
            if (iid) {
                testRunIid.errorLevel = -1;
                testRunIid.errorMsg = "H_I must be nonnegative: " + std::to_string(H_I) + ".";
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "H_I must be nonnegative: " + std::to_string(H_I) + ".";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }

        print_usage();
    }

    if (verbose > 1) printf("Opening file: '%s' (SHA-256 hash %s)\n", file_path, hash);

    if (!read_file(file_path, &data, &testRunNonIid)) {
        printf("Error reading file.\n");

        if (jsonOutput) {
            if (iid) {
                testRunNonIid.errorLevel = -1;
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }

        print_usage();
    }

    if (verbose > 1) printf("Loaded %ld samples made up of %d distinct %d-bit-wide symbols.\n", data.len, data.alph_size, data.word_size);

    if (H_I > data.word_size) {
        printf("H_I (%f) must be at most 'bits_per_symbol' (%d).\n", H_I, data.word_size);
        if (jsonOutput) {
            if(iid) {
                testRunIid.errorLevel = -1;
                testRunIid.errorMsg = "H_I (" + std::to_string(H_I) + ") must be at most 'bits_per_symbol' (" + std::to_string(data.word_size) + ").";
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "H_I (" + std::to_string(H_I) + ") must be at most 'bits_per_symbol' (" + std::to_string(data.word_size) + ").";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }
        free_data(&data);
        exit(-1);
    }

    if (data.alph_size <= 1) {
        printf("Symbol alphabet consists of 1 symbol. No entropy awarded...\n");
        if (jsonOutput) {
            if(iid) {
                testRunIid.errorLevel = -1;
                testRunIid.errorMsg = "Symbol alphabet consists of 1 symbol. No entropy awarded...";
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "Symbol alphabet consists of 1 symbol. No entropy awarded...";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }
        free_data(&data);
        exit(-1);
    }

    if (data.len != MIN_SIZE) {
        printf("\n*** Error: data (len = %ld) does not contain %d samples ***\n\n", data.len, MIN_SIZE);
        if (jsonOutput) {
            if (iid) {
                testRunIid.errorLevel = -1;
                testRunIid.errorMsg = "*** Error: data (len = " + std::to_string(data.len) + ") does not contain " + std::to_string(MIN_SIZE) + " samples ***";
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "*** Error: data (len = " + std::to_string(data.len) + ") does not contain " + std::to_string(MIN_SIZE) + " samples ***";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }
        exit(-1);
    }

    if (verbose > 1) {
        if (data.alph_size < (1 << data.word_size)) printf("\nSymbols have been translated.\n\n");
    }

    rdata = data.symbols;
    cdata = (uint8_t*) malloc(data.len);
    if (cdata == NULL) {
        printf("Error: failure to initialize memory for columns\n");
        if (jsonOutput) {
            if(iid) {
                testRunIid.errorLevel = -1;
                testRunIid.errorMsg = "Error: failure to initialize memory for columns";
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "Error: failure to initialize memory for columns";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }
        exit(-1);
    }

    printf("H_I: %f\n", H_I);

    alpha = 1 - exp(log(0.99) / (r + c));
    X_cutoff = simulateBound(alpha, data.alph_size, H_I, simulation_rounds);
    if (verbose > 0) printf("ALPHA: %.17g, X_cutoff: %d\n", alpha, X_cutoff);

    // get maximum row count
    X_r = 0;
    for (i = 0; i < r; i++) { //row
        memset(counts, 0, 256 * sizeof (int));
        X_i = 0;
        for (j = 0; j < c; j++) {//column
            //[i*r+j] is row i, column j
            //So, we're fixing a row, and then iterate through various columns
            if (++counts[rdata[i * r + j]] > X_i) X_i = counts[rdata[i * r + j]];
        }
        if (X_i > X_r) X_r = X_i;
    }

    // construct column data from row data and get maximum column count
    X_c = 0;
    for (j = 0; j < c; j++) { //columns
        memset(counts, 0, 256 * sizeof (int));
        X_i = 0;
        for (i = 0; i < r; i++) {
            //[i*r+j] is row i, column j
            //So, we're fixing a column and iterating through various rows
            cdata[j * c + i] = rdata[i * r + j];
            if (++counts[cdata[j * c + i]] > X_i) X_i = counts[cdata[j * c + i]];
        }
        if (X_i > X_c) X_c = X_i;
    }

    // perform sanity check on rows and columns of restart data (Section 3.1.4.3)
    X_max = max(X_r, X_c);
    if (verbose > 0) printf("X_max: %d\n", X_max);
    if (X_max > X_cutoff) {
        if (verbose > 0) printf("\n*** Restart Sanity Check Failed ***\n");
        if (jsonOutput) {
            if(iid) {
                testRunIid.errorLevel = -1;
                testRunIid.errorMsg = "Restart Sanity Check Failed.";
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "Restart Sanity Check Failed.";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }
        exit(-1);
    } else if (verbose > 1) printf("\nRestart Sanity Check Passed...\n");


    // Calculate baseline statistics
    int alphabet_size = data.alph_size;
    int sample_size = data.len;

    if ((verbose == 1) || (verbose == 2))
        printf("Calculating baseline statistics...\n");

    calc_stats(&data, rawmean, median);

    // The maximum min-entropy is -log2(1/2^word_size) = word_size
    H_c = data.word_size;
    H_r = data.word_size;
    if (verbose > 0) {
        if (iid) printf("\nRunning IID tests...\n\n");
        else printf("\nRunning non-IID tests...\n\n");

        printf("Running Most Common Value Estimate...\n");
    }

    // Section 6.3.1 - Estimate entropy with Most Common Value

    NonIidTestCase tc631nonIid;
    tc631nonIid.testCaseNumber = "Most Common Value";
    tc631nonIid.data_word_size = data.word_size;

    IidTestCase tc631Iid;
    tc631Iid.testCaseNumber = "Most Common Value";
    tc631Iid.data_word_size = data.word_size;

    ret_min_entropy = most_common(rdata, data.len, data.alph_size, verbose, "Literal");
    if (verbose > 1) printf("\tMost Common Value Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
    tc631nonIid.h_r = ret_min_entropy;
    tc631Iid.h_r = ret_min_entropy;

    H_r = min(ret_min_entropy, H_r);

    ret_min_entropy = most_common(cdata, data.len, data.alph_size, verbose, "Literal");
    if (verbose > 1) printf("\tMost Common Value Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
    tc631nonIid.h_c = ret_min_entropy;
    H_c = min(ret_min_entropy, H_c);

    testRunNonIid.testCases.push_back(tc631nonIid);

    IidTestCase tcOverallIid;
    tcOverallIid.h_r = H_r;
    tcOverallIid.h_c = H_c;
    tcOverallIid.h_i = H_I;
    tcOverallIid.testCaseNumber = "Overall";


    if (!iid) {

        if (data.alph_size == 2) {

            if (verbose > 0) printf("\nRunning Entropic Statistic Estimates (bit strings only)...\n");

            // Section 6.3.2 - Estimate entropy with Collision Test (for bit strings only)

            NonIidTestCase tc632;
            tc632.testCaseNumber = "Collision Test (for bit strings only)";
            tc632.data_word_size = 1;
            ret_min_entropy = collision_test(rdata, data.len, verbose, "Literal");
            if (verbose > 1) printf("\tCollision Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy);
            tc632.h_r = ret_min_entropy;
            H_r = min(ret_min_entropy, H_r);

            ret_min_entropy = collision_test(cdata, data.len, verbose, "Literal");
            if (verbose > 1) printf("\tCollision Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy);
            tc632.h_c = ret_min_entropy;
            H_c = min(ret_min_entropy, H_c);

            testRunNonIid.testCases.push_back(tc632);

            // Section 6.3.3 - Estimate entropy with Markov Test (for bit strings only)

            NonIidTestCase tc633;
            tc633.testCaseNumber = "Markov Test (for bit strings only)";
            tc633.data_word_size = 1;

            ret_min_entropy = markov_test(rdata, data.len, verbose, "Literal");
            if (verbose > 1) printf("\tMarkov Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy);
            tc633.h_r = ret_min_entropy;
            H_r = min(ret_min_entropy, H_r);

            ret_min_entropy = markov_test(cdata, data.len, verbose, "Literal");
            if (verbose > 1) printf("\tMarkov Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy);
            tc633.h_c = ret_min_entropy;
            H_c = min(ret_min_entropy, H_c);

            testRunNonIid.testCases.push_back(tc633);

            // Section 6.3.4 - Estimate entropy with Compression Test (for bit strings only)

            NonIidTestCase tc634;
            tc634.testCaseNumber = "Compression Test (for bit strings only)";
            tc634.data_word_size = 1;

            ret_min_entropy = compression_test(rdata, data.len, verbose, "Literal");
            if (ret_min_entropy >= 0) {
                if (verbose > 1) printf("\tCompression Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy);
                tc634.h_r = ret_min_entropy;
                H_r = min(ret_min_entropy, H_r);
            }

            ret_min_entropy = compression_test(cdata, data.len, verbose, "Literal");
            if (ret_min_entropy >= 0) {
                if (verbose > 1) printf("\tCompression Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy);
                tc634.h_c = ret_min_entropy;
                H_c = min(ret_min_entropy, H_c);
            }
            testRunNonIid.testCases.push_back(tc634);

        }


        if (verbose > 0) printf("\nRunning Tuple Estimates...\n");

        // Section 6.3.5 - Estimate entropy with t-Tuple Test

        NonIidTestCase tc635;
        tc635.testCaseNumber = "T-Tuple Test";
        tc635.data_word_size = data.word_size;

        double row_t_tuple_res, row_lrs_res;
        double col_t_tuple_res, col_lrs_res;
        SAalgs(rdata, data.len, data.alph_size, row_t_tuple_res, row_lrs_res, verbose, "Literal");
        SAalgs(cdata, data.len, data.alph_size, col_t_tuple_res, col_lrs_res, verbose, "Literal");

        if (verbose > 1) printf("\tT-Tuple Test Estimate (Rows) = %f / %d bit(s)\n", row_t_tuple_res, data.word_size);
        tc635.h_r = row_t_tuple_res;
        H_r = min(row_t_tuple_res, H_r);

        if (verbose > 1) printf("\tT-Tuple Test Estimate (Cols) = %f / %d bit(s)\n", col_t_tuple_res, data.word_size);
        tc635.h_c = col_t_tuple_res;
        H_c = min(col_t_tuple_res, H_c);

        testRunNonIid.testCases.push_back(tc635);
        // Section 6.3.6 - Estimate entropy with LRS Test

        NonIidTestCase tc636;
        tc636.testCaseNumber = "LRS Test";
        tc636.data_word_size = data.word_size;

        if (verbose > 1) printf("\tLRS Test Estimate (Rows) = %f / %d bit(s)\n", row_lrs_res, data.word_size);
        tc636.h_r = row_lrs_res;
        H_r = min(row_lrs_res, H_r);

        if (verbose > 1) printf("\tLRS Test Estimate (Cols) = %f / %d bit(s)\n", col_lrs_res, data.word_size);
        tc636.h_c = col_lrs_res;
        H_c = min(col_lrs_res, H_c);

        testRunNonIid.testCases.push_back(tc636);

        if (verbose > 0) printf("\nRunning Predictor Estimates...\n");

        // Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test

        NonIidTestCase tc637;
        tc637.testCaseNumber = "Multi Most Common in Window Test";
        tc637.data_word_size = data.word_size;

        ret_min_entropy = multi_mcw_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc637.h_r = ret_min_entropy;
            H_r = min(ret_min_entropy, H_r);
        }

        ret_min_entropy = multi_mcw_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc637.h_c = ret_min_entropy;
            H_c = min(ret_min_entropy, H_c);
        }

        testRunNonIid.testCases.push_back(tc637);
        // Section 6.3.8 - Estimate entropy with Lag Prediction Test

        NonIidTestCase tc638;
        tc638.testCaseNumber = "Lag Prediction Test";
        tc638.data_word_size = data.word_size;

        ret_min_entropy = lag_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tLag Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc638.h_r = ret_min_entropy;
            H_r = min(ret_min_entropy, H_r);
        }

        ret_min_entropy = lag_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tLag Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc638.h_c = ret_min_entropy;
            H_c = min(ret_min_entropy, H_c);
        }
        testRunNonIid.testCases.push_back(tc638);

        // Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)
        NonIidTestCase tc639;
        tc639.testCaseNumber = "Multi Markov Model with Counting Test (MultiMMC)";
        tc639.data_word_size = data.word_size;

        ret_min_entropy = multi_mmc_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc639.h_r = ret_min_entropy;
            H_r = min(ret_min_entropy, H_r);
        }

        ret_min_entropy = multi_mmc_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc639.h_c = ret_min_entropy;
            H_c = min(ret_min_entropy, H_c);
        }
        testRunNonIid.testCases.push_back(tc639);

        // Section 6.3.10 - Estimate entropy with LZ78Y Test

        NonIidTestCase tc6310;
        tc6310.testCaseNumber = "LZ78Y Test";
        tc6310.data_word_size = data.word_size;

        ret_min_entropy = LZ78Y_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tLZ78Y Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc6310.h_r = ret_min_entropy;
            H_r = min(ret_min_entropy, H_r);
        }

        ret_min_entropy = LZ78Y_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose > 1) printf("\tLZ78Y Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy, data.word_size);
            tc6310.h_c = ret_min_entropy;
            H_c = min(ret_min_entropy, H_c);
        }
        testRunNonIid.testCases.push_back(tc6310);

    } else { /* IID tests */

        // Compute chi square stats
        bool chi_square_test_pass_row = chi_square_tests(rdata, sample_size, alphabet_size, verbose);
        bool chi_square_test_pass_col = chi_square_tests(cdata, sample_size, alphabet_size, verbose);
        bool chi_square_test_pass = chi_square_test_pass_row && chi_square_test_pass_col;

        tcOverallIid.passed_chi_square_tests = chi_square_test_pass;

        if ((verbose == 1) || (verbose == 2)) {
            if (chi_square_test_pass) {
                printf("** Passed chi square tests\n\n");
            }
            else {
                printf("** Failed chi square tests\n\n");
            }
        }
        else if (verbose > 2) {
            if (chi_square_test_pass) {
                printf("Chi square tests: Passed\n");
            }
            else {
                printf("Chi square tests: Failed\n");
            }
        }

        // Compute length of the longest repeated substring stats
        bool len_LRS_test_pass_row = len_LRS_test(rdata, sample_size, alphabet_size, verbose, "Literal");
        bool len_LRS_test_pass_col = len_LRS_test(cdata, sample_size, alphabet_size, verbose, "Literal");
        bool len_LRS_test_pass = len_LRS_test_pass_row && len_LRS_test_pass_col;

        tcOverallIid.passed_longest_repeated_substring_test = len_LRS_test_pass;

        if ((verbose == 1) || (verbose == 2)) {
            if (len_LRS_test_pass) {
                printf("** Passed length of longest repeated substring test\n\n");
            }
            else {
                printf("** Failed length of longest repeated substring test\n\n");
            }
        }
        else if (verbose > 2) {
            if (len_LRS_test_pass) {
                printf("Length of longest repeated substring test: Passed\n");
            }
            else {
                printf("Length of longest repeated substring test: Failed\n");
            }
        }

        // Compute permutation stats
        bool perm_test_pass_row = permutation_tests(&data, rawmean, median, verbose, tcOverallIid);

        data_t data_col;
        memcpy(&data_col, &data, sizeof(data));
        data_col.symbols = rdata;

        bool perm_test_pass_col = permutation_tests(&data_col, rawmean, median, verbose, tcOverallIid);
        bool perm_test_pass = perm_test_pass_row && perm_test_pass_col;

        tcOverallIid.passed_iid_permutation_tests = perm_test_pass;

        if ((verbose == 1) || (verbose == 2)) {
            if (perm_test_pass) {
                printf("** Passed IID permutation tests\n\n");
            }
            else {
                printf("** Failed IID permutation tests\n\n");
            }
        }
        else if (verbose > 2) {
            if (perm_test_pass) {
                printf("IID permutation tests: Passed\n");
            }
            else {
                printf("IID permutation tests: Failed\n");
            }
        }


    }

    if (verbose > 0) {
        printf("\n");
        printf("H_r: %f\n", H_r);
        printf("H_c: %f\n", H_c);
        printf("H_I: %f\n", H_I);
        printf("\n");
    }

    if (min(H_r, H_c) < H_I / 2.0) {
        if (verbose > 0) printf("*** min(H_r, H_c) < H_I/2, Validation Testing Failed ***\n");
        if (jsonOutput) {
            if(iid) {
                testRunIid.errorLevel = -1;
                testRunIid.errorMsg = "min(H_r, H_c) < H_I/2, Validation Testing Failed.";
                ofstream output;
                output.open(outputfilename);
                output << testRunIid.GetAsJson();
                output.close();
            } else {
                testRunNonIid.errorLevel = -1;
                testRunNonIid.errorMsg = "min(H_r, H_c) < H_I/2, Validation Testing Failed.";
                ofstream output;
                output.open(outputfilename);
                output << testRunNonIid.GetAsJson();
                output.close();
            }
        }
        exit(-1);
    }


    testRunIid.testCases.push_back(tcOverallIid);
    testRunIid.errorLevel = 0;

    NonIidTestCase tcOverallNonIid;
    tcOverallNonIid.h_r = H_r;
    tcOverallNonIid.h_c = H_c;
    tcOverallNonIid.h_i = H_I;
    tcOverallNonIid.testCaseNumber = "Overall";
    testRunNonIid.testCases.push_back(tcOverallNonIid);
    testRunNonIid.errorLevel = 0;


    if (jsonOutput) {
        ofstream output;
        output.open(outputfilename);
        if (iid) {
            output << testRunIid.GetAsJson();
        } else {
            output << testRunNonIid.GetAsJson();
        }
        output.close();

    }
    if (verbose > 0) {
        printf("Validation Test Passed...\n\n");
        printf("min(H_r, H_c, H_I): %f\n\n", min(min(H_r, H_c), H_I));
    }
    free(cdata);
    free_data(&data);
    return 0;
}
