#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
//#include "non_iid/non_iid_test_run.h"
//#include "iid/iid_test_run.h"
#include "shared/TestRun.h"
#include "shared/TestRunUtils.h"
#include "non_iid/collision_test.h"
#include "non_iid/lz78y_test.h"
#include "non_iid/multi_mmc_test.h"
#include "non_iid/lag_test.h"
#include "non_iid/multi_mcw_test.h"
#include "non_iid/compression_test.h"
#include "non_iid/markov_test.h"
#include "other/restart_test_run.h"

#include <getopt.h>
#include <limits.h>
#include <iostream>
#include <fstream>

//Each test has a targeted chance of roughly 0.000005, and we need to witness at least 5 failures, so this should be no less than 1000000
#define SIMULATION_ROUNDS 5000000

[[ noreturn ]] void print_usage() {
    printf("Usage is: ea_restart [-i|-n] [-v] [-q] <file_name> [bits_per_symbol] <H_I> [-o output_file.json]\n\n");
    printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (samples),\n");
    printf("\t and in the \"row dataset\" format described in SP800-90B Section 3.1.4.1.\n");
    printf("\t [bits_per_symbol]: Must be between 1-8, inclusive.\n");
    printf("\t <H_I>: Initial entropy estimate.\n");
    printf("\t [-i|-n]: '-i' for IID data, '-n' for non-IID data. Non-IID is the default.\n");
    printf("\t -v: Optional verbosity flag for more output.\n");
    printf("\t -q: Quiet mode, less output to screen.\n");
    printf("\n");
    printf("\t -o: Set Output Type to JSON\n");
    printf("\n");
    printf("\t\t Changes the output format to JSON and sets the file location for the output file.\n");
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
    exit(-1);
}

//Here, we simulate a sort of "worst case" for this test, where there are a maximal number of symbols with maximal probability,
//and the rest is distributed to the other symbols

long int simulateCount(int k, double H_I, uint64_t *xoshiro256starstarState) {
    long int counts[256];
    int current_symbol;
    int k_max;
    long int max_count = 0;
    double p, max_cutoff, p_min, cur_rand;

    assert(k <= 256);

    p = pow(2.0, -H_I);

    k_max = floor(1.0 / p);

    assert(k_max <= k);


    if (k > k_max) {
        max_cutoff = p * k_max;
        p_min = (1.0 - max_cutoff) / (k - k_max);
    } else {
        max_cutoff = 1.0;
        p_min = 0.0;
    }

    for (int j = 0; j < k; j++) counts[j] = 0;

    for (int j = 0; j < 1000; j++) {
        cur_rand = randomUnit(xoshiro256starstarState);
        if (cur_rand < max_cutoff) {
            current_symbol = floor(cur_rand / p);
            assert((current_symbol >= 0) && (current_symbol < k_max));
        } else {
            current_symbol = floor((cur_rand - max_cutoff) / p_min) + k_max;
            assert((current_symbol >= k_max) && (current_symbol < k));
        }
        counts[current_symbol]++;
    }

    for (int j = 0; j < k; j++) {
        if (max_count < counts[j]) max_count = counts[j];
    }

    return max_count;
}

//This returns the bound (cutoff) for the test. Counts equal to this value should pass.
//Larger values should fail.

long int simulateBound(double alpha, int k, double H_I) {
    uint64_t xoshiro256starstarMainSeed[4];
    vector<long int> results(SIMULATION_ROUNDS, -1);
    long int returnIndex;

    assert((k > 1) && (k <= 256));

    seed(xoshiro256starstarMainSeed);

#pragma omp parallel
    {
        uint64_t xoshiro256starstarSeed[4];

        memcpy(xoshiro256starstarSeed, xoshiro256starstarMainSeed, sizeof (xoshiro256starstarMainSeed));
        //Cause the RNG to jump omp_get_thread_num() * 2^128 calls
        xoshiro_jump(omp_get_thread_num(), xoshiro256starstarSeed);

#pragma omp for
        for (int i = 0; i < SIMULATION_ROUNDS; i++) {
            results[i] = simulateCount(k, H_I, xoshiro256starstarSeed);
        }
    }

    sort(results.begin(), results.end());
    assert((results[0] >= (1000 / k)) && (results[0] <= 1000));
    assert((results[SIMULATION_ROUNDS - 1] >= (1000 / k)) && (results[SIMULATION_ROUNDS - 1] <= 1000));
    assert(results[0] <= results[SIMULATION_ROUNDS - 1]);

    returnIndex = ((size_t) floor((1.0 - alpha) * ((double) SIMULATION_ROUNDS))) - 1;
    assert((returnIndex >= 0) && (returnIndex < SIMULATION_ROUNDS));

    return (results[returnIndex]);
}

int main(int argc, char* argv[]) {
    bool iid;
    int verbose = 0;
    bool quietMode = false;
    char *file_path;
    int r = 1000, c = 1000;
    int counts[256];
    long int X_cutoff;
    long i, j, X_i, X_r, X_c, X_max;
    double H_I, H_r, H_c, alpha, ret_min_entropy_row, ret_min_entropy_col;
    byte *rdata, *cdata;
    data_t data;
    int opt;

    iid = false;
    data.word_size = 0;

    bool jsonOutput = false;
    string timestamp = getCurrentTimestamp();
    string outputfilename = timestamp + ".json";

    while ((opt = getopt(argc, argv, "invqo:")) != -1) {
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

    char hash[65];
    sha256_file(file_path, hash);

    RestartTestRun testRunRestart;
    testRunRestart.type = "Restart";
    testRunRestart.timestamp = timestamp;
    testRunRestart.filename = file_path;
    testRunRestart.sha256 = hash;
    testRunRestart.IID = iid;
    
    //NonIidTestRun testRunRestart;
    //testRunRestart.type = "Restart";
    //testRunRestart.timestamp = timestamp;
    //testRunRestart.sha256 = hash;
    //testRunRestart.filename = file_path;

    if (argc == 2) {
        // get bits per word
        data.word_size = atoi(argv[0]);
        if (data.word_size < 1 || data.word_size > 8) {
            printf("Invalid bits per symbol.\n");
            if (jsonOutput) {
                testRunRestart.errorLevel = -1;
                testRunRestart.errorMsg = "Invalid bits per symbol.";
                ofstream output;
                output.open(outputfilename);
                output << testRunRestart.GetAsJson();
                output.close();
            }
            print_usage();
        }
        argv++;
        argc--;
    }

    // get H_I	
    H_I = atof(argv[0]);
    if (H_I < 0) {
        printf("H_I must be nonnegative.\n");
        if (jsonOutput) {
            testRunRestart.errorLevel = -1;
            testRunRestart.errorMsg = "H_I must be nonnegative.";
            ofstream output;
            output.open(outputfilename);
            output << testRunRestart.GetAsJson();
            output.close();
        }

        print_usage();
    }

    if (verbose > 0) 
    	printf("Opening file: '%s'\n", file_path);

    if (!read_file(file_path, &data)) {
        printf("Error reading file.\n");

        if (jsonOutput) {
            testRunRestart.errorLevel = -1;
            testRunRestart.errorMsg = "Error reading file.";
            ofstream output;
            output.open(outputfilename);
            output << testRunRestart.GetAsJson();
            output.close();
        }

        print_usage();
    }

    if (verbose > 0) 
    	printf("Loaded %ld samples made up of %d distinct %d-bit-wide symbols.\n", data.len, data.alph_size, data.word_size);

    if (H_I > data.word_size) {
        printf("H_I must be at most 'bits_per_symbol'.\n");
        if (jsonOutput) {
            testRunRestart.errorLevel = -1;
            testRunRestart.errorMsg = "H_I must be at most 'bits_per_symbol'.";
            ofstream output;
            output.open(outputfilename);
            output << testRunRestart.GetAsJson();
            output.close();
        }

        free_data(&data);
        exit(-1);
    }

    if (data.alph_size <= 1) {
        printf("Symbol alphabet consists of 1 symbol. No entropy awarded...\n");
        if (jsonOutput) {
            testRunRestart.errorLevel = -1;
            testRunRestart.errorMsg = "Symbol alphabet consists of 1 symbol. No entropy awarded...";
            ofstream output;
            output.open(outputfilename);
            output << testRunRestart.GetAsJson();
            output.close();
        }

        free_data(&data);
        exit(-1);
    }

    if (data.len != MIN_SIZE) {
        printf("\n*** Error: data does not contain %d samples ***\n\n", MIN_SIZE);
        if (jsonOutput) {
            testRunRestart.errorLevel = -1;
            testRunRestart.errorMsg = "*** Error: data does not contain " + std::to_string(MIN_SIZE) + " samples ***";
            ofstream output;
            output.open(outputfilename);
            output << testRunRestart.GetAsJson();
            output.close();
        }

        exit(-1);
    }
    if (verbose > 0) {
        if (data.alph_size < (1 << data.word_size)) 
        	printf("\nSymbols have been translated.\n\n");
    }

    rdata = data.symbols;
    cdata = (byte*) malloc(data.len);
    if (cdata == NULL) {
        printf("Error: failure to initialize memory for columns\n");
        if (jsonOutput) {
            testRunRestart.errorLevel = -1;
            testRunRestart.errorMsg = "Error: failure to initialize memory for columns";
            ofstream output;
            output.open(outputfilename);
            output << testRunRestart.GetAsJson();
            output.close();
        }

        exit(-1);
    }

    if (!quietMode)
    	printf("H_I: %f\n", H_I);

    alpha = 1 - exp(log(0.99) / (r + c));
    X_cutoff = simulateBound(alpha, data.alph_size, H_I);
    if (!quietMode) 
    	printf("ALPHA: %.17g, X_cutoff: %ld\n", alpha, X_cutoff);

    // get maximum row count
    X_r = 0;
    for (i = 0; i < r; i++) { //row
        memset(counts, 0, 256 * sizeof (int));
        X_i = 0;
        for (j = 0; j < c; j++) {//column
            
            //[i*r+j] is row i, column j
            //So, we're fixing a row, and then iterate through various columns
            if (++counts[rdata[i * r + j]] > X_i) 
            	X_i = counts[rdata[i * r + j]];
        }

        if (X_i > X_r) 
        	X_r = X_i;
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
            if (++counts[cdata[j * c + i]] > X_i) 
            	X_i = counts[cdata[j * c + i]];
        }

        if (X_i > X_c) 
        	X_c = X_i;
    }

    // perform sanity check on rows and columns of restart data (Section 3.1.4.3)
    X_max = max(X_r, X_c);
    if (!quietMode) 
    	printf("X_max: %ld\n", X_max);

    if (X_max > X_cutoff) {
        printf("\n*** Restart Sanity Check Failed ***\n");
        exit(-1);
    } 
    else if (verbose > 0) 
    	printf("\nRestart Sanity Check Passed...\n");

    // The maximum min-entropy is -log2(1/2^word_size) = word_size
    H_c = data.word_size;
    H_r = data.word_size;
    if (!quietMode) {
        if (iid) 
        	printf("\nRunning IID tests...\n\n");
        else 
        	printf("\nRunning non-IID tests...\n\n");

        printf("Running Most Common Value Estimate...\n");
    }

    // Section 6.3.1 - Estimate entropy with Most Common Value
    ret_min_entropy_row = most_common(rdata, data.len, data.alph_size, verbose, "Literal");
    if (verbose > 0) printf("\tMost Common Value Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy_row, data.word_size);
    H_r = min(ret_min_entropy_row, H_r);

    ret_min_entropy_col = most_common(cdata, data.len, data.alph_size, verbose, "Literal");
    if (verbose > 0) printf("\tMost Common Value Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy_col, data.word_size);
    H_c = min(ret_min_entropy_col, H_c);
/*
    NonIidTestCase tc631;
    tc631.ret_min_entropy_row = ret_min_entropy_row;
    tc631.ret_min_entropy_col = ret_min_entropy_col;
    tc631.data_word_size = data.word_size;
    tc631.h_r = H_r;
    tc631.h_c = H_c;
    tc631.testCaseNumber = "Estimate entropy with Most Common Value";
    testRunRestart.testCases.push_back(tc631);
*/
    if (!iid) {

        if (data.alph_size == 2) {
            if (!quietMode) printf("\nRunning Entropic Statistic Estimates (bit strings only)...\n");

            // Section 6.3.2 - Estimate entropy with Collision Test (for bit strings only)
            ret_min_entropy_row = collision_test(rdata, data.len, verbose, "Literal");
            if (verbose > 0) printf("\tCollision Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy_row);
            H_r = min(ret_min_entropy_row, H_r);

            ret_min_entropy_col = collision_test(cdata, data.len, verbose, "Literal");
            if (verbose > 0) printf("\tCollision Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy_col);
            H_c = min(ret_min_entropy_col, H_c);
/*
            NonIidTestCase tc632;
            tc632.ret_min_entropy_row = ret_min_entropy_row;
            tc632.ret_min_entropy_col = ret_min_entropy_col;
            tc632.data_word_size = data.word_size;
            tc632.h_r = H_r;
            tc632.h_c = H_c;
            tc632.testCaseNumber = "Estimate entropy with Collision Test (for bit strings only)";
            testRunRestart.testCases.push_back(tc632);
*/
            // Section 6.3.3 - Estimate entropy with Markov Test (for bit strings only)
            ret_min_entropy_row = markov_test(rdata, data.len, verbose, "Literal");
            if (verbose > 0) printf("\tMarkov Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy_row);
            H_r = min(ret_min_entropy_row, H_r);

            ret_min_entropy_col = markov_test(cdata, data.len, verbose, "Literal");
            if (verbose > 0) printf("\tMarkov Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy_col);
            H_c = min(ret_min_entropy_col, H_c);
/*
            NonIidTestCase tc633;
            tc633.ret_min_entropy_row = ret_min_entropy_row;
            tc633.ret_min_entropy_col = ret_min_entropy_col;
            tc633.data_word_size = data.word_size;
            tc633.h_r = H_r;
            tc633.h_c = H_c;
            tc633.testCaseNumber = "Estimate entropy with Markov Test (for bit strings only)";
            testRunRestart.testCases.push_back(tc633);
*/
            // Section 6.3.4 - Estimate entropy with Compression Test (for bit strings only)
            ret_min_entropy_row = compression_test(rdata, data.len, verbose, "Literal");
            if (ret_min_entropy_row >= 0) {
                if (verbose > 0) printf("\tCompression Test Estimate (Rows) = %f / 1 bit(s)\n", ret_min_entropy_row);
                H_r = min(ret_min_entropy_row, H_r);
            }

            ret_min_entropy_col = compression_test(cdata, data.len, verbose, "Literal");
            if (ret_min_entropy_col >= 0) {
                if (verbose > 0) printf("\tCompression Test Estimate (Cols) = %f / 1 bit(s)\n", ret_min_entropy_col);
                H_c = min(ret_min_entropy_col, H_c);
            }
/*
            NonIidTestCase tc634;
            tc634.ret_min_entropy_row = ret_min_entropy_row;
            tc634.ret_min_entropy_col = ret_min_entropy_col;
            tc634.data_word_size = data.word_size;
            tc634.h_r = H_r;
            tc634.h_c = H_c;
            tc634.testCaseNumber = "Estimate entropy with Compression Test (for bit strings only)";
            testRunRestart.testCases.push_back(tc634);
*/
        }

        // Section 6.3.5 - Estimate entropy with t-Tuple Test
        if (!quietMode) printf("\nRunning Tuple Estimates...\n");

        double row_t_tuple_res, row_lrs_res;
        double col_t_tuple_res, col_lrs_res;
        SAalgs(rdata, data.len, data.alph_size, row_t_tuple_res, row_lrs_res, verbose, "Literal");
        SAalgs(cdata, data.len, data.alph_size, col_t_tuple_res, col_lrs_res, verbose, "Literal");

        if (verbose > 0) printf("\tT-Tuple Test Estimate (Rows) = %f / %d bit(s)\n", row_t_tuple_res, data.word_size);
        H_r = min(row_t_tuple_res, H_r);

        if (verbose > 0) printf("\tT-Tuple Test Estimate (Cols) = %f / %d bit(s)\n", col_t_tuple_res, data.word_size);
        H_c = min(col_t_tuple_res, H_c);
/*
        NonIidTestCase tc635;
        tc635.ret_min_entropy_row = row_t_tuple_res;
        tc635.ret_min_entropy_col = col_t_tuple_res;
        tc635.data_word_size = data.word_size;
        tc635.h_r = H_r;
        tc635.h_c = H_c;
        tc635.testCaseNumber = "Estimate entropy with t-Tuple Test";
        testRunRestart.testCases.push_back(tc635);
*/
        // Section 6.3.6 - Estimate entropy with LRS Test
        if (verbose > 0) printf("\tLRS Test Estimate (Rows) = %f / %d bit(s)\n", row_lrs_res, data.word_size);
        H_r = min(row_lrs_res, H_r);

        if (verbose > 0) printf("\tLRS Test Estimate (Cols) = %f / %d bit(s)\n", col_lrs_res, data.word_size);
        H_c = min(col_lrs_res, H_c);
/*
        NonIidTestCase tc636;
        tc636.ret_min_entropy_row = row_lrs_res;
        tc636.ret_min_entropy_col = col_lrs_res;
        tc636.data_word_size = data.word_size;
        tc636.h_r = H_r;
        tc636.h_c = H_c;
        tc636.testCaseNumber = "Estimate entropy with LRS Test";
        testRunRestart.testCases.push_back(tc636);
*/

        // Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test
        if (!quietMode) printf("\nRunning Predictor Estimates...\n");

        ret_min_entropy_row = multi_mcw_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_row >= 0) {
            if (verbose > 0) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy_row, data.word_size);
            H_r = min(ret_min_entropy_row, H_r);
        }

        ret_min_entropy_col = multi_mcw_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_col >= 0) {
            if (verbose > 0) printf("\tMulti Most Common in Window (MultiMCW) Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy_col, data.word_size);
            H_c = min(ret_min_entropy_col, H_c);
        }
/*
        NonIidTestCase tc637;
        tc637.ret_min_entropy_row = ret_min_entropy_row;
        tc637.ret_min_entropy_col = ret_min_entropy_col;
        tc637.data_word_size = data.word_size;
        tc637.h_r = H_r;
        tc637.h_c = H_c;
        tc637.testCaseNumber = "Estimate entropy with Multi Most Common in Window Test";
        testRunRestart.testCases.push_back(tc637);
*/
        // Section 6.3.8 - Estimate entropy with Lag Prediction Test
        ret_min_entropy_row = lag_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_row >= 0) {
            if (verbose > 0) printf("\tLag Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy_row, data.word_size);
            H_r = min(ret_min_entropy_row, H_r);
        }

        ret_min_entropy_col = lag_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_col >= 0) {
            if (verbose > 0) printf("\tLag Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy_col, data.word_size);
            H_c = min(ret_min_entropy_col, H_c);
        }
/*
        NonIidTestCase tc638;
        tc638.ret_min_entropy_row = ret_min_entropy_row;
        tc638.ret_min_entropy_col = ret_min_entropy_col;
        tc638.data_word_size = data.word_size;
        tc638.h_r = H_r;
        tc638.h_c = H_c;
        tc638.testCaseNumber = "Estimate entropy with Lag Prediction Test";
        testRunRestart.testCases.push_back(tc638);
*/
        // Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)
        ret_min_entropy_row = multi_mmc_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_row >= 0) {
            if (verbose > 0) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy_row, data.word_size);
            H_r = min(ret_min_entropy_row, H_r);
        }

        ret_min_entropy_col = multi_mmc_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_col >= 0) {
            if (verbose > 0) printf("\tMulti Markov Model with Counting (MultiMMC) Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy_col, data.word_size);
            H_c = min(ret_min_entropy_col, H_c);
        }
/*
        NonIidTestCase tc639;
        tc639.ret_min_entropy_row = ret_min_entropy_row;
        tc639.ret_min_entropy_col = ret_min_entropy_col;
        tc639.data_word_size = data.word_size;
        tc639.h_r = H_r;
        tc639.h_c = H_c;
        tc639.testCaseNumber = "Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)";
        testRunRestart.testCases.push_back(tc639);
*/
        // Section 6.3.10 - Estimate entropy with LZ78Y Test
        ret_min_entropy_row = LZ78Y_test(rdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_row >= 0) {
            if (verbose > 0) printf("\tLZ78Y Prediction Test Estimate (Rows) = %f / %d bit(s)\n", ret_min_entropy_row, data.word_size);
            H_r = min(ret_min_entropy_row, H_r);
        }

        ret_min_entropy_col = LZ78Y_test(cdata, data.len, data.alph_size, verbose, "Literal");
        if (ret_min_entropy_col >= 0) {
            if (verbose > 0) printf("\tLZ78Y Prediction Test Estimate (Cols) = %f / %d bit(s)\n", ret_min_entropy_col, data.word_size);
            H_c = min(ret_min_entropy_col, H_c);
        }
/*
        NonIidTestCase tc6310;
        tc6310.ret_min_entropy_row = ret_min_entropy_row;
        tc6310.ret_min_entropy_col = ret_min_entropy_col;
        tc6310.data_word_size = data.word_size;
        tc6310.h_r = H_r;
        tc6310.h_c = H_c;
        tc6310.testCaseNumber = "Estimate entropy with LZ78Y Test";
        testRunRestart.testCases.push_back(tc6310);
*/
    }
    if (!quietMode) {
        printf("\n");
        printf("H_r: %f\n", H_r);
        printf("H_c: %f\n", H_c);
        printf("H_I: %f\n", H_I);
        printf("\n");
    }

    RestartTestCase tcOverall;
    tcOverall.h_r = H_r;
    tcOverall.h_c = H_c;
    //tcOverall.h_i = H_I;
    tcOverall.testCaseNumber = "Overall";
    testRunRestart.testCases.push_back(tcOverall);
    testRunRestart.errorLevel = 0;

    if (jsonOutput) {
        ofstream output;
        output.open(outputfilename);
        if (iid) {
            output << testRunRestart.GetAsJson();
        } else {
            output << testRunRestart.GetAsJson();
        }
        output.close();

    }
    if (!quietMode) {
        if (min(H_r, H_c) < H_I / 2.0) printf("*** min(H_r, H_c) < H_I/2, Validation Testing Failed ***\n");
        else {
            printf("Validation Test Passed...\n\n");
            printf("min(H_r, H_c, H_I): %f\n\n", min(min(H_r, H_c), H_I));
        }
    }
    
    free(cdata);
    free_data(&data);
    return 0;
}
