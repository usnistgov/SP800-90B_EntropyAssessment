/* VERSION information is kept in utils.h. Please update when a new version is released */


#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "iid/iid_test_run.h"
#include "shared/TestRunUtils.h"
#include "iid/permutation_tests.h"
#include "iid/chi_square_tests.h"
#include <openssl/sha.h>
#include <omp.h>
#include <getopt.h>
#include <limits.h>

#include <iostream>
#include <fstream>


[[ noreturn ]] void print_usage() {
    printf("Usage is: ea_iid [-i|-c] [-a|-t] [-v] [-q] [-l <index>,<samples> ] <file_name> [bits_per_symbol]\n\n");
    printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (samples).\n");
    printf("\t [bits_per_symbol]: Must be between 1-8, inclusive. By default this value is inferred from the data.\n");
    printf("\t [-i|-c]: '-i' for initial entropy estimate, '-c' for conditioned sequential dataset entropy estimate. The initial entropy estimate is the default.\n");
    printf("\t [-a|-t]: '-a' produces the 'H_bitstring' assessment using all read bits, '-t' truncates the bitstring used to produce the `H_bitstring` assessment to %d bits. Test all data by default.\n", MIN_SIZE);
    printf("\t Note: When testing binary data, no `H_bitstring` assessment is produced, so the `-a` and `-t` options produce the same results for the initial assessment of binary data.\n");
    printf("\t -v: Optional verbosity flag for more output. Can be used multiple times.\n");
    printf("\t -q: Quiet mode, less output to screen. This will override any verbose flags.\n");
    printf("\t -l <index>,<samples>\tRead the <index> substring of length <samples>.\n");
    printf("\n");
    printf("\t Samples are assumed to be packed into 8-bit values, where the least significant 'bits_per_symbol'\n");
    printf("\t bits constitute the symbol.\n");
    printf("\n");
    printf("\t -i: Initial Entropy Estimate (Section 3.1.3)\n");
    printf("\n");
    printf("\t\t Computes the initial entropy estimate H_I as described in Section 3.1.3\n");
    printf("\t\t (not accounting for H_submitter) using the entropy estimators specified in\n");
    printf("\t\t Section 6.3.  If 'bits_per_symbol' is greater than 1, the samples are also\n");
    printf("\t\t converted to bitstrings and assessed to create H_bitstring; for multi-bit symbols,\n");
    printf("\t\t two entropy estimates are computed: H_original and H_bitstring.\n");
    printf("\t\t Returns min(H_original, bits_per_symbol X H_bitstring). The initial entropy\n");
    printf("\t\t estimate H_I = min(H_submitter, H_original, bits_per_symbol X H_bitstring).\n");
    printf("\n");
    printf("\t -c: Conditioned Sequential Dataset Entropy Estimate (Section 3.1.5.2)\n");
    printf("\n");
    printf("\t\t Computes the entropy estimate per bit h' for the conditioned sequential dataset if the\n");
    printf("\t\t conditioning function is non-vetted. The samples are converted to a bitstring.\n");
    printf("\t\t Returns h' = min(H_bitstring).\n");
    printf("\n");
    printf("\t -o: Set Output Type to JSON\n");
    printf("\n");
    printf("\t\t Changes the output format to JSON and sets the file location for the output file.\n");
    printf("\n");
    printf("\t --version: Prints tool version information");
    printf("\n");
    exit(-1);
}

int main(int argc, char* argv[]) {

    bool initial_entropy, all_bits;
    int verbose = 1; //verbose 0 is for JSON output, 1 is the normal mode, 2 is the NIST tool verbose mode, and 3 is for extra verbose output
    double rawmean, median;
    char* file_path;
    data_t data;
    int opt;
    unsigned long subsetIndex = ULONG_MAX;
    unsigned long subsetSize = 0;
    unsigned long long inint;
    char *nextOption;

    data.word_size = 0;
    initial_entropy = true;
    all_bits = true;

    bool quietMode = false;
    bool jsonOutput = false;
    string timestamp = getCurrentTimestamp();
    string outputfilename;
    string commandline = recreateCommandLine(argc, argv);

    IidTestRun testRun;
    testRun.timestamp = timestamp;
    testRun.commandline = commandline;
    
    for (int i = 0; i < argc; i++) {
        std::string Str = std::string(argv[i]);
        if ("--version" == Str) {
            printVersion("iid");
            exit(0);
        }
    }

    while ((opt = getopt(argc, argv, "icatvl:qo:")) != -1) {
        switch (opt) {
            case 'i':
                initial_entropy = true;
                break;
            case 'c':
                initial_entropy = false;
                break;
            case 'a':
                all_bits = true;
                break;
            case 't':
                all_bits = false;
                break;
            case 'v':
                verbose++;
                break;
            case 'l':
                inint = strtoull(optarg, &nextOption, 0);
                if ((inint > ULONG_MAX) || (errno == EINVAL) || (nextOption == NULL) || (*nextOption != ',')) {

                    testRun.errorLevel = -1;
                    testRun.errorMsg = "Error on index/samples.";

                    if (jsonOutput) {
                        ofstream output;
                        output.open(outputfilename);
                        output << testRun.GetAsJson();
                        output.close();
                    }
                    print_usage();
                }
                subsetIndex = inint;

                nextOption++;

                inint = strtoull(nextOption, NULL, 0);
                if ((inint > ULONG_MAX) || (errno == EINVAL)) {
                    testRun.errorLevel = -1;
                    testRun.errorMsg = "Error on index/samples.";

                    if (jsonOutput) {
                        ofstream output;
                        output.open(outputfilename);
                        output << testRun.GetAsJson();
                        output.close();
                    }
                    print_usage();
                }

                subsetSize = inint;
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
    if ((argc != 2) && (argc != 1)) {
        printf("Incorrect usage.\n");
        print_usage();
    }

    // If quiet mode is enabled, force minimum verbose
    if (quietMode) {
        verbose = 0;
    }

    // get filename
    file_path = argv[0];

    testRun.filename = file_path;

    if (argc == 2) {
        // get bits per word
        data.word_size = atoi(argv[1]);
        if (data.word_size < 1 || data.word_size > 8) {

            testRun.errorLevel = -1;
            testRun.errorMsg = "Invalid bits per symbol: " + std::to_string(data.word_size) + ".";

            if (jsonOutput) {
                ofstream output;
                output.open(outputfilename);
                output << testRun.GetAsJson();
                output.close();
            }

            printf("Invalid bits per symbol: %d.\n", data.word_size);
            print_usage();
        }
    }

    // Record hash of input file
    char hash[2*SHA256_DIGEST_LENGTH+1];
    sha256_file(file_path, hash);
    testRun.sha256 = hash;

    if (verbose > 1) {
        if (subsetSize == 0) {
            printf("Opening file: '%s' (SHA-256 hash %s)\n", file_path, hash);
        } else {
            printf("Opening file: '%s' (SHA-256 hash %s), reading block %ld of size %ld\n", file_path, hash, subsetIndex, subsetSize);
        }
    }
    if (!read_file_subset(file_path, &data, subsetIndex, subsetSize, &testRun)) {
        if (jsonOutput) {
            ofstream output;
            output.open(outputfilename);
            output << testRun.GetAsJson();
            output.close();
        }

        printf("Error reading file.\n");
        print_usage();
    }

    if (verbose > 1) printf("Loaded %ld samples of %d distinct %d-bit-wide symbols\n", data.len, data.alph_size, data.word_size);

    if (data.alph_size <= 1) {

        testRun.errorLevel = -1;
        testRun.errorMsg = "Symbol alphabet consists of 1 symbol. No entropy awarded...";

        if (jsonOutput) {
            ofstream output;
            output.open(outputfilename);
            output << testRun.GetAsJson();
            output.close();
        }

        printf("Symbol alphabet consists of 1 symbol. No entropy awarded...\n");
        free_data(&data);
        exit(-1);
    }

    if (!all_bits && (data.blen > MIN_SIZE)) data.blen = MIN_SIZE;

    if ((verbose > 1) && ((data.alph_size > 2) || !initial_entropy)) printf("Number of Binary samples: %ld\n", data.blen);
    if (data.len < MIN_SIZE) printf("\n*** Warning: data contains less than %d samples ***\n\n", MIN_SIZE);
    if (verbose > 1) {
        if (data.alph_size < (1 << data.word_size)) printf("\nSamples have been translated\n");
    }

    // Calculate baseline statistics
    int alphabet_size = data.alph_size;
    int sample_size = data.len;

    if ((verbose == 1) || (verbose == 2))
        printf("Calculating baseline statistics...\n");

    calc_stats(&data, rawmean, median);

    if (verbose == 2) {
        printf("\tRaw Mean: %f\n", rawmean);
        printf("\tMedian: %f\n", median);
        printf("\tBinary: %s\n\n", (alphabet_size == 2 ? "true" : "false"));
    } else if (verbose > 2) {
        printf("Raw Mean = %.17g\n", rawmean);
        printf("Median = %.17g\n", median);
        printf("Binary = %s\n", (alphabet_size == 2 ? "true" : "false"));
    }

    IidTestCase tc;
    tc.mean = rawmean;
    tc.median = median;
    tc.binary = (alphabet_size == 2);

    double H_original = data.word_size;
    double H_bitstring = 1.0;

    // Compute the min-entropy of the dataset
    if (initial_entropy) {
        H_original = most_common(data.symbols, sample_size, alphabet_size, verbose, "Literal");
    }
    tc.h_original = H_original;

    if (((data.alph_size > 2) || !initial_entropy)) {
        H_bitstring = most_common(data.bsymbols, data.blen, 2, verbose, "Bitstring");
    }
    tc.h_bitstring = H_bitstring;

    double h_assessed = data.word_size;
    if ((verbose == 1) || (verbose == 2)) {
        if (initial_entropy) {
            printf("H_original: %f\n", H_original);
            if (data.alph_size > 2) {
                printf("H_bitstring: %f\n", H_bitstring);
                printf("min(H_original, %d X H_bitstring): %f\n", data.word_size, min(H_original, data.word_size * H_bitstring));
            }
        } else {
            printf("h': %f\n", H_bitstring);
        }
    } else if (verbose > 2) {
        h_assessed = data.word_size;

        if ((data.alph_size > 2) || !initial_entropy) {
            h_assessed = min(h_assessed, H_bitstring * data.word_size);
            printf("H_bitstring = %.17g\n", H_bitstring);
            printf("H_bitstring Per Symbol = %.17g\n", H_bitstring * data.word_size);
        }

        if (initial_entropy) {
            h_assessed = min(h_assessed, H_original);
            printf("H_original = %.17g\n", H_original);
        }

        printf("Assessed min entropy: %.17g\n", h_assessed);
    }
    tc.h_assessed = h_assessed;

    // Compute chi square stats
    bool chi_square_test_pass = chi_square_tests(data.symbols, sample_size, alphabet_size, verbose);
    tc.passed_chi_square_tests = chi_square_test_pass;

    if ((verbose == 1) || (verbose == 2)) {
        if (chi_square_test_pass) {
            printf("** Passed chi square tests\n\n");
        } else {
            printf("** Failed chi square tests\n\n");
        }
    } else if (verbose > 2) {
        if (chi_square_test_pass) {
            printf("Chi square tests: Passed\n");
        } else {
            printf("Chi square tests: Failed\n");
        }
    }

    // Compute length of the longest repeated substring stats
    bool len_LRS_test_pass = len_LRS_test(data.symbols, sample_size, alphabet_size, verbose, "Literal");
    tc.passed_longest_repeated_substring_test = len_LRS_test_pass;

    if ((verbose == 1) || (verbose == 2)) {
        if (len_LRS_test_pass) {
            printf("** Passed length of longest repeated substring test\n\n");
        } else {
            printf("** Failed length of longest repeated substring test\n\n");
        }
    } else if (verbose > 2) {
        if (len_LRS_test_pass) {
            printf("Length of longest repeated substring test: Passed\n");
        } else {
            printf("Length of longest repeated substring test: Failed\n");
        }
    }

    // Compute permutation stats
    bool perm_test_pass = permutation_tests(&data, rawmean, median, verbose, tc);
    tc.passed_iid_permutation_tests = perm_test_pass;

    if ((verbose == 1) || (verbose == 2)) {
        if (perm_test_pass) {
            printf("** Passed IID permutation tests\n\n");
        } else {
            printf("** Failed IID permutation tests\n\n");
        }
    } else if (verbose > 2) {
        if (perm_test_pass) {
            printf("IID permutation tests: Passed\n");
        } else {
            printf("IID permutation tests: Failed\n");
        }
    }

    testRun.testCases.push_back(tc);
    testRun.errorLevel = 0;

    if (jsonOutput) {
        ofstream output;
        output.open(outputfilename);
        output << testRun.GetAsJson();
        output.close();
    }

    free_data(&data);
    return 0;
}
