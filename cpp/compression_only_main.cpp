/* Compression Test Only - Fast entropy assessment using NIST SP800-90B compression test */

#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/TestRunUtils.h"
#include "non_iid/compression_test.h"

#include <getopt.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <openssl/sha.h>

[[ noreturn ]] void print_usage() {
    printf("Usage is: ea_compression_only [-i|-c] [-a|-t] [-v] [-q] [-l <index>,<samples> ] <file_name> [bits_per_symbol]\n\n");
    printf("\t <file_name>: Must be relative path to a binary file with at least 1 million entries (samples).\n");
    printf("\t [bits_per_symbol]: Must be between 1-8, inclusive. By default this value is inferred from the data.\n");
    printf("\t [-i|-c]: '-i' for initial entropy estimate, '-c' for conditioned sequential dataset entropy estimate. The initial entropy estimate is the default.\n");
    printf("\t [-a|-t]: '-a' produces the 'H_bitstring' assessment using all read bits, '-t' truncates the bitstring used to produce the `H_bitstring` assessment to %d bits. Test all data by default.\n", MIN_SIZE);
    printf("\t -v: Optional verbosity flag for more output. Can be used multiple times.\n");
    printf("\t -q: Quiet mode, less output to screen. This will override any verbose flags.\n");
    printf("\t -l <index>,<samples>\tRead the <index> substring of length <samples>.\n");
    printf("\n");
    printf("\t This tool runs ONLY the compression test from NIST SP800-90B for fast entropy assessment.\n");
    printf("\t For full validation, use the complete ea_non_iid tool.\n");
    printf("\n");
    printf("\t --version: Prints tool version information\n");
    printf("\n");
    exit(-1);
}

int main(int argc, char* argv[]) {

    bool initial_entropy, all_bits;
    int verbose = 1;
    bool quietMode = false;
    char *file_path;
    double H_original, H_bitstring, ret_min_entropy, h_assessed;
    data_t data;
    int opt;
    unsigned long subsetIndex = ULONG_MAX;
    unsigned long subsetSize = 0;
    unsigned long long inint;
    char *nextOption;

    data.word_size = 0;
    initial_entropy = true;
    all_bits = true;

    for (int i = 0; i < argc; i++) {
        std::string Str = std::string(argv[i]);
        if ("--version" == Str) {
            printVersion("compression_only");
            exit(0);
        }
    }

    while ((opt = getopt(argc, argv, "icatvql:")) != -1) {
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
            case 'q':
                quietMode = true;
                break;
            case 'l':
                inint = strtoull(optarg, &nextOption, 0);
                if ((inint > ULONG_MAX) || (errno == EINVAL) || (nextOption == NULL) || (*nextOption != ',')) {
                    printf("Error on index/samples.\n");
                    print_usage();
                }
                subsetIndex = inint;

                nextOption++;

                inint = strtoull(nextOption, NULL, 0);
                if ((inint > ULONG_MAX) || (errno == EINVAL)) {
                    printf("Error on index/samples.\n");
                    print_usage();
                }
                subsetSize = inint;
                break;
            default:
                print_usage();
        }
    }

    argc -= optind;
    argv += optind;

    // Parse args
    if ((argc != 1) && (argc != 2)) {
        printf("Incorrect usage.\n");
        print_usage();
    }

    // If quiet mode is enabled, force minimum verbose
    if (quietMode) {
        verbose = 0;
    }

    // get filename
    file_path = argv[0];

    char hash[2*SHA256_DIGEST_LENGTH+1];
    sha256_file(file_path, hash);

    if (argc == 2) {
        // get bits per word
        inint = atoi(argv[1]);
        if (inint < 1 || inint > 8) {
            printf("Invalid bits per symbol.\n");
            print_usage();
        } else {
            data.word_size = inint;
        }
    }

    if (verbose > 1) {
        if (subsetSize == 0) printf("Opening file: '%s' (SHA-256 hash %s)\n", file_path, hash);
        else printf("Opening file: '%s' (SHA-256 hash %s), reading block %ld of size %ld\n", file_path, hash, subsetIndex, subsetSize);
    }

    // Read the file (note: using NULL for testRun since we're not using JSON output)
    if (!read_file_subset(file_path, &data, subsetIndex, subsetSize, NULL)) {
        printf("Error reading file.\n");
        exit(-1);
    }

    if (verbose > 1) printf("Loaded %ld samples of %d distinct %d-bit-wide symbols\n", data.len, data.alph_size, data.word_size);

    if (data.alph_size <= 1) {
        printf("Symbol alphabet consists of 1 symbol. No entropy awarded...\n");
        free_data(&data);
        exit(-1);
    }

    if (!all_bits && (data.blen > MIN_SIZE)) data.blen = MIN_SIZE;

    if ((verbose > 1) && ((data.alph_size > 2) || !initial_entropy)) printf("Number of Binary Symbols: %ld\n", data.blen);
    if (data.len < MIN_SIZE) printf("\n*** Warning: data contains less than %d samples ***\n\n", MIN_SIZE);

    // The maximum min-entropy is -log2(1/2^word_size) = word_size
    // The maximum bit string min-entropy is 1.0
    H_original = data.word_size;
    H_bitstring = 1.0;

    if ((verbose == 1) || (verbose == 2)) {
        printf("\nRunning Compression Test Only...\n\n");
    }

    // Section 6.3.4 - Estimate entropy with Compression Test (for bit strings only)
    if (((data.alph_size > 2) || !initial_entropy)) {
        ret_min_entropy = compression_test(data.bsymbols, data.blen, verbose, "Bitstring");
        if (ret_min_entropy >= 0) {
            if (verbose >= 1) printf("Compression Test Estimate (bit string) = %f / 1 bit(s)\n", ret_min_entropy);
            H_bitstring = min(ret_min_entropy, H_bitstring);
        } else {
            printf("Compression test failed for bitstring data\n");
            free_data(&data);
            exit(-1);
        }
    }

    if (initial_entropy && (data.alph_size == 2)) {
        ret_min_entropy = compression_test(data.symbols, data.len, verbose, "Literal");
        if (ret_min_entropy >= 0) {
            if (verbose >= 1) printf("Compression Test Estimate = %f / 1 bit(s)\n", ret_min_entropy);
            H_original = min(ret_min_entropy, H_original);
        } else {
            printf("Compression test failed for literal data\n");
            free_data(&data);
            exit(-1);
        }
    }

    // Calculate assessed entropy
    h_assessed = data.word_size;

    if ((data.alph_size > 2) || !initial_entropy) {
        h_assessed = min(h_assessed, H_bitstring * data.word_size);
    }

    if (initial_entropy) {
        h_assessed = min(h_assessed, H_original);
    }

    // Output results
    if ((verbose == 1) || (verbose == 2)) {
        printf("\n=== COMPRESSION TEST RESULTS ===\n");
        if (initial_entropy) {
            printf("H_original: %f\n", H_original);
            if (data.alph_size > 2) {
                printf("H_bitstring: %f\n", H_bitstring);
                printf("min(H_original, %d X H_bitstring): %f\n", data.word_size, min(H_original, data.word_size * H_bitstring));
            }
        } else {
            printf("h': %f\n", H_bitstring);
        }
        printf("\nAssessed Min-Entropy (Compression Test Only): %f\n", h_assessed);
        printf("=================================\n");
    } else if (verbose > 2) {
        if ((data.alph_size > 2) || !initial_entropy) {
            printf("H_bitstring = %.17g\n", H_bitstring);
            printf("H_bitstring Per Symbol = %.17g\n", H_bitstring * data.word_size);
        }

        if (initial_entropy) {
            printf("H_original = %.17g\n", H_original);
        }

        printf("Assessed min entropy (compression only): %.17g\n", h_assessed);
    } else if (verbose == 0) {
        // Quiet mode - just output the final result
        printf("%.6f\n", h_assessed);
    }

    free_data(&data);
    return 0;
}