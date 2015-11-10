# DRAFT NIST SP 800-90B (August 2012) Section 9.3 tests
#
# Estimating the Min-Entropy of non-IID Sources
#
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness.
#
# T. A. Hall
# tim.hall@nist.gov
#
# 8 September 2014

import sys
import time

from util90b import get_parser, to_dataset
from noniid_collision import collision_test
from partial_collection import partial_collection_test
from markov import markov_test
from maurer import maurer_universal_statistic
from frequency import frequency_test
from sanity_checks import compression_sanity_check, collision_sanity_check
from SP90Bv2_predictors import MultiMCW, Lag, MultiMMC, LZ78Y


##################
# main program
##################
if __name__ == '__main__':
    start_time = time.time()

    # get command line arguments
    args = get_parser('non-IID').parse_args()
    datafile = args.datafile
    bits_per_symbol = int(args.bits_per_symbol)
    if args.use_bits:
        use_bits = int(args.use_bits)
    else:
        use_bits = bits_per_symbol
    verbose = bool(args.verbose)

    with open(datafile, 'rb') as file:
        # Read in raw bytes and convert to list of output symbols
        bytes_in = bytearray(file.read())
        dataset = to_dataset(bytes_in, bits_per_symbol)
        if verbose:
            # print file and dataset details
            print ("Read in file %s, %d bytes long." % (datafile, len(bytes_in)))
            print ("Dataset: %d %d-bit symbols." % (len(dataset), bits_per_symbol))
            print ("Output symbol values: min = %d, max = %d" % (min(dataset), max(dataset)))

        # 
        if use_bits < bits_per_symbol:
            mask = (2**use_bits) - 1
            dataset = [s & mask for s in dataset]
            if verbose:
                print ("* Using only low %d bits out of %d." % (use_bits, bits_per_symbol))
                print ("* Using output symbol values: min = %d, max = %d" % (min(dataset), max(dataset)))

        n = 2 ** use_bits
        minEntropy = float(use_bits)

        # Section 6.3.2 The Collision Estimate
        pmax, minH = collision_test(dataset, n)
        if verbose:
            print("- Collision test          : p(max) = %g, min-entropy = %g" % (pmax, minH))
        minEntropy = min(minH, minEntropy)

        # Section 6.3.3 The Markov Estimate
        # If more than 6 bits per symbol, map down to 6 bits per symbol and run Markov test
        if use_bits > 6:
            pmax, minH = markov_test([s&63 for s in dataset], 64, 0.95)
            if verbose:
                print("- Markov test (map 6 bits): p(max) = %g, min-entropy = %g" % (pmax, minH))
        else:
            pmax, minH = markov_test(dataset, n, 0.99)
            if verbose:
                print("- Markov test             : p(max) = %g, min-entropy = %g" % (pmax, minH))
        minEntropy = min(minH, minEntropy)

        # Section 9.3.6 The Compression Test
        valid, pmax, minH = maurer_universal_statistic(dataset, n)
        if valid:
            if verbose:
                print("- Compression test        : p(max) = %g, min-entropy = %g" % (pmax, minH))
            minEntropy = min(minH, minEntropy)
        else:
            print("- Compression test *not valid* for this data set.")

        # Section 9.3.7 The Frequency Test
        pmax, minH = frequency_test(dataset, n, 0.95)
        if verbose:
            print("- Frequency test          : p(max) = %g, min-entropy = %g" % (pmax, minH))
        minEntropy = min(minH, minEntropy)

        # frequency test can give a negative result in extreme cases, so do following:
        minEntropy = max(0.0, minEntropy)

        
        pmax, minH = LZ78Y(dataset, verbose)
        if verbose:
            print("- LZ78Y predictor: p(max) = %g, min-entropy = %g" % (pmax, minH))
        minEntropy = min(minH, minEntropy)
        print("min-entropy = %g" % (minEntropy))

        # run_sanity_check:
        if minEntropy > 0.0:
            if compression_sanity_check(dataset, minEntropy, verbose)[0] == 'pass' and\
                collision_sanity_check(dataset, minEntropy, verbose) == 'pass':
                print("sanity check = PASS")
            else:
                print("sanity check = FAIL")
        else:
            # if min-entropy is 0.0 then no sanity check is needed
            print("sanity check = N/A")

    if verbose:
        print("time: (%g sec)" % (time.time() - start_time))
