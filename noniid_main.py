# DRAFT NIST SP 800-90B (January 2016) Sections 6.2 and 6.3 
#
# Estimating the Min-Entropy of non-IID Sources
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# Oriniginal code by T. A. Hall
# tim.hall@nist.gov
#
# Updated by Kerry McKay
# March 3, 2016


import sys
import time
import math

from util90b import get_parser, to_dataset, mapData
from mostCommonValue import most_common
from noniid_collision import collision_test
from markov import markov_test
from maurer import maurer_universal_statistic
from SP90Bv2_predictors import MultiMCW, Lag, MultiMMC, LZ78Y
from tuple import t_tuple
from LRS import LRS_estimate


##################
# main program
##################
if __name__ == '__main__':

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
        k = len(set(dataset))
        if verbose:
            # print file and dataset details
            print ("Read in file %s, %d bytes long." % (datafile, len(bytes_in)))
            print ("Dataset: %d %d-bit symbols, %d symbols in alphabet." % (len(dataset), bits_per_symbol, k))
            print ("Output symbol values: min = %d, max = %d" % (min(dataset), max(dataset)))

        # 
        if use_bits < bits_per_symbol:
            mask = (2**use_bits) - 1
            dataset = [s & mask for s in dataset]
            k = len(set(dataset))
            if verbose:
                print ("* Using only low %d bits out of %d. %d symbols in reduced alphabet." % (use_bits, bits_per_symbol, k))
                print ("* Using output symbol values: min = %d, max = %d" % (min(dataset), max(dataset)))

        # some estimates require mapping the dataset to sequential outputs starting at 0
        mapped = mapData(dataset)

        # max min-entropy is -log_2(1/k)
        minEntropy = float(math.log(k,2))


        if verbose:
            print ("Running entropic statistic estimates:")
            
        # Section 6.3.1 The Most Common Value Estimate
        startTime = time.time()
        pmax, minH = most_common(dataset)
        endTime = time.time()
        if verbose:
            print("- Most Common Value Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, endTime - startTime))
        minEntropy = min(minH, minEntropy)

        # Section 6.3.2 The Collision Estimate
        startTime = time.time()
        pmax, minH = collision_test(dataset, k)
        endTime = time.time()
        if verbose:
            print("- Collision Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)

        # Section 6.3.3 The Markov Estimate
        # If more than 6 bits per symbol, map down to 6 bits per symbol and run Markov test
        if use_bits > 6:
            startTime = time.time()
            pmax, minH = markov_test([s&63 for s in dataset], 64, 0.99)
            endTime = time.time()
            if verbose:
                print("- Markov Estimate (map 6 bits): p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        else:
            startTime = time.time()
            pmax, minH = markov_test(mapped, k, 0.99)
            endTime = time.time()
            if verbose:
                print("- Markov Estimate: p(max) = %g, min-entropy = %g  (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)

        # Section 6.3.4 The Compression Estimate
        startTime = time.time()
        pmax, minH = maurer_universal_statistic(mapped, k, verbose)
        endTime = time.time()
        if verbose:
            print("- Compression Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)

        # Section 6.3.5 The t-Tuple Estimate
        startTime = time.time()
        pmax, minH = t_tuple(dataset, verbose)
        endTime = time.time()
        if verbose:
            print("- t-Tuple Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)
        
        # Section 6.3.6 The LRS Estimate
        startTime = time.time()
        pmax, minH = LRS_estimate(dataset, verbose)
        endTime = time.time()
        if verbose:
            print("- LRS Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)


        # Section 6.3.7 Multi Most Common in Window prediction estimate
        startTime = time.time()
        pmax, minH = MultiMCW(dataset, verbose)
        endTime = time.time()
        if verbose:
            print("MultiMCW Prediction Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)

        # Section 6.3.8 Lag prediction estimate
        startTime = time.time()
        pmax, minH = Lag(dataset, verbose)
        endTime = time.time()
        if verbose:
            print("Lag Prediction Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)


        # Section 6.3.9 MultiMMC prediction estimate
        startTime = time.time()
        pmax, minH = MultiMMC(dataset, verbose)
        endTime = time.time()
        if verbose:
            print("MultiMMC Prediction Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)
        

        # Section 6.3.10 LZ78Y prediction estimate
        startTime = time.time()
        pmax, minH = LZ78Y(dataset, verbose)
        endTime = time.time()
        if verbose:
            print("LZ78Y Prediction Estimate: p(max) = %g, min-entropy = %g (%f s)" % (pmax, minH, (endTime - startTime)))
        minEntropy = min(minH, minEntropy)

        print("-----------------------")
        print("min-entropy = %g" % (minEntropy))


        print("\nDon't forget to run the sanity check on a restart dataset using H_I = %g" % minEntropy )
