# DRAFT NIST SP 800-90B (January 2016) Section 5
#
# Estimating the Min-Entropy of IID Sources
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness.  See user guide for full disclaimer.
#
# Tim Hall
# tim.hall@nist.gov
#
# Kerry McKay
# kerry.mckay@nist.gov
#
# February 2016

import time
import sys


from util90b import get_parser, to_dataset
from permutation_tests import permutation_test
from chi_square_tests import pass_chi_square_tests
from mostCommonValue import most_common
from LRS import lenLRS


##################
# main program
##################
if __name__ == '__main__':
    start_time = time.time()

    # get command line arguments
    args = get_parser('IID').parse_args()
    datafile = args.datafile
    bits_per_symbol = int(args.bits_per_symbol)
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
            print ("Output symbol values: min = %d, max = %d\n" % (min(dataset), max(dataset)))

        # Run cheap and easy tests first.

        # run LRS test
        if lenLRS(dataset, verbose):
            if verbose:
                print ("** Passed LRS test")
        else:
            print ("** Failed LRS test")
            print ("IID = False")
            sys.exit(0)

        # run chi-square tests
        if pass_chi_square_tests(dataset, verbose):
            if verbose:
                print ("** Passed chi square tests")
        else:
            if verbose:
                print ("** Failed chi square tests")
            print ("IID = False")
            sys.exit(0)

        # Calculate min-entropy
        pmax, minH = most_common(dataset)
        print("min-entropy = %g" % (minH))

        # Run permutation test.
        # Run this last, because it is memory-intensive and time-consuming.
        if permutation_test(dataset, verbose):
            if verbose:
                print ("** Passed IID permutation tests")
        else:
            if verbose:
                print ("** Failed IID permutation tests")
            print ("IID = False")
            sys.exit(0)

        # If we made it this far, all tests have been passed.
        if verbose:
            print ("\nIID = True")

        print("\nDon't forget to run the sanity check on a restart dataset using H_I = %g" % minH )
