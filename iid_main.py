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

import os
import time
import sys


from util90b import get_parser, to_dataset_packed, to_dataset_unpacked
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
    start_position = int (args.start_position)
    read_amount = int (args.read_amount)
    packed = bool (args.packed)

    with open(datafile, 'rb') as file:
        # Read in raw bytes and convert to list of output symbols
        file.seek(start_position, os.SEEK_SET)
        bytes_in = bytearray(file.read(read_amount))
        if packed:
          dataset = to_dataset_packed (bytes_in, bits_per_symbol)
        else:
          dataset = to_dataset_unpacked (bytes_in, bits_per_symbol)
        k = len(set(dataset))
        if verbose:
            # print file and dataset details
            print ("Read in file %s, %d bytes long." % (datafile, len(bytes_in)))
            print ("Dataset: %d %d-bit symbols, %d symbols in alphabet." % (len(dataset), bits_per_symbol, k))
            print ("Output symbol values: min = %d, max = %d\n" % (min(dataset), max(dataset)))

        #######################################
        # STEP 1: Determine if Dataset is IID #
        #######################################
        # determine if dataset is IID using shuffle and Chi-square tests
        passed_permutation_tests = permutation_test(list(dataset), verbose)

        if passed_permutation_tests:
            if verbose:
                print ("** Passed IID permutation tests")
        else:
            if verbose:
                print ("** Failed IID permutation tests")
            print ("IID = False")
            sys.exit(0)

        # run chi-square tests on dataset
        if pass_chi_square_tests(dataset, verbose):
            if verbose:
                print ("** Passed chi square tests")
        else:
            if verbose:
                print ("** Failed chi square tests")
            print ("IID = False")
            sys.exit(0)

        # run LRS test
        if lenLRS(dataset, verbose):
            if verbose:
                print ("** Passed LRS test")
            print ("\nIID = True")
        else:
            print ("** Failed LRS test")
            print ("IID = False")
            sys.exit(0)

        ############################################
        # STEP 2: Calculate min-entropy of dataset #
        ############################################
        pmax, minH = most_common(dataset)
        print("min-entropy = %g" % (minH))


        print("\nDon't forget to run the sanity check on a restart dataset using H_I = %g" % minH )
