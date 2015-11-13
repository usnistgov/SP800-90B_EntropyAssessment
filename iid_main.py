# DRAFT NIST SP 800-90B (August 2012) Section 9.1 and 9.2 tests
#
# Estimating the Min-Entropy of IID Sources
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness.
#
# Tim Hall
# tim.hall@nist.gov
#
# 3 September 2014

import time
import sys


from util90b import get_parser, to_dataset
from shuffle_tests import shuffle_test_scores, get_shuffle_tests
from chi_square_tests import pass_chi_square_tests
from sanity_checks import compression_sanity_check, collision_sanity_check
from mostCommonValue import most_common


##################
# main program
##################
if __name__ == '__main__':
    start_time = time.time()

    # get command line arguments
    args = get_parser('IID').parse_args()
    datafile = args.datafile
    bits_per_symbol = int(args.bits_per_symbol)
    number_of_shuffles = int(args.number_of_shuffles)
    verbose = bool(args.verbose)

    with open(datafile, 'rb') as file:
        # Read in raw bytes and convert to list of output symbols
        bytes_in = bytearray(file.read())
        dataset = to_dataset(bytes_in, bits_per_symbol)
        if verbose:
            # print file and dataset details
            print ("Read in file %s, %d bytes long." % (datafile, len(bytes_in)))
            print ("Dataset: %d %d-bit symbols." % (len(dataset), bits_per_symbol))
            print ("Output symbol values: min = %d, max = %d\n" % (min(dataset), max(dataset)))

        #######################################
        # STEP 1: Determine if Dataset is IID #
        #######################################
        # determine if dataset is IID using shuffle and Chi-square tests
        # divide the dataset into 10 subsets of equal length
        sublength = len(dataset) // 10
        data_subsets = [dataset[i*sublength : (i+1)*sublength] for i in range(10)]

        ranks = [shuffle_test_scores((data_subsets[i], number_of_shuffles)) for i in range(10)]

        # for each shuffle test, find rank of original subsets in shuffled subsets
        # e.g., for 1,000 shuffles, the limits are 50 and 950
        llimit = number_of_shuffles // 20
        ulimit = number_of_shuffles - llimit
        passed_shuffle_tests = True

        for test in get_shuffle_tests():
            if verbose:
                # print column headers: test name, scores on unshuffled dataset
                # and ranks of scores against scores on shuffled datasets
                sys.stdout.write("\n{0} Test:".format(test.name.capitalize()).ljust(24))
                print("\n\tScores".ljust(24) + "   Ranks")

            bad_rank_count = [0 for s in test.scores]

            # for each of the 10 data subsets
            for i in range(10):
                ranks_i = [r for r in ranks[i][test.name][1]]

                if verbose:
                    # print out scores and ranks for this test
                    scores = (test.format % ranks[i][test.name][0]).ljust(24)
                    # identify ranks in lowest or highest 5% with asterisk (*)
                    rfmt = lambda r, l, u: "%4d "%r if r > l and r < u else "%4d*"%r
                    rankings = [rfmt(r,llimit,ulimit) for r in ranks_i]
                    rankings = ' '.join(rankings)
                    print("\t" + scores + rankings)

                for j in range(len(ranks_i)):
                    if ranks_i[j] <= llimit or ranks_i[j] >= ulimit:
                        bad_rank_count[j] += 1

            if verbose:
                # print out the ranking results
                rr = [" %4d"%brc for brc in bad_rank_count]
                print("\t ".ljust(24) + ("  --- ")*len(rr))
                print("\t ".ljust(24) + " ".join(rr))

            # 8 or more out of 10 ranks in top or bottom 5% for a given
            # score means that subset fails iid
            passed_this_test = True
            for brc in bad_rank_count:
                if brc > 7:
                    passed_this_test = False
                    passed_shuffle_tests = False

            if verbose:
                if passed_this_test:
                    print("Passed %s Test" % test.name.capitalize())
                else:
                    print("Failed %s Test" % test.name.capitalize())

        if passed_shuffle_tests:
            if verbose:
                print("\n** Passed iid shuffle tests")
                print()
        else:
            if verbose:
                print("\n** Failed iid shuffle tests")
            print("IID = False")
            sys.exit(0)

        # run chi-square tests on dataset
        if pass_chi_square_tests(dataset, data_subsets, verbose):
            print("IID = True")
        else:
            print("IID = False")
            sys.exit(0)

        ############################################
        # STEP 2: Calculate min-entropy of dataset #
        ############################################
        minH = most_common(dataset)
        print("min-entropy = %g" % (minH))

        #######################################
        # STEP 3: Run Sanity Check on dataset #
        #######################################
        if compression_sanity_check(dataset, minH, verbose)[0] == 'pass' and\
            collision_sanity_check(dataset, minH, verbose) == 'pass':
            print("sanity check = PASS")
        else:
            print("sanity check = FAIL")

    if verbose:
        print("time: (%g sec)" % (time.time() - start_time))
