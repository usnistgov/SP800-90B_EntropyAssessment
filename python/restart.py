# DRAFT NIST SP 800-90B (January 2016) Section 3.1.4 
#
# Sanity checks on restart datasets
#
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# Kerry McKay
# March 3, 2016

import sys
import time
import math

from util90b import get_parser, to_dataset, get_z
from mostCommonValue import most_common_restart


# The input file is the row dataset. This program will derive the colomn dataset from the row dataset.


##################
# main program
##################
if __name__ == '__main__':

    # get command line arguments
    args = get_parser('restart').parse_args()
    datafile = args.datafile
    bits_per_symbol = int(args.bits_per_symbol)
    verbose = bool(args.verbose)
    H_I = float(args.H_I)

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

        # max min-entropy is -log_2(1/k)
        minEntropy = float(math.log(k,2))


        # make sure that the dataset is the expected size
        # (otherwise, the column dataset won't be created correctly)
        assert len(dataset) == 1000000

    
        # 1. Let alpha be 0.01/(2000k), where k is the sample size
        alpha = 0.01/(2000*k)

        # 2. For each row of the matrix, find the frequency of the most
        #    common sample value. Let F_R be the maximum.
        if verbose:
            print ("\nRunning sanity check on row dataset:")
        f = [0 for i in range(1000)]
        for i in range(1000):
            f[i] = most_common_restart(dataset[i*1000:(i+1)*1000])
        F_R = max(f)
        if verbose:
            print ("- F_R: %d" % F_R)
            
        # 3. repeat the sample process for the column matrix. Let F_C be the
        #    maximum.
        if verbose:
            print ("Running sanity check on column dataset:")
        f = [0 for i in range(1000)]
        for i in range(1000):
            column = []
            for j in range(1000):
                column.append(dataset[j*1000+i])
            f[i] = most_common_restart(column)
        F_C = max(f)
        if verbose:
            print ("- F_C: %d" % F_C)

        # 4. Let F = max(F_R, F_C)
        F = max(F_R, F_C)

        # 5. Let p = 2**-H_I. Find the upper bound U of the (1-alpha)% confidence
        #    interval
        p = math.pow(2, -H_I)
        Z = get_z(alpha)
        U = 1000*p + Z*math.sqrt(1000*p*(1-p))
        if verbose: 
##            print ("alpha:", alpha )         
##            print ("z:",Z)
            print ("U: %f" % U)

        # 6. If F is greater than U, the test fails
        if F > U:
            print ("Failed the restart tests")
            print ("*** Validation failed. No entropy estimate awarded.")
        else:
            print ("Passed the restart tests")
            print ("*** Final entropy estimate: %f" % H_I)

        

