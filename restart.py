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

from collections import Counter

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
    start_position = int (args.start_position)
    read_amount = int (args.read_amount)

    with open(datafile, 'rb') as file:
        # Read in raw bytes and convert to list of output symbols
        bytes_in = bytearray(file.read(start_position))
        bytes_in = bytearray(file.read(read_amount))
        dataset = to_dataset(bytes_in, bits_per_symbol)
        k = len(set(dataset))
        if verbose:
            # print file and dataset details
            print ("Read in file %s, %d bytes long." % (datafile, len(bytes_in)))
            print ("Dataset: %d %d-bit symbols, %d symbols in alphabet." % (len(dataset), bits_per_symbol, k))
            print ("Output symbol values: min = %d, max = %d" % (min(dataset), max(dataset)))
            for cti in range (10, 1000000, 50000):
              print ("sample dataset value %d" % (dataset [cti]))
            tcount = Counter (dataset)
            tcount = tcount.most_common (10)
            for pair in tcount:
              print ("val= %d  count= %d" % (pair [0], pair [1]))

        # max min-entropy is -log_2(1/k)
        minEntropy = float(math.log(k,2))


        # make sure that the dataset is the expected size
        # (otherwise, the column dataset won't be created correctly)
        dataset = dataset [:1000000]
        assert len(dataset) == 1000000

    
        # 1. Let alpha be 0.01/(2000k), where k is the sample size
        alpha = 0.01/(2000*k)

        # 2. For each row of the matrix, find the frequency of the most
        #    common sample value. Let F_R be the maximum.
        #    This is an inefficient way to determine the max frequency,
        #    but it creates and preserves the list of row frequecies
        #    that can be used for future analysis, so leave it.  
        #    Also, optimization yields little benefit here.
        if verbose:
            print ("\nRunning sanity check on row dataset:")
        fr = [0] * 1000
        maxf = 1
        for i in range(1000):
            fr [i] = most_common_restart(dataset[i*1000:(i+1)*1000])
            if fr [i] >= maxf:
              maxf = fr [i]
              if verbose:
                print ("i %d  max freq %d" % (i, maxf))
        F_R = max (fr)
        if verbose:
            print ("- F_R: %d" % F_R)
            
        # 3. repeat the sample process for the column matrix. Let F_C be the
        #    maximum.
        #    This is an inefficient way to determine the max frequency,
        #    but it creates and preserves the list of column frequecies
        #    that can be used for future analysis, so leave it.  
        #    Also, optimization yields little benefit here.
        if verbose:
            print ("Running sanity check on column dataset:")
        fc = [0] * 1000
        maxf = 1
        for i in range(1000):
            column = []
            for j in range(1000):
                column.append(dataset[j*1000+i])
            fc [i] = most_common_restart(column)
            if fc [i] >= maxf:
              maxf = fc [i]
              if verbose:
                print ("i %d  max freq %d" % (i, maxf))
        F_C = max (fc)
        if verbose:
            print ("- F_C: %d" % F_C)

        # 4. Let F = max(F_R, F_C)
        F = max(F_R, F_C)

        # 5. Let p = 2**-H_I. Find the upper bound U of the (1-alpha)% confidence
        #    interval
        p = math.pow (2, -H_I)
        Z = get_z (alpha)
        U = (1000 * p) + (Z * math.sqrt (1000 * p * (1-p)))
        if verbose: 
            print ("    p: %8.5e" % (p))         
            print ("alpha: %8.5e"% (alpha))         
            print ("    Z: %8.5e" % (Z))
            print ("    U: %f" % (U))

        # 6. If F is greater than U, the test fails
        if F > U:
            print ("Failed the restart tests")
            print ("*** Validation failed. No entropy estimate awarded.")
        else:
            print ("Passed the restart tests")
            print ("*** Final entropy estimate: %f" % H_I)

        

