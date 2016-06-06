# DRAFT NIST SP 800-90B (January 2016) 
#
# Utility and helper functions for NIST SP 800-90B Testing
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# 2012 code by Tim Hall
# tim.hall@nist.gov
#
# Updated for January 2016 draft by Kerry McKay
# June 6, 2016

import argparse
import bisect

# add the command line option parsing details for either the
# IID, non-IID, or restart tests
# Defined here so as not to clutter up main program files
def get_parser(test):
    if test != 'IID' and test != 'restart':
        test = 'non-IID'

    descr = 'Run the Draft NIST SP 800-90B (January 2016) %s Tests' % test
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument(dest='datafile', metavar='datafile', help='dataset on which to run tests')
    parser.add_argument(dest='bits_per_symbol',metavar='bits_per_symbol', help='number of bits used to represent sample output values')

    if test == 'non-IID':
        parser.add_argument('-u', '--usebits', dest='use_bits', metavar='use_bits', help='use only the N lowest order bits per sample')
    elif test == 'restart':
        parser.add_argument(dest='H_I', metavar='H_I', help='initial entropy estimate')
        
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='verbose mode: show detailed test results')

    return parser

# Convert a bytearray of raw bytes to a list containing the dataset
# Use bits-per-symbol to do the mapping
#
# Note that 1 bit per symbol assumes that 1 bit per byte is used, the data is
# not packed (i.e., 8 1-bit output samples are not packed into 1 byte)
#
# Does not handle > 8 bits per symbol.
def to_dataset(bytes, bits_per_symbol):
    assert bits_per_symbol > 0
    assert bits_per_symbol <= 8

    # mask for relevant bits
    mask = 2**bits_per_symbol - 1

    N = len(bytes)
    print ("reading %d bytes of data" % N)

    if bits_per_symbol <= 8: # includes 1-bit per symbol
        return [b & mask for b in bytes]
    else:
        return list()


# Map dataset to list 0..n. This is necessary for datasets here k is not
# a power of 2 or not all symbols in the range are present in the dataset.
# For example, a trinary source that outputs 13,14, and 15 would then
# convert the sequence to 0, 1, and 2 so that it will play nicely with
# the list indices used in the methods.
def mapData(seq):
    sortedData = sorted(list(set(seq))) #create new sorted list of unique symbols
    return [sortedData.index(seq[i]) for i in range(len(seq))]



# Find the z-value for restart tests
# Given alpha, find z corresponding to upper bound of confidence interval
# alpha_min= 1.95312E-08, alpha_max=2.5E-06
def get_z(alpha):
    z_values = {2.5E-06: 4.708129716,
                2.45E-06: 4.712247134,
                2.4E-06: 4.716446023,
                2.35E-06: 4.72072975,
                2.3E-06: 4.725101895,
                2.25E-06: 4.729566273,
                2.2E-06: 4.734126952,
                2.15E-06: 4.738788277,
                2.1E-06: 4.743554897,
                2.05E-06: 4.748431792,
                2E-06: 4.753424309,
                1.95E-06: 4.758538192,
                1.9E-06: 4.763779631,
                1.85E-06: 4.769155304,
                1.8E-06: 4.77467243,
                1.75E-06: 4.780338834,
                1.7E-06: 4.78616301,
                1.65E-06: 4.792154206,
                1.6E-06: 4.798322513,
                1.55E-06: 4.804678972,
                1.5E-06: 4.811235695,
                1.45E-06: 4.818006016,
                1.4E-06: 4.825004652,
                1.35E-06: 4.832247908,
                1.3E-06: 4.839753915,
                1.25E-06: 4.847542912,
                1.2E-06: 4.855637585,
                1.15E-06: 4.864063487,
                1.1E-06: 4.872849536,
                1.05E-06: 4.882028643,
                1E-06: 4.891638476,
                9.5E-07: 4.901722437,
                9E-07: 4.912330888,
                8.5E-07: 4.923522727,
                8E-07: 4.935367445,
                7.5E-07: 4.947947838,
                7E-07: 4.96136365,
                6.5E-07: 4.975736557,
                6E-07: 4.99121714,
                5.5E-07: 5.007994871,
                5E-07: 5.026312837,
                4.5E-07: 5.046490185,
                4E-07: 5.06895775,
                3.5E-07: 5.094317391,
                3E-07: 5.123447001,
                2.5E-07: 5.157701313,
                2E-07: 5.199337583,
                1.5E-07: 5.252559482,
                1E-07: 5.326723889,
                5E-08: 5.451310443,
                1.95312E-08: 5.61610279}
    
    if alpha in z_values.keys():
        return z_values[alpha]
    else:
        if alpha < min(z_values.keys()):
            return z_values[min(z_values.keys())]
        elif alpha > max(z_values.keys()):
            return z_values[max(z_values.keys())]
        else:
            #find closest match
            best = 10 #just a number that's too high
            for k in z_values.keys():
                if abs(k-alpha) < abs(best-alpha):
                    best = k
            return z_values[best]
                
