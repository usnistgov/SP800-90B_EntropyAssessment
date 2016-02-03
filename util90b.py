# DRAFT NIST SP 800-90B (August 2012) 
#
# Utility and helper functions for NIST SP 800-90B Testing
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness.
#
# Tim Hall
# tim.hall@nist.gov
#
# 28 August 2014

import argparse
import bisect

# add the command line option parsing details for either the
# IID or non-IID tests
# Defined here so as not to clutter up main program files
def get_parser(test):
    if test != 'IID':
        test = 'non-IID'

    descr = 'Run the Draft NIST SP 800-90B (August 2012) %s Tests' % test
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument(dest='datafile', metavar='datafile', help='dataset on which to run tests')
    parser.add_argument(dest='bits_per_symbol',metavar='bits_per_symbol', help='number of bits used to represent sample output values')

    if test == 'non-IID':
        parser.add_argument('-u', '--usebits', dest='use_bits', metavar='use_bits', help='use only the N lowest order bits per sample')
        
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='verbose mode: show detailed test results')

    return parser

# Convert a bytearray of raw bytes to a list containing the dataset
# Use bits-per-symbol to do the mapping
#
# Note that 1 bit per symbol assumes that 1 bit per byte is used, the data is
# not packed (i.e., 8 1-bit output samples are not packed into 1 byte)
#
# Does not handle > 32 bits per symbol.
##def to_dataset(bytes, bits_per_symbol):
##    assert bits_per_symbol > 0
##    assert bits_per_symbol <= 32
##
##    if bits_per_symbol <= 8: # includes 1-bit per symbol
##        return list(bytes)
##    elif bits_per_symbol <= 16:
##        return [bytes[i]*256 + bytes[i+1] for i in range(0,len(bytes),2)]
##    elif bits_per_symbol <= 24:
##        return [bytes[i]*(256*256) + bytes[i+1]*256 + bytes[i+2] for i in range(0,len(bytes),3)]
##    elif bits_per_symbol <= 32:
##        return [bytes[i]*(256*256*256) + bytes[i+1]*(256*256) + bytes[i+2]*256 + bytes[i+3] for i in range(0,len(bytes),4)]
##    else:
##        return list()

def to_dataset(bytes, bits_per_symbol):
    assert bits_per_symbol > 0
    assert bits_per_symbol <= 32

    # mask for relevant bits
    mask = 2**bits_per_symbol - 1

    # Ignore incomplete symbols. E.g., if there are 17 bytes of 16-bit symbols,
    # only read first 16 bytes and disgard remaining 4 bits.
    N = len(bytes)
    print "reading %d bytes of data" % N

    if bits_per_symbol <= 8: # includes 1-bit per symbol
##        return list(bytes)
        return [b & mask for b in bytes]
    elif bits_per_symbol <= 16:
        return [(bytes[i]*256 + bytes[i+1]) & mask for i in range(0,N,2)]
    elif bits_per_symbol <= 24:
        return [(bytes[i]*(256*256) + bytes[i+1]*256 + bytes[i+2])  & mask for i in range(0,N,3)]
    elif bits_per_symbol <= 32:
        return [(bytes[i]*(256*256*256) + bytes[i+1]*(256*256) + bytes[i+2]*256 + bytes[i+3]) & mask for i in range(0,N,4)]
    else:
        return list()

# Rank(S) function in Section 9.1.2
# S - score on original, unshuffled data set
# L - list of scores on 1,000 shuffled data sets
#
# NOTE: I think the Rank(S) function in the DRAFT SP 800-90B (Aug 2012)
# needs some work.  A few inconsistencies. -- TAH
def Rank(S, L):
    length = len(L)
    midpoint = length // 2
    
    if S == L[midpoint]:
        rank = midpoint
    elif S > L[midpoint]:
        rank = bisect.bisect_left(L, S)
    else:
        rank = bisect.bisect_right(L, S)

    return min(rank + 1, length)


# Map dataset to list 0..n. This is necessary for datasets here k is not
# a power of 2 or not all symbols in the range are present in the dataset.
# For example, a trinary source that outputs 13,14, and 15 would then
# convert the sequence to 0, 1, and 2 so that it will play nicely with
# the list indices used in the methods.
def mapData(seq):
    sortedData = sorted(list(set(seq))) #create new sorted list of unique symbols
    return [sortedData.index(seq[i]) for i in range(len(seq))]


