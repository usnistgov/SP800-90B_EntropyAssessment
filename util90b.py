# DRAFT NIST SP 800-90B (January 2016) 
#
# Utility and helper functions for NIST SP 800-90B Testing
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness.
#
# 2012 code by Tim Hall
# tim.hall@nist.gov
#
# Updated for January 2016 draft by Kerry McKay
# kerry.mckay@nist.gov

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
# Does not handle > 32 bits per symbol.
def to_dataset(bytes, bits_per_symbol):
    assert bits_per_symbol > 0
    assert bits_per_symbol <= 8

    # mask for relevant bits
    mask = 2**bits_per_symbol - 1

    # Ignore incomplete symbols. E.g., if there are 17 bytes of 16-bit symbols,
    # only read first 16 bytes and disgard remaining 4 bits.
    N = len(bytes)
    print "reading %d bytes of data" % N

    if bits_per_symbol <= 8: # includes 1-bit per symbol
        return [b & mask for b in bytes]
##    elif bits_per_symbol <= 16:
##        return [(bytes[i]*256 + bytes[i+1]) & mask for i in range(0,N,2)]
##    elif bits_per_symbol <= 24:
##        return [(bytes[i]*(256*256) + bytes[i+1]*256 + bytes[i+2])  & mask for i in range(0,N,3)]
##    elif bits_per_symbol <= 32:
##        return [(bytes[i]*(256*256*256) + bytes[i+1]*(256*256) + bytes[i+2]*256 + bytes[i+3]) & mask for i in range(0,N,4)]
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
# alpha_min= 0.00000002, alpha_max=0.0000025
def get_z(alpha):
    z_values = {0.0000250: 4.21479967,
                0.0000245: 4.219357057,
                0.0000240: 4.2240038,
                0.0000235: 4.228743578,
                0.0000230: 4.233580304,
                0.0000225: 4.238518145,
                0.0000220: 4.243561545,
                0.0000215: 4.248715249,
                0.0000210: 4.253984334,
                0.0000205: 4.259374237,
                0.0000200: 4.264890794,
                0.0000195: 4.270540276,
                0.0000190: 4.276329438,
                0.0000185: 4.282265566,
                0.0000180: 4.288356537,
                0.0000175: 4.294610884,
                0.0000170: 4.301037872,
                0.0000165: 4.307647584,
                0.0000160: 4.314451022,
                0.0000155: 4.321460219,
                0.0000150: 4.328688376,
                0.0000145: 4.33615002,
                0.0000140: 4.343861183,
                0.0000135: 4.351839624,
                0.0000130: 4.360105083,
                0.0000125: 4.368679593,
                0.0000120: 4.377587847,
                0.0000115: 4.386857647,
                0.0000110: 4.396520453,
                0.0000105: 4.406612056,
                0.0000100: 4.417173413,
                0.0000095: 4.428251702,
                0.0000090: 4.439901645,
                0.0000085: 4.452187228,
                0.0000080: 4.465183916,
                0.0000075: 4.478981593,
                0.0000070: 4.493688505,
                0.0000065: 4.509436651,
                0.0000060: 4.526389321,
                0.0000055: 4.544751891,
                0.0000050: 4.56478773,
                0.0000045: 4.586842452,
                0.0000040: 4.611382362,
                0.0000035: 4.639058488,
                0.0000030: 4.670819821,
                0.0000025: 4.708129716,
                0.0000020: 4.753424309,
                0.0000015: 4.811235695,
                0.0000010: 4.891638475,
                0.0000005: 5.026312835,
                0.0000002: 5.61610279}
    
    if alpha in z_values.keys():
        return z_values[alpha]
    else:
        print "need to find alpha"
        if alpha < min(z_values.keys()):
            return z_values[min(z_values.keys())]
