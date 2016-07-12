
# DRAFT NIST SP 800-90B (January 2016)
#
# Section 5.2 - Additional Chi-square Statistical Tests
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# Tim Hall
# tim.hall@nist.gov
# 3 September 2014
#
# Updated by Kerry McKay
# March 3, 2106

import math
from collections import OrderedDict, Counter
from operator import itemgetter
import itertools

# does the dataset pass the chi-square tests?
def pass_chi_square_tests(dataset, verbose=False):

    # chi-square independence test
    is_binary = (min(dataset) == 0 and max(dataset) == 1)

    if is_binary:
        score, df = binary_chi_square_independence(dataset)
    else:
        score, df = chi_square_independence(dataset)
    cutoff = chi_square_cutoff(df)
    if verbose:
        print("\nChi square independence\n\tscore = %g, degrees of freedom = %d, cut-off = %g" % (score, df, cutoff))

    if score < cutoff:
        if verbose:
            print("** Passed chi-square independence test")
    else:
        if verbose:
            print("** Failed chi-square independence tests")
            return False

    # chi-square goodness-of-fit test
    # divide the dataset into 10 subsets of equal length
    sublength = len(dataset) // 10
    data_subsets = [dataset[i*sublength : (i+1)*sublength] for i in range(10)]
    if is_binary:
        score, df = binary_goodness_of_fit(data_subsets)
    else:
        score, df = goodness_of_fit(data_subsets)
    cutoff = chi_square_cutoff(df)
    if verbose:
        print("\nChi square goodness-of-fit\n\tscore = %g, degrees of freedom = %d cut-off = %g" % (score, df, cutoff))

    if score < cutoff:
        if verbose:
            print("** Passed chi-square goodness-of-fit test\n")
    else:
        if verbose:
            print("** Failed chi-square goodness-of-fit tests")
            return False

    return True


# helper function - note that it expects a sequence (e.g., list) of
# datasets.  When passing in entire dataset, pass in as one-element list
def _internal_get_symbol_counts(datasets):
    # Construct a count of every value that appears in the datasets.
    values = {}
    count = 0
    for xl in datasets:
        for x in xl:
            count +=1
            values[x] = values.get(x,0)+1
            
    return values, count


# Testing Independence for Non-Binary Data - Section 5.2.1
def chi_square_independence(s):
    A = set(s)
    L = len(s)

    p = {} #proportions
    e = {} #expected values
    pair_counts = {} #occurances of each pair
    
    # 1. Find the proportion p_i of each x_i in S. Calculate the expected
    #    number of occurances of each possible pair as (p_i)(p_j)(L-1)
    for a in A:
        p[a]= float(s.count(a))/L
    pair_iterator = itertools.product(A,A)
    more = True
    while more:
        try:
            pair = pair_iterator.next()
            e[pair] = p[pair[0]]*p[pair[1]]*(L-1)
        except:
            more = False

    for i in range(L-1):
        pair = tuple(s[i:i+2])
        pair_counts[pair] =  pair_counts.get(pair,0)+1

    # 2. Allocate the possible pairs, starting from the smallest expected value, 
    #    into bins s.t. the expected value of each bin is at least 5. The
    #    expected value of a bin is equal to the sum of the expected values of
    #    the pairs that are included in the bin. After allocating all pairs, if
    #    the expected value of the last bin is less than 5, merge the last two
    #    bins. Let q be the number of bins constructed after this procedure.
    #
    #    For each pair, count the number of observed values for each bin
    
    bins = []
    e_sorted = OrderedDict(sorted(e.items(), key=itemgetter(1)))
    q = 0
    
    # bins[q][0] is pairs, bins[q][1] is frequency, bins[q][2] is expected
    bins.append(['',0,0])
    
    for pair in e_sorted:
        bins[q][0] = bins[q][0] + str(pair)
        bins[q][1] = bins[q][1] + pair_counts.get(pair,0)
        bins[q][2] = bins[q][2] + e[pair]

        if e[pair] >= 5:
            q += 1
            bins.append(['',0,0])

    if bins[q][0] == '':            
        bins.pop() #remove empty bin
        q -= 1

    #merge last two bins if expected value of last bin < 5
    if bins[q][2] < 5:
        bins[q-1][0] = bins[q-1][0] + bins[q][0]
        bins[q-1][1] = bins[q-1][1] + bins[q][1]
        bins[q-1][2] = bins[q-1][2] + bins[q][2]
        bins.pop()
        q -= 1

    # caclulate the test statistic, T
    T = 0
    for i in range(q):
        T += float((bins[i][1] - bins[i][2])**2)/bins[i][2]

    # return statistic with q-1 df (since our indices start at 0, df is q)
    return T, q



# HELPER FUNCTIONS

def hamming_weight(x):
    return bin(x).count('1')

def kbit_element(bits):
    k = len(bits)
    return sum([bits[i] << (k-i-1) for i in range(k)])


# Testing independence for Binary Data - Section 5.2.3
def binary_chi_square_independence(s):
    L = len(s)

    # 1. Let p0 and p1 be the proportion of zeroes and ones in S
    p0 = s.count(0)/float(L)
    p1 = s.count(1)/float(L)

    # 2. Find the maximum integer m such that p0**m > 5/L and p1**m < 5/L.
    #    If m is greater than 11, assign m to 11. If m is 1, the test fails.
    m = 11
    while m > 1:
        if math.pow(p0,m)>(5.0/L) and math.pow(p1,m)>(5.0/L):
            break
        else:
            m -= 1

    # Apply the test if m >= 2
    if m > 1:
        # Initialize T to 0.
        T = 0

        # 2. for each possible m-bit tuple
        # 2a. Let o be the number of times that the pattern occurs in the
        #     input sequence S. Note that the tuples are allowed to overlap.
        m_tuples = itertools.product([0,1], repeat=m)
        o = {}
        for i in range(L-m+1):
            t = tuple(s[i:i+m])
            o[t] = o.get(t, 0) + 1

        for t in m_tuples:
            # 2b. Let w be the number of ones in the tuple
            w = sum(t)
            
            # 2c. Let e = p_1**w * p_0**(m-w) * (L-m+1)
            e = math.pow(p1,w)*math.pow(p0,m-w)*(L-m+1)

            #update T
            T += math.pow(o.get(tuple(t),0)-e,2)/float(e)
                
        return T, math.pow(2,m)-1
    
    else:
        return None, None



# Testing Goodness-of-fit for Binary Data - Section 5.2.4
# This test checks the distribution of the number of ones in non-overlapping
# intervals of the input data to determine whether the distribution of the
# ones remains the same throughout the sequence.
def binary_goodness_of_fit(subsets):
    assert len(subsets) == 10

    # 1. Let p be the proportion of ones in S, i.e., p = (the number of ones in S)/ L.
    c = [sum(s) for s in subsets] # number of ones in all 10 subsets
    L = sum([len(s) for s in subsets]) # number of elements in all 10 subsets
    p = sum(c)/float(L) #proportion

    
    # 2. Partition S into ten non-overlapping sub-sequences of length floor(L/10).
    #   (Data already split before callint this function)

    # 3.Initialize T to 0.
    T = 0

    # 4. Let the expected number of ones in each sub-sequence S_d be e=p*floor(L/10)
    e = p*math.floor(L/10)

    # 5. for d = 0 to 9
    for d in range(10):
        # 5a. Let o_i be the number of ones in S_d
        o = sum(subsets[d])

        # 5b. update the test statistic, T
        T += float(o - e)**2/e

    # T is a chi-square random variable with 9 degrees of freedom.
    return T, 9
    

# Test for Goodness of Fit for non-binary data - Section 5.2.2
def goodness_of_fit(subsets):
    assert len(subsets) == 10

    # 1. Let c_i be the number of occurrances of x_i in s, and let e_i = c_i/10.
    e, N = _internal_get_symbol_counts(subsets)
    for xi in e.keys():
        e[xi] /= 10.0

    # 2. Let List[i] be the sample value with the ith smallest e_i
    List = sorted(e.items())


    # 3. Starting from List[1], allocate the sample values into bins. Assign
    #    consecutive List[i] values to a bin until the sum of the e_i for those
    #    binned items is at least five, then begin assigning the following
    #    List[i] values to the next bin. If the expected value of the last bin
    #    is less than five, merge the last two bins. Let q be the number of bins
    #    constructed after this procedure.
    #
    # 4. Let E_i be the expected number of sample values in Bin i 
    E = [0]
    bins = [[]]
    q = 0
    
    for i in range(len(List)):
        bins[q].append(List[i][0])
        E[q] += List[i][1]
        
        if E[q] >= 5:
            q += 1
            E.append(0)
            bins.append([])

    if E[q] == 0:           
        #remove empty bin
        E.pop()
        q -= 1

    # The chi-square goodness-of-fit test is executed as follows:
    T = 0
    for d in range(10):
        o = [0 for i in range(q+1)]
        for i in range(q+1):
            # Let o_i be the number of sample values from bin i in the data subset s_d
            for x in bins[i]:
                o[i] += subsets[d].count(x)
            T += float(o[i]-E[i])**2/E[i]
            
    # return statistic and 9*(q-1) df
    # (since our bin indices start at 0 instead of 1, df=9q)
    return T, 9*q
    
   


################################################
# CHI-SQUARE CUTOFF VALUE (I.E., CRITICAL VALUE)
################################################

# Return the chi-square cutoff (or critical value) for the given degrees
# of freedom for p = 0.001.
def chi_square_cutoff(df):
    assert df > 0

    if df < 101:
        return critical_value[df]
    else:
        return calc_chi_square_cutoff(df)


# Table of Chi-square critical values.
#
# Values for nu (i.e., df) = 1 to 100 taken from:
# http://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm
#
critical_value =\
    {1:10.828, 2:13.816, 3:16.266, 4:18.467, 5:20.515, 6:22.458,\
    7:24.322, 8:26.125, 9:27.877, 10:29.588, 11:31.264, 12:32.91,\
    13:34.528, 14:36.123, 15:37.697, 16:39.252, 17:40.79, 18:42.312,\
    19:43.82, 20:45.315, 21:46.797, 22:48.268, 23:49.728, 24:51.179,\
    25:52.62, 26:54.052, 27:55.476, 28:56.892, 29:58.301, 30:59.703,\
    31:61.098, 32:62.487, 33:63.87, 34:65.247, 35:66.619, 36:67.985,\
    37:69.347, 38:70.703, 39:72.055, 40:73.402, 41:74.745, 42:76.084,\
    43:77.419, 44:78.75, 45:80.077, 46:81.4, 47:82.72, 48:84.037,\
    49:85.351, 50:86.661, 51:87.968, 52:89.272, 53:90.573, 54:91.872,\
    55:93.168, 56:94.461, 57:95.751, 58:97.039, 59:98.324, 60:99.607,\
    61:100.888, 62:102.166, 63:103.442, 64:104.716, 65:105.988, 66:107.258,\
    67:108.526, 68:109.791, 69:111.055, 70:112.317, 71:113.577, 72:114.835,\
    73:116.092, 74:117.346, 75:118.599, 76:119.85, 77:121.1, 78:122.348,\
    79:123.594, 80:124.839, 81:126.083, 82:127.324, 83:128.565, 84:129.804,\
    85:131.041, 86:132.277, 87:133.512, 88:134.746, 89:135.978, 90:137.208,\
    91:138.438, 92:139.666, 93:140.893, 94:142.119, 95:143.344, 96:144.567,\
    97:145.789, 98:147.01, 99:148.23, 100:149.449}



# Use a modified Wilson-Hilferty transformation to approximate the
# chi-square cutoff (critical value)
#
# Refer to Abramowitz and Stegun, "Approximations to the Chi-Square
# Distribution for large v (nu)" and "Approximations for the Inverse
# Function for large v (nu)" both on p.941 in Chapter 26.
#
# http://people.math.sfu.ca/~cbm/aands/page_941.htm
#
# "Handbook of Mathematical Functions", M. Abramowitz and I. A. Stegun
#
def calc_chi_square_cutoff(df):

    # For use with degrees of freedom > 30
    assert df > 30

    # 1. x_p s.t. Q(x_p) = 1 - P(x_p) = p = 0.001
    #    P(x) is cdf for standard normal distribution
    #    critical value x_p = +3.090
    #    See http://www.itl.nist.gov/div898/handbook/eda/section3/eda3671.htm
    x_p = 3.090

    # 2. Compute h_v using Abramowitz and Stegun 26.4.15 (p. 941)
    #    A and S has h60 = 0.0043 for x = 3.0
    #    so I figure h60 is about 0.0048 for 3.09
    h60 = 0.0048
    h_v = (60.0/df) * h60

    # 3. Approximate using Abramowitz and Stegun 26.4.18 (p. 941)
    term = 2.0/(9.0 * df)
    chisqr = df * pow(1.0 - term + (x_p - h_v)*math.sqrt(term), 3)

    return chisqr


