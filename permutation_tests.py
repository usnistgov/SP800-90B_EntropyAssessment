# Permutation Tests on Independence
# DRAFT NIST SP 800-90B (January 2016) Section 5.1
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# T. A. Hall, 2012
# tim.hall@nist.gov
#
# Updated by Kerry McKay
# March 3, 2016


import bz2
import sys
from collections import namedtuple
from random import randrange


# Use numpy shuffle function if available
# Much faster
try:
    import numpy as np
    shuffle = np.random.shuffle
except:
    import random
    shuffle = random.shuffle
   

# Calculate the mean and median and identify if binary or not
def calc_stats(dataset):
    mean = sum(dataset)/float(len(dataset))

    if min(dataset) == 0 and max(dataset) == 1:
        is_binary = True
        median = 0.5
    else:
        is_binary = False
        sd = sorted(dataset)
        halfway = len(sd) // 2
        if len(sd) % 2 == 0:
            median = (sd[halfway] + sd[halfway-1]) / 2.0
        else:
            median = sd[halfway]

    return (mean, median, is_binary)



###########################################
# Section 5.1 conversions for binary data #
###########################################

# conversion I
# partitions the sequences into 8-bit non-overlapping blocks, and counts the number of ones in each block.
def conversion1(s):
    sp = [s[i:i+8].count(1) for i in range(0,int(len(s)/8),8)]

    # handle incomplete block
    if len(s)%8 > 0:
        sp.append(s[len(sp)*8:].count(1))
    
    return sp

# conversion II
# partitions the sequences into 8-bit non-overlapping blocks, and calculates the integer value of each block. 
def conversion2(s):
    sp = str(s)[1:-1]
    sp = sp.replace(', ','')
    padded = sp + '0' * (8 - len(s) % 8)
    sp_bytes = [int(padded[i:i+8], 2) for i in range(0,int(len(padded))/8,8)]
    
    return sp_bytes


###############################
# Section 5.1 test statistics #
###############################

# Section 5.1.1 - Excursion Test Statistic
def excursion(s, mean):
    L = len(s)
    # 1. calculate the average of the sample values
    #    (passed as argument)
    # 2. For i=1 to L, calculate d_i
    # 3. T = max(d_1, ..., d_L)

    di = 0
    maximum = 0
    runningsum = 0
    for i in range(L):
        runningsum += s[i]
        di = runningsum-(i+1)*mean

        if di > maximum:
            maximum = di
        elif -di > maximum:
            maximum = -di
       
    return maximum

# Section 5.1.2 - Number of Directional Runs        
def numDirectionalRuns(s):
    L = len(s)

    # 1. Construct the sequence s'
    # We've given s' as input. just double check
    if set(s) != set([1, -1]):
        print ("error - input string is not sequence of 1's and -1's")
        quit()

    # 2. The test statistic T is the number of runs in s'
    numruns = 1
    for i in range(1,L):
        if s[i] != s[i-1]:
            numruns += 1
            
    return numruns

#Section 5.1.3 - Length of Directional Runs
def lenDirectionalRuns(s):
    L = len(s)

    # 1. Construct the sequence s'
    # We've given s' as input. just double check
    if not set(s).issubset(set([1, -1])):#set(s) > set([1, -1]):
        print set(s)
        print ("error - input string is not sequence of 1's and -1's")
        quit()
    
    # 2. The test statistic T is the length of the longest run in s'
    maxrun = 0
    run = 1
    for i in range(1,L):
        if s[i] == s[i-1]:
            run += 1
        else:
            if run > maxrun:
                maxrun  = run
            run = 1

    # handle the last observed run
    if run > maxrun:
        maxrun  = run
    
    return maxrun

#Section 5.1.4 - Number of Increases and Decreases
def numIncreasesDecreases(s):
    L = len(s)

    # 1. Construct the sequence s'
    # We've given s' as input. just double check
    if set(s) != set([1, -1]):
        print ("error - input string is not sequence of 1's and -1's")
        quit()

        
    #2. calculate the number of +1 and -1 in s'; the test statisitic T is
    #   the maximum of these two numbers
    pos = s.count(1) 
    return max(pos, L-pos)

#Section 5.1.5 - Number of Runs Based on the Median
def numRunsMedian(sp):

    # Steps 1 and 2 found in function altSequence2

    # 3. The test statistic is the number of runs in s'
    numruns = 1
    for i in range(1,len(sp)):
        if sp[i] != sp[i-1]:
            numruns += 1
            
    return numruns

# Section 5.1.6 - Length of Runs Based on Median
def lenRunsMedian(sp):
    
    # Steps 1 and 2 found in function altSequence2

    # 3. The test statistic is the length of the longest run in s'
    maxrun = 0
    run = 1
    for i in range(1,len(sp)):
        if sp[i] == sp[i-1]:
            run += 1
        else:
            if run > maxrun:
                maxrun  = run
            run = 1
            
    # handle the last observed run
    if run > maxrun:
        maxrun  = run
    
    return maxrun


# Section 5.1.7 - Average Collision Test Statistic
# takes list of collisions, C, as input 
def avgCollision(C):
    #Steps 1 through 3 found in function findCollisions

     #4. The test statistic T is the average of all values in the list C.       
    return sum(C)/float(len(C))


# Section 5.1.8 - Maximum Collision Test Statistic
# takes list of collisions, C, as input 
def maxCollision(C):
    #Steps 1 through 3 found in function findCollisions
        
    #4. The test statistic T is the maximum value in the list C.
    return max(C)


#Section 5.1.9 - Periodicity Test Statistic
def periodicity(s, p):
    L = len(s)

    #1. Initialize T to 0
    T = 0

    #2. For i=1 to L-p
    #       If (s_i = s_(i+p)), increment T by one
    i = 0
    while i < (L-p):
        if s[i] == s[i+p]:
            T += 1
        i += 1    
        
    return T


#Section 5.1.10 - Covariance Test Statistic
def covariance(s,p):
    L = len(s)

    #1. Initialize T to 0
    T = 0

    #2. For i=1 to L-p
    #       T = T+s_i*s_(i+p)
    i = 0
    while i < (L-p):
        T = T+(s[i]*s[i+p])
        i += 1
        
    return T


#Section 5.1.11 - Compression Test Statistics
def compression(s):
    # 1. Samples in the data subset are encoded as a character string
    #    containing a list of values separated by a single space, e.g.,
    #    (144,21,139,0,0,15) becomes "144 21 139 0 0 15"
    encodeS = " ".join([str(i) for i in s])

    # 2. The character string is processed with the BZ2 compression algorithm.
    compressed_string = bz2.compress(encodeS.encode())

    # 3. The score returned is the length of the compressed string, in bytes.
    return (len(compressed_string),)


def get_test_names():
    
    test_names = ('excursion', 'numDirectionalRuns', 'lenDirectionalRuns', 'numIncreasesDecreases', \
                  'numRunsMedian','lenRunsMedian','avgCollision','maxCollision','periodicity(1)',\
                  'periodicity(2)','periodicity(8)','periodicity(16)','periodicity(32)',\
                  'covariance(1)','covariance(2)','covariance(8)','covariance(16)','covariance(32)',\
                  'compression')

    return test_names


# create alternate sequence for several statistics
def altSequence1(s):
    L = len(s)
    sp = []
    for i in range(L-1):
        if s[i] > s[i+1]:
            sp.append(-1)
        else:
            sp.append(1)
    return sp


# create alternate sequence for median-based statistics
def altSequence2(s, X):
    L = len(s)

    #1. Find the median X of S (argument into function)
    #2. construct a temporary sequence s'
    sp = []
    for i in range(L-1):
        if s[i] < X:
            sp.append(-1)
        else:
            sp.append(1)
    return sp


# Finds list of collisions for collision statistics
def findCollisions(s):
    L = len(s)
    # 1. Let C be a list of the number of the samples observed to find two occurrences
    #    of the same value in the input sequence. C is initially empty.
    C = []

    #2. Let i = 1
    i = 0 #adjusted for index starting at 0

    #3. While i<L
    #   3a. Find the smallest j such that (s_i,...s_(i+j-1)) contains two identical values.
    #       If no such j exists, break out of the while loop.
    #   3b. Add j to the list C
    #   3c. i = i+j +1
    i = 0
    for ell, index in enumerate(s):
        if index in s[i:ell]:
            C.append(ell - i + 1) # j = ell-i+1
            i = ell + 1 
    
    return C

# Permutation testing 
def permutation_test(s, verbose=False):

    mean, median, is_binary = calc_stats(s)

    # 1. For each test i
    # 1.1 Assign the counters to zero
    C = {}
    for i in get_test_names():
        C[i] = [0,0]

    if verbose:
        print ("Calculating statistics on original sequence")

    # 1.2 Calculate the test statistic T_i on s; denote the result as t[i]
    t = {}

    # conversions for binary data
    if is_binary:
        cs1 = conversion1(s)
        cs2 = conversion2(s)

    t['excursion'] = excursion(s, mean)
    
    if is_binary:
        # conversion 1 required for binary data for directional run and
        # increases/decreases statistics
        sp1 = altSequence1(cs1)
        t['numDirectionalRuns'] = numDirectionalRuns(sp1)
        t['lenDirectionalRuns'] = lenDirectionalRuns(sp1)
        t['numIncreasesDecreases'] = numIncreasesDecreases(sp1)
    else:
        sp1 = altSequence1(s)
        t['numDirectionalRuns'] = numDirectionalRuns(sp1)
        t['lenDirectionalRuns'] = lenDirectionalRuns(sp1)
        t['numIncreasesDecreases'] = numIncreasesDecreases(sp1)

    if is_binary:
        sp2= altSequence2(cs2, 0.5)
    else:
        sp2= altSequence2(s, median)
    t['numRunsMedian'] = numRunsMedian(sp2)
    t['lenRunsMedian'] = lenRunsMedian(sp2)

    if is_binary:
        # conversion 2 required for collision statistics on binary data
        collision_list = findCollisions(cs2)
    else:
        collision_list = findCollisions(s)
    t['avgCollision'] = avgCollision(collision_list)
    t['maxCollision'] = maxCollision(collision_list)

    if is_binary:
        # conversion 1 required for periodicity statistics on binary data
        t['periodicity(1)'] = periodicity(cs1,1)
        t['periodicity(2)'] = periodicity(cs1,2)
        t['periodicity(8)'] = periodicity(cs1,8)
        t['periodicity(16)'] = periodicity(cs1,16)
        t['periodicity(32)'] = periodicity(cs1,32)
    else:
        t['periodicity(1)'] = periodicity(s,1)
        t['periodicity(2)'] = periodicity(s,2)
        t['periodicity(8)'] = periodicity(s,8)
        t['periodicity(16)'] = periodicity(s,16)
        t['periodicity(32)'] = periodicity(s,32)
        
    if is_binary:
        # conversion 1 required for covariance statistics on binary data
        t['covariance(1)'] = covariance(cs1,1)
        t['covariance(2)'] = covariance(cs1,2)
        t['covariance(8)'] = covariance(cs1,8)
        t['covariance(16)'] = covariance(cs1,16)
        t['covariance(32)'] = covariance(cs1,32)
    else:
        t['covariance(1)'] = covariance(s,1)
        t['covariance(2)'] = covariance(s,2)
        t['covariance(8)'] = covariance(s,8)
        t['covariance(16)'] = covariance(s,16)
        t['covariance(32)'] = covariance(s,32)
        
    t['compression'] = compression(s)

    # 2. For j=1 to 10000
        # 2.1 Permute the input data using Fisher-Yates
    
    if verbose:
        print ("Calculating statistics on permuted sequences")
    for j in range(10000):
        if verbose:
            sys.stdout.write("\rpermutation tests:\t%.2f percent complete" % (float(j)/100))
            sys.stdout.flush()
        shuffle(s)

        # 2.2 for each test i
        # 2.2.1 Calcualte the test statistic; denote the result as tp[i]
        tp = {}

        # conversions for binary data
        if is_binary:
            cs1 = conversion1(s)
            cs2 = conversion2(s)

        tp['excursion'] = excursion(s, mean)
        
        if is_binary:
            # conversion 1 required for binary data for directional run and
            # increases/decreases statistics
            sp1 = altSequence1(cs1)
            tp['numDirectionalRuns'] = numDirectionalRuns(sp1)
            tp['lenDirectionalRuns'] = lenDirectionalRuns(sp1)
            tp['numIncreasesDecreases'] = numIncreasesDecreases(sp1)
        else:
            sp1 = altSequence1(s)
            tp['numDirectionalRuns'] = numDirectionalRuns(sp1)
            tp['lenDirectionalRuns'] = lenDirectionalRuns(sp1)
            tp['numIncreasesDecreases'] = numIncreasesDecreases(sp1)

        if is_binary:
            sp2= altSequence2(cs2, 0.5)
        else:
            sp2= altSequence2(s, median)
        tp['numRunsMedian'] = numRunsMedian(sp2)
        tp['lenRunsMedian'] = lenRunsMedian(sp2)

        if is_binary:
            # conversion 2 required for collision statistics on binary data
            collision_list = findCollisions(cs2)
        else:
            collision_list = findCollisions(s)
        tp['avgCollision'] = avgCollision(collision_list)
        tp['maxCollision'] = maxCollision(collision_list)

        if is_binary:
            # conversion 1 required for periodicity statistics on binary data
            tp['periodicity(1)'] = periodicity(cs1,1)
            tp['periodicity(2)'] = periodicity(cs1,2)
            tp['periodicity(8)'] = periodicity(cs1,8)
            tp['periodicity(16)'] = periodicity(cs1,16)
            tp['periodicity(32)'] = periodicity(cs1,32)
        else:
            tp['periodicity(1)'] = periodicity(s,1)
            tp['periodicity(2)'] = periodicity(s,2)
            tp['periodicity(8)'] = periodicity(s,8)
            tp['periodicity(16)'] = periodicity(s,16)
            tp['periodicity(32)'] = periodicity(s,32)
            
        if is_binary:
            # conversion 1 required for covariance statistics on binary data
            tp['covariance(1)'] = covariance(cs1,1)
            tp['covariance(2)'] = covariance(cs1,2)
            tp['covariance(8)'] = covariance(cs1,8)
            tp['covariance(16)'] = covariance(cs1,16)
            tp['covariance(32)'] = covariance(cs1,32)
        else:
            tp['covariance(1)'] = covariance(s,1)
            tp['covariance(2)'] = covariance(s,2)
            tp['covariance(8)'] = covariance(s,8)
            tp['covariance(16)'] = covariance(s,16)
            tp['covariance(32)'] = covariance(s,32)
            
        tp['compression'] = compression(s)
        
        # 2.2.2 If (tp[i] > t[i]), increment C[i,0]. If (tp[i] = t[i]), increment C[i,1].
        for i in get_test_names():
            if tp[i] > t[i]:
                C[i][0] += 1
            elif tp[i] == t[i]:
                C[i][1] += 1

    if verbose:
        print ("\nstatistic\t\t\tC[i][0]  C[i][1]")
        for i in get_test_names():
            print ("%25s %8d %8d" % (i, C[i][0], C[i][1]))

    # 3. If (C[i][0]+C[i][1] <= 5) or (C[i][0] >= 9995) for any i, reject the IID assumption
    #    Else assume the noise source produces IID output.
    for i in get_test_names():
        if ((C[i][0]+C[i][1]) <= 5) or (C[i][0] >= 9995):
            return (False, C[i][0], C[i][1])
        else:
            return (True, C[i][0], C[i][1])
