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
    # handle complete blocks
    sp = [s[i:i+8].count(1) for i in range(0,len(s)-len(s)%8,8)]

    # handle incomplete block
    if len(s)%8 > 0:
        sp.append(s[len(sp)*8:].count(1))
    
    return sp

# conversion II
# partitions the sequences into 8-bit non-overlapping blocks, and calculates the integer value of each block. 
def conversion2(s):
    sp = str(s)[1:-1]
    sp = sp.replace(', ','')
    
    if len(s) % 8 == 0:
        padding = 0
    else:
        padding = 8 - len(s) % 8

    padded = sp + '0' * padding
    sp_bytes = [int(padded[i:i+8], 2) for i in range(0,len(padded),8)]

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
    for i in range(L):
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
    if verbose:
      print ("t['excursion'] %s" % (str (t['excursion'])))
    
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
    if verbose:
      print ("t['numDirectionalRuns'] %s" % (str (t['numDirectionalRuns'])))
      print ("t['lenDirectionalRuns'] %s" % (str (t['lenDirectionalRuns'])))
      print ("t['numIncreasesDecreases'] %s" % (str (t['numIncreasesDecreases'])))

    if is_binary:
        sp2= altSequence2(s, 0.5)
    else:
        sp2= altSequence2(s, median)
    t['numRunsMedian'] = numRunsMedian(sp2)
    t['lenRunsMedian'] = lenRunsMedian(sp2)
    if verbose:
      print ("t['numRunsMedian'] %s" % (str (t['numRunsMedian'])))
      print ("t['lenRunsMedian'] %s" % (str (t['lenRunsMedian'])))

    if is_binary:
        # conversion 2 required for collision statistics on binary data
        collision_list = findCollisions(cs2)
    else:
        collision_list = findCollisions(s)
    t['avgCollision'] = avgCollision(collision_list)
    t['maxCollision'] = maxCollision(collision_list)
    if verbose:
      print ("t['avgCollision'] %s" % (str (t['avgCollision'])))
      print ("t['maxCollision'] %s" % (str (t['maxCollision'])))

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
    if verbose:
      print ("t['periodicity(1)'] %s" % (str (t['periodicity(1)'])))
      print ("t['periodicity(2)'] %s" % (str (t['periodicity(2)'])))
      print ("t['periodicity(8)'] %s" % (str (t['periodicity(8)'])))
      print ("t['periodicity(16)'] %s" % (str (t['periodicity(16)'])))
      print ("t['periodicity(32)'] %s" % (str (t['periodicity(32)'])))
        
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
    if verbose:
      print ("t['covariance(1)'] %s" % (str (t['covariance(1)'])))
      print ("t['covariance(2)'] %s" % (str (t['covariance(2)'])))
      print ("t['covariance(8)'] %s" % (str (t['covariance(8)'])))
      print ("t['covariance(16)'] %s" % (str (t['covariance(16)'])))
      print ("t['covariance(32)'] %s" % (str (t['covariance(32)'])))
        
    t['compression'] = compression(s)
    if verbose:
      print ("t['compression'] %s\n" % (str (t['compression'])))

    # 2. For j=1 to 10000
        # 2.1 Permute the input data using Fisher-Yates

    if verbose:
        print ("Calculating statistics on permuted sequences")
    for j in range(10000):
        shuffle(s)

        # 2.2 for each test i
        # 2.2.1 Calcualte the test statistic; denote the result as tp[i]
        tp = {}

        # conversions for binary data
        if is_binary:
            cs1 = conversion1(s)
            cs2 = conversion2(s)

        if ((C['excursion'][0]+C['excursion'][1]) <= 5) or (((j + 1) - C['excursion'][0]) <= 5):
          tp['excursion'] = excursion(s, mean)
        
        if is_binary:
            # conversion 1 required for binary data for directional run and
            # increases/decreases statistics
            sp1 = altSequence1(cs1)
            if ((C['numDirectionalRuns'][0]+C['numDirectionalRuns'][1]) <= 5) or (((j + 1) - C['numDirectionalRuns'][0]) <= 5):
              tp['numDirectionalRuns'] = numDirectionalRuns(sp1)
            if ((C['lenDirectionalRuns'][0]+C['lenDirectionalRuns'][1]) <= 5) or (((j + 1) - C['lenDirectionalRuns'][0]) <= 5):
              tp['lenDirectionalRuns'] = lenDirectionalRuns(sp1)
            if ((C['numIncreasesDecreases'][0]+C['numIncreasesDecreases'][1]) <= 5) or (((j + 1) - C['numIncreasesDecreases'][0]) <= 5):
              tp['numIncreasesDecreases'] = numIncreasesDecreases(sp1)
        else:
            sp1 = altSequence1(s)
            if ((C['numDirectionalRuns'][0]+C['numDirectionalRuns'][1]) <= 5) or (((j + 1) - C['numDirectionalRuns'][0]) <= 5):
              tp['numDirectionalRuns'] = numDirectionalRuns(sp1)
            if ((C['lenDirectionalRuns'][0]+C['lenDirectionalRuns'][1]) <= 5) or (((j + 1) - C['lenDirectionalRuns'][0]) <= 5):
              tp['lenDirectionalRuns'] = lenDirectionalRuns(sp1)
            if ((C['numIncreasesDecreases'][0]+C['numIncreasesDecreases'][1]) <= 5) or (((j + 1) - C['numIncreasesDecreases'][0]) <= 5):
              tp['numIncreasesDecreases'] = numIncreasesDecreases(sp1)

        if is_binary:
            sp2= altSequence2(s, 0.5)
        else:
            sp2= altSequence2(s, median)
        if ((C['numRunsMedian'][0]+C['numRunsMedian'][1]) <= 5) or (((j + 1) - C['numRunsMedian'][0]) <= 5):
          tp['numRunsMedian'] = numRunsMedian(sp2)
        if ((C['lenRunsMedian'][0]+C['lenRunsMedian'][1]) <= 5) or (((j + 1) - C['lenRunsMedian'][0]) <= 5):
          tp['lenRunsMedian'] = lenRunsMedian(sp2)

        if is_binary:
            # conversion 2 required for collision statistics on binary data
            if (((C['avgCollision'][0]+C['avgCollision'][1]) <= 5) or (((j + 1) - C['avgCollision'][0]) <= 5)
              or ((C['maxCollision'][0]+C['maxCollision'][1]) <= 5) or (((j + 1) - C['maxCollision'][0]) <= 5)):
              collision_list = findCollisions(cs2)
        else:
            if (((C['avgCollision'][0]+C['avgCollision'][1]) <= 5) or (((j + 1) - C['avgCollision'][0]) <= 5)
              or ((C['maxCollision'][0]+C['maxCollision'][1]) <= 5) or (((j + 1) - C['maxCollision'][0]) <= 5)):
              collision_list = findCollisions(s)
        if ((C['avgCollision'][0]+C['avgCollision'][1]) <= 5) or (((j + 1) - C['avgCollision'][0]) <= 5):
          tp['avgCollision'] = avgCollision(collision_list)
        if ((C['maxCollision'][0]+C['maxCollision'][1]) <= 5) or (((j + 1) - C['maxCollision'][0]) <= 5):
          tp['maxCollision'] = maxCollision(collision_list)

        if is_binary:
            # conversion 1 required for periodicity statistics on binary data
            if ((C['periodicity(1)'][0]+C['periodicity(1)'][1]) <= 5) or (((j + 1) - C['periodicity(1)'][0]) <= 5):
              tp['periodicity(1)'] = periodicity(cs1,1)
            if ((C['periodicity(2)'][0]+C['periodicity(2)'][1]) <= 5) or (((j + 1) - C['periodicity(2)'][0]) <= 5):
              tp['periodicity(2)'] = periodicity(cs1,2)
            if ((C['periodicity(8)'][0]+C['periodicity(8)'][1]) <= 5) or (((j + 1) - C['periodicity(8)'][0]) <= 5):
              tp['periodicity(8)'] = periodicity(cs1,8)
            if ((C['periodicity(16)'][0]+C['periodicity(16)'][1]) <= 5) or (((j + 1) - C['periodicity(16)'][0]) <= 5):
              tp['periodicity(16)'] = periodicity(cs1,16)
            if ((C['periodicity(32)'][0]+C['periodicity(32)'][1]) <= 5) or (((j + 1) - C['periodicity(32)'][0]) <= 5):
              tp['periodicity(32)'] = periodicity(cs1,32)
        else:
            if ((C['periodicity(1)'][0]+C['periodicity(1)'][1]) <= 5) or (((j + 1) - C['periodicity(1)'][0]) <= 5):
              tp['periodicity(1)'] = periodicity(s,1)
            if ((C['periodicity(2)'][0]+C['periodicity(2)'][1]) <= 5) or (((j + 1) - C['periodicity(2)'][0]) <= 5):
              tp['periodicity(2)'] = periodicity(s,2)
            if ((C['periodicity(8)'][0]+C['periodicity(8)'][1]) <= 5) or (((j + 1) - C['periodicity(8)'][0]) <= 5):
              tp['periodicity(8)'] = periodicity(s,8)
            if ((C['periodicity(16)'][0]+C['periodicity(16)'][1]) <= 5) or (((j + 1) - C['periodicity(16)'][0]) <= 5):
              tp['periodicity(16)'] = periodicity(s,16)
            if ((C['periodicity(32)'][0]+C['periodicity(32)'][1]) <= 5) or (((j + 1) - C['periodicity(32)'][0]) <= 5):
              tp['periodicity(32)'] = periodicity(s,32)
            
        if is_binary:
            # conversion 1 required for covariance statistics on binary data
            if ((C['covariance(1)'][0]+C['covariance(1)'][1]) <= 5) or (((j + 1) - C['covariance(1)'][0]) <= 5):
              tp['covariance(1)'] = covariance(cs1,1)
            if ((C['covariance(2)'][0]+C['covariance(2)'][1]) <= 5) or (((j + 1) - C['covariance(2)'][0]) <= 5):
              tp['covariance(2)'] = covariance(cs1,2)
            if ((C['covariance(8)'][0]+C['covariance(8)'][1]) <= 5) or (((j + 1) - C['covariance(8)'][0]) <= 5):
              tp['covariance(8)'] = covariance(cs1,8)
            if ((C['covariance(16)'][0]+C['covariance(16)'][1]) <= 5) or (((j + 1) - C['covariance(16)'][0]) <= 5):
              tp['covariance(16)'] = covariance(cs1,16)
            if ((C['covariance(32)'][0]+C['covariance(32)'][1]) <= 5) or (((j + 1) - C['covariance(32)'][0]) <= 5):
              tp['covariance(32)'] = covariance(cs1,32)
        else:
            if ((C['covariance(1)'][0]+C['covariance(1)'][1]) <= 5) or (((j + 1) - C['covariance(1)'][0]) <= 5):
              tp['covariance(1)'] = covariance(s,1)
            if ((C['covariance(2)'][0]+C['covariance(2)'][1]) <= 5) or (((j + 1) - C['covariance(2)'][0]) <= 5):
              tp['covariance(2)'] = covariance(s,2)
            if ((C['covariance(8)'][0]+C['covariance(8)'][1]) <= 5) or (((j + 1) - C['covariance(8)'][0]) <= 5):
              tp['covariance(8)'] = covariance(s,8)
            if ((C['covariance(16)'][0]+C['covariance(16)'][1]) <= 5) or (((j + 1) - C['covariance(16)'][0]) <= 5):
              tp['covariance(16)'] = covariance(s,16)
            if ((C['covariance(32)'][0]+C['covariance(32)'][1]) <= 5) or (((j + 1) - C['covariance(32)'][0]) <= 5):
              tp['covariance(32)'] = covariance(s,32)
            
        if ((C['compression'][0]+C['compression'][1]) <= 5) or (((j + 1) - C['compression'][0]) <= 5):
          tp['compression'] = compression(s)
        
        # 2.2.2 If (tp[i] > t[i]), increment C[i,0]. If (tp[i] = t[i]), increment C[i,1].
        pass_flag = True
        for i in get_test_names():
            if ((C[i][0]+C[i][1]) <= 5) or (((j + 1) - C[i][0]) <= 5):
              if tp[i] > t[i]:
                  C[i][0] += 1
              elif tp[i] == t[i]:
                  C[i][1] += 1
              pass_flag = False
        if verbose:
          if ((C['excursion'][0]+C['excursion'][1]) <= 5) or (((j + 1) - C['excursion'][0]) <= 5):
            print ("t tp ['excursion'] %d  %d  %d" % (t['excursion'], tp['excursion'], t['excursion'] - tp['excursion']))
          if ((C['numDirectionalRuns'][0]+C['numDirectionalRuns'][1]) <= 5) or (((j + 1) - C['numDirectionalRuns'][0]) <= 5):
            print ("t tp ['numDirectionalRuns'] %d  %d  %d" % (t['numDirectionalRuns'], tp['numDirectionalRuns'],
              t['numDirectionalRuns'] - tp['numDirectionalRuns']))
          if ((C['lenDirectionalRuns'][0]+C['lenDirectionalRuns'][1]) <= 5) or (((j + 1) - C['lenDirectionalRuns'][0]) <= 5):
            print ("t tp ['lenDirectionalRuns'] %d  %d  %d" % (t['lenDirectionalRuns'], tp['lenDirectionalRuns'],
              t['lenDirectionalRuns'] - tp['lenDirectionalRuns']))
          if ((C['numIncreasesDecreases'][0]+C['numIncreasesDecreases'][1]) <= 5) or (((j + 1) - C['numIncreasesDecreases'][0]) <= 5):
            print ("t tp ['numIncreasesDecreases'] %d  %d  %d" % (t['numIncreasesDecreases'], tp['numIncreasesDecreases'],
              t['numIncreasesDecreases'] - tp['numIncreasesDecreases']))
          if ((C['numRunsMedian'][0]+C['numRunsMedian'][1]) <= 5) or (((j + 1) - C['numRunsMedian'][0]) <= 5):
            print ("t tp ['numRunsMedian'] %d  %d  %d" % (t['numRunsMedian'], tp['numRunsMedian'], t['numRunsMedian'] - tp['numRunsMedian']))
          if ((C['lenRunsMedian'][0]+C['lenRunsMedian'][1]) <= 5) or (((j + 1) - C['lenRunsMedian'][0]) <= 5):
            print ("t tp ['lenRunsMedian'] %d  %d  %d" % (t['lenRunsMedian'], tp['lenRunsMedian'], t['lenRunsMedian'] - tp['lenRunsMedian']))
          if ((C['avgCollision'][0]+C['avgCollision'][1]) <= 5) or (((j + 1) - C['avgCollision'][0]) <= 5):
            print ("t tp ['avgCollision'] %8.5e  %8.5e  %8.5e" % (t['avgCollision'], tp['avgCollision'], t['avgCollision'] - tp['avgCollision']))
          if ((C['maxCollision'][0]+C['maxCollision'][1]) <= 5) or (((j + 1) - C['maxCollision'][0]) <= 5):
            print ("t tp ['maxCollision'] %d  %d  %d" % (t['maxCollision'], tp['maxCollision'], t['maxCollision'] - tp['maxCollision']))
          if ((C['periodicity(1)'][0]+C['periodicity(1)'][1]) <= 5) or (((j + 1) - C['periodicity(1)'][0]) <= 5):
            print ("t tp ['periodicity(1)'] %d  %d  %d" % (t['periodicity(1)'], tp['periodicity(1)'], t['periodicity(1)'] - tp['periodicity(1)']))
          if ((C['periodicity(2)'][0]+C['periodicity(2)'][1]) <= 5) or (((j + 1) - C['periodicity(2)'][0]) <= 5):
            print ("t tp ['periodicity(2)'] %d  %d  %d" % (t['periodicity(2)'], tp['periodicity(2)'], t['periodicity(2)'] - tp['periodicity(2)']))
          if ((C['periodicity(8)'][0]+C['periodicity(8)'][1]) <= 5) or (((j + 1) - C['periodicity(8)'][0]) <= 5):
            print ("t tp ['periodicity(8)'] %d  %d  %d" % (t['periodicity(8)'], tp['periodicity(8)'], t['periodicity(8)'] - tp['periodicity(8)']))
          if ((C['periodicity(16)'][0]+C['periodicity(16)'][1]) <= 5) or (((j + 1) - C['periodicity(16)'][0]) <= 5):
            print ("t tp ['periodicity(16)'] %d  %d  %d" % (t['periodicity(16)'], tp['periodicity(16)'], 
              t['periodicity(16)'] - tp['periodicity(16)']))
          if ((C['periodicity(32)'][0]+C['periodicity(32)'][1]) <= 5) or (((j + 1) - C['periodicity(32)'][0]) <= 5):
            print ("t tp ['periodicity(32)'] %d  %d  %d" % (t['periodicity(32)'], tp['periodicity(32)'], 
              t['periodicity(32)'] - tp['periodicity(32)']))
          if ((C['covariance(1)'][0]+C['covariance(1)'][1]) <= 5) or (((j + 1) - C['covariance(1)'][0]) <= 5):
            print ("t tp ['covariance(1)'] %d  %d  %d" % (t['covariance(1)'], tp['covariance(1)'], t['covariance(1)'] - tp['covariance(1)']))
          if ((C['covariance(2)'][0]+C['covariance(2)'][1]) <= 5) or (((j + 1) - C['covariance(2)'][0]) <= 5):
            print ("t tp ['covariance(2)'] %d  %d  %d" % (t['covariance(2)'], tp['covariance(2)'], t['covariance(2)'] - tp['covariance(2)']))
          if ((C['covariance(8)'][0]+C['covariance(8)'][1]) <= 5) or (((j + 1) - C['covariance(8)'][0]) <= 5):
            print ("t tp ['covariance(8)'] %d  %d  %d" % (t['covariance(8)'], tp['covariance(8)'], t['covariance(8)'] - tp['covariance(8)']))
          if ((C['covariance(16)'][0]+C['covariance(16)'][1]) <= 5) or (((j + 1) - C['covariance(16)'][0]) <= 5):
            print ("t tp ['covariance(16)'] %d  %d  %d" % (t['covariance(16)'], tp['covariance(16)'], t['covariance(16)'] - tp['covariance(16)']))
          if ((C['covariance(32)'][0]+C['covariance(32)'][1]) <= 5) or (((j + 1) - C['covariance(32)'][0]) <= 5):
            print ("t tp ['covariance(32)'] %d  %d  %d" % (t['covariance(32)'], tp['covariance(32)'], t['covariance(32)'] - tp['covariance(32)']))
          if ((C['compression'][0]+C['compression'][1]) <= 5) or (((j + 1) - C['compression'][0]) <= 5):
            print ("t tp ['compression'] %d  %d  %d" % (t['compression'][0], tp['compression'][0], t['compression'][0] - tp['compression'][0]))
          for i in get_test_names():
            print ("%25s %8d %8d" % (i, C[i][0], C[i][1]))
          print ("\n")
          print("permutation tests:\t%.2f percent complete\n" % (float(j + 1)/100))
          print ("\n")
        if pass_flag:  # the permutation criteria have already been met
          return True

    if verbose:
        print("\n                statistic  C[i][0]  C[i][1]")
        print("-------------------------------------------")
        for i in get_test_names():
            if ((C[i][0]+C[i][1]) <= 5) or (C[i][0] >= 9995):
                print ("%24s* %8d %8d" % (i, C[i][0], C[i][1]))
            else:
                print ("%25s %8d %8d" % (i, C[i][0], C[i][1]))
        print("(* denotes failed test)")

    # 3. If (C[i][0]+C[i][1] <= 5) or (C[i][0] >= 9995) for any i, reject the IID assumption
    #    Else assume the noise source produces IID output.
    for i in get_test_names():
        if ((C[i][0]+C[i][1]) <= 5) or (C[i][0] >= 9995):
            return False
    return True
