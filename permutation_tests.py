## @package shuffle_tests
#
# Shuffling Tests on Independence and Stability
# DRAFT NIST SP 800-90B (August 2012) Section 9.1.2
#
# Tim Hall
# tim.hall@nist.gov
#
# 2 September 2014

import bz2
import sys
from collections import namedtuple
from util90b import Rank

# Use numpy shuffle function if available
# Much faster
try:
    import numpy as np
    shuffle = np.random.shuffle
except:
    import random
    shuffle = random.shuffle



def shuffle_test_scores(shuffle_set):
    dataset = shuffle_set[0]
    number_of_shuffles = shuffle_set[1]

    # find mean, median, is dataset binary
    stats = calc_stats(dataset)
    
    tests = get_shuffle_tests()

    # compute all scores on original, unshuffled data set
    original_scores = {}
    shuffle_scores = {}
    for test in tests:
        scores = test.fcn(dataset, stats)
        original_scores[test.name] = scores
        shuffle_scores[test.name] = []

    # compute all scores on shuffled data sets
    for n in range(number_of_shuffles):
        shuffled = list(dataset)
        shuffle(shuffled)

        for test in tests:
            scores = test.fcn(shuffled, stats)
            shuffle_scores[test.name].append(scores)

    # for all scores from all tests, find the rank of the original dataset's
    # score in the shuffled datasets scores
    original_scores_and_ranks = {}
    for test in tests:
        # test may return more than one score
        r = []
        for i in range(len(test.scores)):
            score_i = [s[i] for s in shuffle_scores[test.name]]
            r.append(Rank(original_scores[test.name][i], sorted(score_i)))

        original_scores_and_ranks[test.name] = (original_scores[test.name], r)

    # return the original dataset scores + ranks
    return original_scores_and_ranks

    

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


# return a list of six shuffle tests, including their:
#  1. name
#  2. function pointer
#  3. score description
#  4. format strings for printing scores
def get_shuffle_tests():
    # store a test's: name, function pointer, score descriptions and score format
    Test = namedtuple('Test','name fcn scores format')
    
    # list of six shuffle tests
    tests = []
    tests.append(Test('compression',compression_test, ['bytes'], "%d"))
    tests.append(Test("over/under", over_under_runs_score, ['longest run', 'number of runs'], "%d  %d"))
    tests.append(Test("excursion", excursion_score, ['score'], "%g"))
    tests.append(Test("directional runs", directional_runs_scores, ['number of runs', 'longest run','max up/down'], "%d  %4d  %d"))
    tests.append(Test("covariance", covariance_score, ['score'], "%g"))
    tests.append(Test("collision", collision_score, ['min counts', 'avg counts','max counts'], "%d  %g  %d"))

    return tests


# Compression score - Section 9.1.2.1 of DRAFT NIST SP 800-90B
def compression_test(data, stats='not used'):
    # 1. Samples in the data subset are encoded as a character string
    #    containing a list of values separated by commas, e.g.,
    #    "144,21,139,0,0,15"
    encodeS = ",".join([str(s) for s in data])

    # 2. The character string is processed with the BZ2 compression algorithm.
    compressed_string = bz2.compress(encodeS.encode())

    # 3. The score returned is the length of the compressed string, in bytes.
    return (len(compressed_string),)


# Over/Under Runs Scores - Section 9.1.2.2
def over_under_runs_score(data_subset, stats):
    # 1. The data subset is used to compute a median for the values in this data
    # subset.  If the data subset is binary, then the median is 0.5
    median = stats[1]

    # 2. A temporary data subset is constructed for each of the original and
    # shuffled data subsets. For each element in the data subset:
    #   a. If element is larger than median, append +1
    #   b. If element is smaller than median, append -1
    #   c. If element is the same as median, do not append anything
    temp_subset = [-1 if x < median else 1 for x in data_subset if x != median]

    # 3. Determine the longest run of -1 or +1 in temporary data subset, let
    #    the length of this longest run be the first score.
    # 4. The number of runs of -1 and +1 values in the temporary data subset
    #    is reported as the second score.
    longest_run = 0
    number_of_runs = 0
    previous = temp_subset[0]
    runlength = 1
    for current in temp_subset[1:]:
        if previous == current:
            runlength += 1
        else:
            if runlength > longest_run: longest_run = runlength
            number_of_runs += 1
            runlength = 1
        previous = current

    # handle final run in subset
    if runlength > longest_run: longest_run = runlength
    number_of_runs += 1

    return longest_run, number_of_runs


# Excursion Score - Section 9.1.2.3
# This test measures the greatest distance the dataset strays away from its
# average behavior.
#
# This test is useful for detecting cases where the dataset has some trends,
# positive correlations between nearby samples, or where some underlying hidden
# variable is changing the average value of samples from time to time.  Unlike
# the over/under runs score, this test can differentiate between small and large
# deviations from the average behavior of the source.  
def excursion_score(input_sequence, stats):
    mean = stats[0]
    # 1. For j = 1 to <subset-length>,
    #    dj = absolute value of (sum of first j samples - j * mu)

    # 2. The score returned is the largest of the dj values
    dj = 0.0
    maximum = 0.0

    for x in input_sequence:
        dj += x - mean

        if dj > maximum:
            maximum = dj
        elif -dj > maximum:
            maximum = -dj

    return (maximum,)


# Directional Runs Scores - Section 9.1.2.4
def directional_runs_scores(dataset, stats):
    # 1. The temporary data subset temp is produced from the data subset
    #    as follows:
    #    a. If the input is not binary:
    #       For i = 0 to (the length of original data subset) - 2:
    #            If s_i < s_i+1, then temp_i = 1
    #            Else if s_i > s_i+1, then temp_i = -1
    #            Else temp_i = 0
    #    b. If the input is binary:
    #       ...first require processing in which bits are combined into bytes.
    #       Then, a new data subset is created from the Hamming weights of the
    #       successive bytes and then the temporary dataset is generated from
    #       this.
    #       For i to (length of original data subset)/8 - 1:
    #            W_i = hamming_weight(s_i,...,s_i+7)
    #       For i = 0 to (length of sequence W) - 2:
    #            If s_i < s_i+1, then temp_i = 1
    #            Else if s_i > s_i+1, then temp_i = -1
    #            Else temp_i = 0
    is_binary = stats[2]
    s = dataset
    if not is_binary:
        temp = [1 if s[i] < s[i+1] else -1 if s[i] > s[i+1] else 0 for i in range(len(s)-1)]
    else:
        W = [sum(s[i:i+8]) for i in range(0, len(s), 8)]
        temp = [1 if W[i]<W[i+1] else -1 if W[i]>W[i+1] else 0 for i in range(len(W)-1)]

    # 2. Calculate the scores on the temporary data subset
    runtype = 0 # runtype of 0 to handle 'temp' starting with a run of '0's
    current_run = 0
    longest_run = 0
    upcount = 0
    downcount = 0
    runs = 0
    for current in temp:
        # current symbol is a '1'
        if current == 1:
            upcount += 1
            if runtype == 1:
                current_run += 1
            else:
                # start a new run of '1's
                if current_run > longest_run:
                    longest_run = current_run
                runtype = 1
                current_run = 1
                runs += 1
        # current symbol is a '-1'
        elif current == -1:
            downcount += 1
            if runtype == -1:
                current_run += 1
            else:
                # start a new run of '-1's
                if current_run > longest_run:
                    longest_run = current_run
                runtype = -1
                current_run = 1
                runs += 1
        # current symbol is a '0'
        # 'runtype = 0' means we are still in a run of '0's
        # at the start of the 'temp' data subset
        elif runtype != 0:
            # increment length of current run, whether '1's or '-1's
            current_run += 1
        
    # handle last run
    if current_run > longest_run:
        longest_run = current_run

    return runs, longest_run, max(upcount, downcount)


# Covariance Score - Section 9.1.2.5
def covariance_score(data_subset, stats):
    # 1. count = 0
    count = 0.0

    # 2. mu = mean s_0,s_1,...,s_(N/10-1)
    mu = stats[0]

    # 3. for i = 1 to floor(N/10):
    #        count = count + (s_i - mu)(s_i-1 - mu)
    s = data_subset
    for i in range(1, len(s)):
        count += (s[i] - mu) * (s[i-1] - mu)

    # 4. Score = count/(floor(N/10) - 1)
    score = count/(len(s)-1.0)

    return (score,)


# convert a sequence of 8 bits into a byte value
def to_byte(bits):
    return sum([bits[i] << (7-i) for i in range(8)])

# Collision score - Section 9.1.2.6
def collision_score(data_subset, stats):
    # In case of binary data, a data subset consisting of binary data is
    # converted into a sequence of 8-bit bytes before being subjected to
    # this test, so that each sequence of eight bits in original subset becomes
    # one byte in the modified data subset
    is_binary = stats[2]
    if is_binary:
        N = (len(data_subset) // 8) * 8
        s = [to_byte(data_subset[i:i+8]) for i in range(0, N, 8)]
    else:
        s = data_subset

    # 1. Counts is a list of counts of samples needed to see a duplicate; the
    # list is initially empty
    counts = []

    # 2. pos = 0
    pos = 0

    # 3. While pos < (length of data subset):
    #    a. Find the smallest j s.t. s_pos .. s_pos+j contains one duplicate
    #       value.
    #       i. If no such j exists, break out of while loop
    for pos_j,s_pos_j in enumerate(s):
        if s_pos_j in s[pos:pos_j]:
            counts.append(pos_j-pos) # difference is j 
            pos = pos_j + 1

    # 4. Return the following values as scores: min,avg,max of counts
    return min(counts), sum(counts)/float(len(counts)), max(counts)
