# Longest Repeated Substring functions
# Sections 5.2.5 and 6.3.6 of DRAFT NIST SP 800-90B (January 2016)
#
# Kerry McKay
# CSD/ITL/NIST

from collections import Counter
import math
import itertools
from tuple import find_tuples


# Length of the Longest Repeated Substring Test - Section 5.2.5
# This test checks the IID assumption using the length of the longest repeated
# substring. If this length is significantly longer than the expected value,
# then the test invalidates the IID assumption. The test can be applied to
# binary and non-binary inputs.
def lenLRS(s, verbose=True):
    if verbose:
        print "LRS test"
    L = len(s)

    # 1. Find the proportion p_i of each possible input value x_i in S
    c = Counter(s)
    p = [x/float(L) for x in c.values()]

    # 2. calculate the collision probability
    p_col = sum(map(lambda x: x ** 2, p))

    # 3. find the length of the longest repeated substring W
    W = find_v(s,1,2)

    # 4. Calculate the number of overlapping subsequences of length W in S and
    #    the number of pairs as (L-W+1) choose 2
    N = (L-W+1)*(L-W) / 2 

    # 5. Let E be a binomially distributed random variable with parameters
    #    N = L-W+1 choose 2 and a probability of success p_col. Calculate the
    #    probability that E is greater than or equal to 1.
    PrE = 1-math.pow(1-p_col,N)
    if verbose:
        print "\tW:",W,"Pr(E>=1):",PrE

    if PrE < 0.001:
        return False
    else:
        return True



#helper function to find u in LRS estimate
def find_u(s, threshold):
    count = threshold+1 #loop control
    u = 0
    while count >= threshold:
        u += 1
        c = Counter(find_tuples(s, u))
        count = c.most_common(1)[0][1]
    return u


# helper function to find v in LRS estimate
# starts at u
def find_v(s, u, threshold):
    L = len(s)
    count = threshold+1 #loop control
    v = u
    while count >= threshold:
        v += 1
        c = Counter(find_tuples(s, v))
        count = c.most_common(1)[0][1]
    return v-1


# Longest Repeated Substring (LRS) Estimate - Section 6.3.6
def LRS_estimate(s, verbose=False):
    L = len(s)

    # 1. Find the smallest u such that the number of occurrences of the most
    #    common u-tuple in S is less than 20.
    u = find_u(s, 20)
    if verbose:
        print "u:",u

    # 2. Find the largest v such that the number of occurrences of the most
    #    common v-tuple in S is at least 2 and the most common (v+1)-tuple in S
    #    occurs once. 
    v = find_v(s, u, 2)
    if verbose:
        print "v:",v

    # 3. For W=u to v, compute the estimated W-tuple collision probability
    # 4. For each P_w, compute the estimated average collision probability per
    #    string symbol
    P = []
    for W in range(u,v+1):
        C = Counter(find_tuples(s, W))
        numerator = 0
        for c in C.values():
            #sum occurences of unique W-tuples. Only occurance counts >=2 contribute
            if c == 2:
                #save some computation time. 2 choose 2 is 1
                numerator += 1
            elif c == 3:
                #save some computation time. 3 choose 2 is 3
                numerator += 3
            elif c > 3:
                numerator += math.factorial(c)/2/math.factorial(c-2)
        denom = (L-W+1)*(L-W) / 2 
        P.append(math.pow(float(numerator)/denom, 1.0/W))
        

    # The entropy estimate is calculated as -log_2 max(Pmax)
    return max(P), -math.log(max(P),2)

