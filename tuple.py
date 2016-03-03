# Non-iid t-tuple estimate defined in Draft SP 800-90B (January 2016)
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# Kerry McKay
# CSD/ITL/NIST
# February 2016

import math
from collections import Counter

def find_tuples(s, t):
    return zip(*[s[i:] for i in range(t)])

def t_tuple(s, verbose=False):
    L = len(s)
  
    # 1. Find the largest t such that the number of occurrences of the most
    # common t-tuple in S is at least 35.
    #
    # 2. Let Q[i] store the number of occurrences of the most common i-tuple in
    # S for i=1, ...,t. For example, in S=(2, 2, 0, 1, 0, 2, 0, 1, 2, 1, 2, 0, 1,
    # 2, 1, 0, 0, 1, 0, 0, 0), Q[1] = max(#0's,#1's,#2's)=9 and Q[2]=4 is
    # obtained by the number of 01's. 
    Q = []
    i = 0
    while True:
        i += 1
        c = Counter(find_tuples(s, i))
        if verbose:
          print "\ntuples (and counts):", c
          print "max occurances:",c.most_common(1)
        count = c.most_common(1)[0][1]
        Q.append(c.most_common(1)[0][1])

        if Q[i-1]<35:
            break

    Q.pop()
    t = i-1

    if verbose:
        print "\nt:", t
        print "Q:", Q


    # 3. For i=1 to t, an estimate for pmax is computed as:
    #   Let P[i]=Q[i] / (L-i+1), and compute an estimate on the maximum
    #   individual sample value probability as P_max[i]=P[i]^(1/i).
    P = [Q[i]/float(L-i) for i in range(t)]
    Pmax = [P[i]**(1.0/(i+1)) for i in range(t)]

    if verbose:
        print "P:", P
        print "Pmax:", Pmax
    
    # 4. The entropy estimate is calculated as -log_2 max(P_max[1],...,P_max[t]).
    min_h = -math.log(max(Pmax),2)
    if verbose:
      print "min-entropy:",min_h

    return max(Pmax), min_h


