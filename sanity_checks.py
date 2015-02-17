#
# Sanity Checks Against Entropy Estimates
# DRAFT NIST SP 800-90B (August 2012) Section 9.4
#
# Tim Hall
# tim.hall@nist.gov
# 18 Mar 2013
#

import sys
import bz2
import math


# Compression Sanity Check - Section 9.4.1 of DRAFT NIST SP 800-90B
#
# Let S be a collected dataset defined in Section 7.1.  Divide S into ten
# equal-length, non-overlapping data subsets S1,...,S10.  Let Si be one of the
# ten data subsets, and let H be the entropy estimate for the data source, as
# calculated in Section 9.2 or Section 9.3 (depending upon the test track
# taken).  The test is performed as follows on the ten data subsets of samples:
#
def compression_sanity_check(S, H, verbose=False):
    # Find length of a data subset, equal to floor(N/10)
    # Used to determine Si and in Step 4
    L = len(S) // 10

    if verbose:
        print("\nCompression sanity check...")
    # For i = 1 to 10:
    for i in range(10):
        Si = S[i*L:(i+1)*L]

        # 1. EncodeS = the encoding of Si as a string
        EncodeS = ",".join([str(s) for s in Si])

        # 2. C = compress(EncodeS)
        C = bz2.compress(EncodeS.encode())

        # 3. Len = the length of C in bits
        Len =  8 * len(C)

        # 4. If Len < H * floor(N/10), the dataset fails the test
        if verbose:
            sys.stdout.write("\tdataset %d compressed length = %d, cutoff = %g..." % (i+1,Len,H * L))
        if Len < (H * L):
            if verbose:
                print("Fail\n")
            return 'fail', i+1
        else:
            if verbose:
                print("Pass\n")

    return 'pass', None


# Function used in Step 1 below:
#
# b = Round((2 log2(N/b) -1)/H)
def bfunc(b, N, H):
    term = 2.0*math.log(float(N)/b, 2.0) - 1.0
    term /= float(H)
    return int(round(term))

# Poisson CDF
# http://www.itl.nist.gov/div898/handbook/eda/section3/eda366j.htm
# note - variable name is 'lambd' because 'lambda' is a keyword in Python
def Poisson_cdf(x, lambd):
    # term for i=0 = lambd^0/0!
    term_i = 1.0
    CDF = term_i
    # calc term and cumulative sum for i = 1 to x
    for i in range(1, x+1):
        term_i *= lambd/float(i)
        CDF += term_i

    CDF *= math.exp(-lambd)

    return CDF
    
        
# Collision Sanity Check - Section 9.4.2 of DRAFT NIST SP 800-90B
#
def collision_sanity_check(dataset, H, verbose=False):

    # 1. Using the complete dataset collected as specified in Section 7.1,
    #    determine the size of the 'tuples' to be examined for collisions.
    #    This test will examine each successive b-tuple of samples from the
    #    dataset.  To determine b, numerically solve using the following
    #    equation for b:
    #
    #        b = Round((2 log2(N/b) -1)/H)
    N = len(dataset)
    bp = 1
    b = bfunc(bp, N, H)
    while b > bp:
        bp += 1
        b = bfunc(bp, N, H)

    if verbose:
        print("\nCollision sanity check...")
        print("\tDividing dataset into %d-tuples" % (b))

    # 2. Divide the dataset into successive b-tuples and count the number of
    #    occurrences of each value of a b-tuple.
    count = {}
    for i in range(N//b):
        btuple = tuple(dataset[i*b:(i+1)*b])
        count[btuple] = count.get(btuple,0) + 1

    # 3. Determine whether the noise source and entropy estimate should be
    #    rejected

    #    Rule 1 - If three or more b-tuples have the same value, the source is
    #    rejected
    if verbose:
        sys.stdout.write("\tCheck rule 1 - do three or more %d-tuples have the same value?..." % (b))
    if max(count.values()) > 2:
        most = sorted(count.values())
        most.reverse()
        most = most[:3]
        if verbose:
            print("Fail")
        for btuple in count.keys():
            if count[btuple] > 2 and count[btuple] in most:
                if verbose:
                    print("\t\t%r occurs %d times" % (btuple, count[btuple]))
                most.remove(count[btuple])
        return 'fail'
    else:
        if verbose:
            print("Pass")

    #    Rule 2 - A p-value on the number of collisions is determined,
    #    approximating the probability of x collisions as a Poisson distribution
    #    with lambda = 2^(-Hb)(floor(N/b)^2 / 2) and choosing as the cutoff value
    #    the smallest number of collisions such that the probability of that
    #    many or more collisions is less than 1/10,000
    if verbose:
        print("\tCheck rule 2 - probability of number of collisions below cutoff")
    collisions = [1 for btuple in count.keys() if count[btuple] > 1]
    x = sum(collisions)
    
    # Use Poisson cdf to determine P(X > x)
    lambd = pow(2.0, -H * b) * (pow(N // b, 2)/ 2.0)
    cdf = Poisson_cdf(x, lambd)

    if verbose:
        sys.stdout.write("\t\tnumber of collisions = %d, cutoff = %g..." % (x, lambd))
    if (1.0 - cdf) < 1.0/10000.0:
        if verbose:
            print("Fail\n")
        return 'fail'
    else:
        if verbose:
            print("Pass\n")

    return 'pass'
