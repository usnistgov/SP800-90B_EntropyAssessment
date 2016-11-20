# The Collision Estimate
# Section 6.3.2 of DRAFT NIST SP 800-90B (2016)
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# T. A. Hall
# CSD/ITL/NIST
# 25 Feb 2013
#
# Updated by Kerry McKay, Nov 2015

import math
from is_close import isclose
import numpy as np
import sys

# Continued fraction for F(1/z) derived from Eq. 8.9.2 at http://dlmf.nist.gov/8.9
# Used in Step 9 in Section 6.3.2
def F(n, zInv):
    z = 1.0/zInv
    denom = 1.0 + n/z

    for i in range(1, n):
        denom = z + -i/denom
        denom = 1.0 + (n-i)/denom

    denom = z + -n/denom
    return 1.0/denom

# Expected value of statistic based on one-parameter family of prob distributions
# Used in step 9 of Section 6.3.2 (right side of equation)
def calcEpS(p, k):
    q = (1.0 - p)/float(k - 1)
    i_k = 1.0 / k
    ip = 1.0/p
    iq = 1.0/q
    iq2 = iq * iq
    Ep = p * iq2 * (1.0 + i_k * (ip - iq)) * F(k, q) - (p * iq * i_k * (ip - iq))
    return Ep


# Binary search for p that solves equation in step 9 of Section 6.3.2
#note, this presumes a decreasing function
def solve_for_p(mu_bar, n, tolerance=1e-09):
    assert n > 0

    if mu_bar > calcEpS(1.0/float(n), n):
        print ("Sought value is greater than largest possible value. Giving up.")
        return False, 0.0

    #This is a hackish way of checking to see if the difference is within approximately 4 ULPs
    absEpsilon = 4.0 * max((np.nextafter(mu_bar, mu_bar+1.0) - mu_bar), (mu_bar - np.nextafter(mu_bar, mu_bar-1.0)))
    #If we don't have numpy, then this will work for most of the ranges we're concerned with
    #absEpsilon = sys.float_info.epsilon

    ldomain = 1.0/float(n)
    hdomain = 1.0

    lbound = ldomain
    lvalue = float("inf")
    hbound = hdomain
    hvalue = float("-inf")

    #Note that the bounds are in the interval [0, 1], so underflows
    #are an issue, but overflows are not
    center = (lbound + hbound) / 2
    assert (center > ldomain) and (center < hdomain)

    centerVal = calcEpS(center, n)

    for j in range(1, 1076):
        if isclose(mu_bar, centerVal, tolerance, absEpsilon):
            return True, center

        if lbound >= hbound:
            print ("Bounds have converged after %d rounds and target was not found" % j)
            return False, 0.0

        if mu_bar < centerVal:
            lbound = center
            lvalue = centerVal
        else:
            hbound = center
            hvalue = centerVal

        if (mu_bar > lvalue) or (mu_bar < hvalue):
            print ("Target is not within the search interval after %d rounds" % j)
            return False, 0.0

        center = (lbound + hbound) / 2.0

        if (center <= ldomain) or (center >= hdomain):
            print ("The next center is outside of the proscribed domain after %d rounds" % j)
            return False, 0.0

        centerVal = calcEpS(center, n)

    if isclose(mu_bar, centerVal, tolerance, absEpsilon):
        return True, p
    else:
        print ("Binary search failed to converge, or the sought value is outside the interval.")
        return False, 0.0

# Section 6.3.2- Collision Estimate
def collision_test(s, n):
    # 1. Set v=1, index=1
    # 2. Beginning with s_index, step through the dataset until any observed
    #    value is repeated; i.e., find smallest j s.t. si = sj for some
    #    i with 1 <= i < j
    # 3. Set t_v = j-index+1, v=v+1, and indes=j+1
    # 4. Repeat steps 2-3 until the end of the dataset is reached
    # 5. set v = v - 1
    # 6. If v < 1000, the noise source outputs will be mapped down based on
    #    the ranking provided, and the data will be retested.

    k=len(set(s))
    index = 0
    t = [0]
    for ell, i in enumerate(s):
        if i in s[index:ell]:
            t.append(ell + 1) # account for index starting at 1
            index = ell + 1
    diff_t = [t[i]-t[i-1] for i in range(1,len(t))]
    v = float(len(diff_t)) # float() so it works in Python 2.6,2.7
    if v < 1000:
        print ("Must map down data for collision estimate, then retest")
        return 0.0, 0.0

    # 7. Calculate sample mean, mu, and sample stddev, sigma of the differences
    # of collision times
    mu = sum(diff_t)/v
    sigma = sum([(ti-mu)**2 for ti in diff_t]) / v
    sigma = math.sqrt(sigma)

    # 8. Compute the lower-bound of the confidence interval for the mean
    # based on a normal distribution with confidence level alpha=0.95
    mu_bar = mu -(2.576 * sigma / math.sqrt(v))

    # 9. Using a binary search, solve for the parameter p s.t. Ep equals mu_bar
    valid, p = solve_for_p(mu_bar, n)

    # 10. The min-entropy is the negative logarithm of the parameter p:
    # min-entropy = -log2(p)
     # If the search does not yield a solution, then estimate max min-entropy  
    if not valid:
        # No solution to equation. Assume max min-entropy.
        return 1.0/k, math.log(k,2)
    else:
        min_entropy = -math.log(p, 2.0)
        return p, min_entropy
