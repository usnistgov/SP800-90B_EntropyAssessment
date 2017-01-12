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

# Used in Step 9 in Section 6.3.2
#revised to be more consistent with the standard function in 8.9.2 here: http://dlmf.nist.gov/8.9
#Note that if we call the continued fraction from the DLMF as L(a, z), we are interested in
#F(k, z) = L(k+1, 1/z). Thus (in the notation of http://dlmf.nist.gov/3.10)
#b_0 = 0, b_j = 1 (j>0)
#a_1 = z, a_2 = (-k) z, a_3 = z, a_4 = (1-k)z, a_5 = 2z, etc.
def F(k, z):

    #first, correct "a" so that it means the same thing as in 8.9.2
    a = k+1

    tiny = 1E-30
    epsilon = sys.float_info.epsilon
    maxiter = 500

    #Now apply the modified Lentz's algorithm (see Numerical Recipes, section 5.2)
    #f_0 = tiny, as b_0 = 0
    #C_0 = tiny
    #D_0 = 0

    #j = 1 
    #a_1 = z
    #D_1 = 1 + a_1 D_0 = 1.0
    #C_1 = 1 + a_1 / C_0 (or tiny if this yields 0)
    D = 1.0
    C = 1.0 + z/tiny
    if C == 0.0:
        C = tiny

    delta = C * D
    f = tiny*delta

    for j in range(1, maxiter):
        #2j
        aj = (j-a)*z
        D = 1.0 + aj * D
        if D == 0.0:
            D = tiny
        C = 1.0 + aj/C
        if C == 0.0:
            C = tiny
        D = 1.0 / D
        delta = C * D
        f = f*delta

        #2j+1
        aj = j*z
        D = 1.0 + aj * D
        if D == 0.0:
            D = tiny
        C = 1.0 + aj/C
        if C == 0.0:
            C = tiny
        D = 1.0 / D
        delta = C * D
        f = f*delta

        if abs(delta - 1.0) < epsilon:
            #print('F: return %.17g after %d rounds' % (f, 2*j+1))
            return f

    print('Fell through after maximum iterations')
    return f

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
    #absEpsilon = 4.0 * sys.float_info.epsilon

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

        if mu_bar < centerVal:
            lbound = center
            lvalue = centerVal
        else:
            hbound = center
            hvalue = centerVal

#We now verify that ldomain <= lbound < center < hbound <= hdomain
#and that target in [ hvalue, lvalue ]
        if lbound >= hbound:
            print ("Bounds have converged after %d rounds and target was not found. Returning largest bound." % rounds)
            return True, min(max(lbound, hbound), hdomain)

        if (lbound < ldomain) or (lbound > hdomain) or (hbound < ldomain) or (hbound > hdomain):
            print ("The current search interval is not a subset of the domain after %d rounds and target was not found." % rounds)
            return False, 0.0

        if (mu_bar > lvalue) or (mu_bar < hvalue):
            print ("Target is not within the search interval after %d rounds" % j)
            return False, 0.0

        lastCenter = center
        center = (lbound + hbound) / 2.0

        if (center <= lbound) or (center >= hbound):
            print ("The next center is outside of the search interval after %d rounds" % j)
            return False, 0.0

        if lastCenter == center:
            print ("Detected cycle after %d rounds. Returning upper bound." % j)
            return True, hbound

        centerVal = calcEpS(center, n)

        #invariant: if this isn't true, then this isn't loosely monotonic
        if (centerVal < hvalue) or (centerVal > lvalue):
            print ("CenterVal is not within the search value interval after %d rounds. Returning upper bound." % j)
            return True, hbound

    if isclose(mu_bar, centerVal, tolerance, absEpsilon):
        return True, p
    else:
        print ("Binary search failed to converge. Returning upper bound.")
        return True, min(hbound, hdomain)

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
    mu_bar = mu -(2.5758293035489008 * sigma / math.sqrt(v))

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
