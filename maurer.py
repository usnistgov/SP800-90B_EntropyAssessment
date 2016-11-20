# Non-iid compression test using the Maurer universal statistic as defined
# in Draft SP 800-90B (January 2016) Section 6.3.4
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# T. A. Hall
# tim.hall@nist.gov
# 8 September 2014
#
# Updated by Kerry McKay
# March 3, 2016

from math import sqrt
from math import log
from math import floor

import numpy as np
import sys

from is_close import isclose
#from multiprocessing import Pool


d = 1000

# G(p) as defined in Step 8 of Section 9.3.6.3
def G(p, v):
    Gsum = 0.0
    N = v + d
    # The inner sum for s = 1 -> d...compute once...use in every iteration of outer loop
    inSum = p * p * sum([log(s, 2.0) * pow(1.0-p, s-1) for s in range(1,(d+1))])
    Gsum = inSum * (N - float(d))
    
    # log2(s) * (1-p)^(s-1) used in next two terms
    st = [log(s, 2.0) * pow(1.0-p, s-1) for s in range((d+1), N+1)]

    # The additional 's < t' sum term for s = (d+1) -> N-1
    Gsum += p * p * sum([(N-i-(d+1)) * st[i] for i in range(N-(d+1))])
    
    # The 's = t' term for s = (d+1) -> N
    Gsum += p * sum(st)

    return Gsum / v



# Expected value of Maurer Statistic using the one-param family of probability
# distributions in Step 7. of 9.3.6.3
def EppM(p, n, v):
    q = (1.0 - p)/(n - 1.0)
    return G(p, v) + (n - 1) * G(q, v)


# Binary search algorithm to find value of p s.t. the expected value
# of the Maurer Universal Statistic equals mu_bar within tolerance
#presumes a decreasing function
def solve_for_p(mu_bar, n, v, tolerance=1e-09):
    assert n > 0

    #This is a hackish way of checking to see if the difference is within approximately 4 ULPs
    absEpsilon = 4.0 * max((np.nextafter(mu_bar, mu_bar+1.0) - mu_bar), (mu_bar - np.nextafter(mu_bar, mu_bar-1.0)))
    #If we don't have numpy, then this will work for most of the ranges we're concerned with
    #absEpsilon = sys.float_info.epsilon

    if mu_bar > EppM(1.0/float(n), n, v):
        return False, 0.0

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

    centerVal = EppM(center, n, v)

    for rounds in range(1076):
        if isclose(mu_bar, centerVal, tolerance, absEpsilon):
            return True, center;

        if lbound >= hbound:
            print ("Bounds have converged after %d rounds and target was not found" % rounds)
            return False, 0.0

        if mu_bar < centerVal:
            lbound = center
            lvalue = centerVal
        else:
            hbound = center
            hvalue = centerVal

        if (mu_bar > lvalue) or (mu_bar < hvalue):
            print ("Target is not within the search interval after %d rounds" % rounds)
            return False, 0.0

        center = (lbound + hbound) / 2.0

        if (center <= ldomain) or (center >= hdomain):
            print ("The next center is outside of the proscribed domain after %d rounds" % rounds)
            return False, 0.0
        
        centerVal = EppM(center, n, v)

    if isclose(mu_bar, centerVal, tolerance, absEpsilon):
        return True, p
    else:
        return False, 0.0

# Maurer Universal Statistic compression test

# NOTE on index values - SP 800-90B Sec. 9.3.6.x uses one-based index values,
# whereas Python uses zero-based index values.  Therefore, I always add one
# to the index value returned from enumerate.  It only affects the result if an
# observation value is encountered once in the test data subset but not in the
# 1,000 sample dictionary subset.  In this case, an Ai value is an absolute
# index value and not a difference.

# Sec. 6.3.4 The Compression Estimate
def maurer_universal_statistic(dataset, k, verbose=False):
    # 1. Partition the dataset into two disjoint groups. These two groups
    #    will form the dictionary and the test data.
    # 1.a. Create the dictionary from the first d=1000 observations,
    # 1.b. Use the remaining v = L - d observations for testing
    #
    # 2. Initialize the dictionary, dict to an all zero array of size k.
    #    For i from 1 to d, let dict[s_i]=i
    v = len(dataset) - d
    dictionary = [0 for i in range(k)]
    for index, obs in enumerate(dataset[:d]):
        dictionary[obs] = index + 1
     
    # 3. Run the test data against the dictionary created in Step 2
    # 3.a. Let Di be a list of length v. 
    # 3.b. For i from d+1 to L:
    # 3.b.i If dict[si] is non-zero, then Di-d = i - dict[si]. Update the
    #       dictionary with the index of the most recent observation,
    #       dict[si] = i.
    # 3.b.ii If dict[si] is zero, add that value to the dictionary, i.e.,
    #       dict[si]=i. Let Di-d = i. 
    D = []
    for index, obs in enumerate(dataset[d:]):
        D.append(index + (d+1) - dictionary[obs])
        dictionary[obs] = index + (d+1)

    # 4. Let b = the number of bits needed to represent the largest symbol
    #    in the out alphabet. Calculate the sample mean, mu, and the sample
    #    standard dev, sigma, of (log2(Di)) for i = 1 -> v 
    #
    #    where the Di are the calculated differences from Step 3b
    b = floor(log(max(dataset),2.0))+1
    D = [log(Di, 2.0) for Di in D]
    mu = sum(D) / v
    #compute the variance of the dataset
    c = 0.7 - 0.8/b + (4+32.0/b)*(v**(-3.0/b))/15
    sigma = sum([Di * Di for Di in D])/v - (mu * mu)
    sigma = sqrt(sigma)
    sigma = c*sigma

    # 5. Compute the lower-bound of the 99% confidence interval for the mean 
    #    based on a normal distribution (student-T distribution is the nominal
    #    distribution, but the sample size is expected to be large.)
    mu_bar = mu - (2.5758293035489008*sigma)/sqrt(v)
    if verbose: 
        print("\tmu=%g, sigma=%g, mu_bar=%g" % (mu, sigma, mu_bar))
    
    # 6. Using a binary search, solve for parameter p
    valid, p = solve_for_p(mu_bar, k, v)
    if not valid:
        # No solution to equation. Assume max min-entropy.
        return 1.0/k, log(k,2.0)

    # 9. The min-entropy is the negative logarithm of the parameter, p
    min_entropy = -log(p, 2.0)

    return p, min_entropy

