import math

# The Collision Estimate
# Section 6.3.2 of DRAFT NIST SP 800-90B (2015)
#
# T. A. Hall
# CSD/ITL/NIST
# 25 Feb 2013
#
# Updated by Kerry McKay, 2015


# Continued fraction for F(1/z) derived from Eq. 8.9.2 at http://dlmf.nist.gov/8.9
# Used in Step 12 in Section 9.3.3.3
def F(n, zInv):
    z = 1.0/zInv
    denom = 1.0 + n/z

    for i in range(1, n):
        denom = z + -i/denom
        denom = 1.0 + (n-i)/denom

    denom = z + -n/denom
    return 1.0/denom


# Expected value of statistic based on one-parameter family of prob distributions
# Steps 11, 12 of Section 9.3.3.3
def calcEpS(p, n):
    q = (1.0 - p)/float(n - 1)
    i_n = 1.0 / n
    ip = 1.0/p
    iq = 1.0/q
    iq2 = iq * iq
    Ep = p * iq2 * (1.0 + i_n * (ip - iq)) * F(n, q) - (p * iq * i_n * (ip - iq))
    return Ep


# Binary search for p s.t. Ep(p) equals mu_bar
def solve_for_p(mu_bar, n):
    p_c = 0.5
    adj = 0.5
    Ep = calcEpS(p_c, n)
    Ep_maxvalid = calcEpS(1.0/float(n),n)
    minp = 1.0/float(n) 
    if mu_bar > Ep_maxvalid:
        print("\tcollision test not run:")
        print("\tmu_bar = %g, max valid value for this test and model = %g" % (mu_bar, Ep_maxvalid))
        return False, 0.0
    while abs(mu_bar - Ep) > .0001:
        adj /= 2.0
        if mu_bar < Ep:
            p_c += adj
        else:
            p_c -= adj
            #occasionally dips below lowest possible pmax. This is to fix that
            if p_c < minp:
                p_c = minp
        Ep = calcEpS(p_c, n)
        #print('\tp = %g, Ep = %g, mu_bar = %g' % (p_c, Ep, mu_bar))
    return True, p_c

# Section 9.3.3.3 - Collision Test Details
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
        print "Must map down data for collision estimate, then retest"
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
    min_entropy = -math.log(p, 2.0)
    # If the search does not yield a solution, then estimate max min-entropy  
    if not valid:
        # No solution to equation. Assume max min-entropy.
        return 1.0/k, math.log(k,2)

    return p, min_entropy
