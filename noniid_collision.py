import math

# The Collision Test
# Section 9.3.3 of DRAFT NIST SP 800-90B (August 2012)
#
# T. A. Hall
# CSD/ITL/NIST
# 25 Feb 2013


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
        Ep = calcEpS(p_c, n)
        #print('\tp = %g, Ep = %g, mu_bar = %g' % (p_c, Ep, mu_bar))
    return True, p_c

# Section 9.3.3.3 - Collision Test Details
def collision_test(dataset, n):
    # 1. Beginning with s1, step through the dataset until any observed
    # value is repeated; i.e., find smallest j s.t. si = sj for some
    # i with 1 <= i < j
    # 2. Define a sequence t. Set t0 = 0, t1 = j
    # 3. Set v = 1
    # 4. Starting with sample j, step through the remaining samples,
    # sj+1,sJ=2,...,sN, until there is a k and an l s.t. sk = sl,
    # with j < k < l
    # 5. v = v + 1
    # 6. Define tv = l (ell)
    # 7. Set j = l (ell)
    # 8. Repeat steps 4-7 until the end of the dataset is reached. A sequence
    # {t0,t1,...,tv} is generated.  If v < 1000, the noise source outputs will
    # be mapped down based on the ranking provided, and the data will be
    # retested.

    j = 0
    t = [0]
    for ell, s in enumerate(dataset):
        if s in dataset[j:ell]:
            t.append(ell + 1) # account for index starting at 1
            j = ell + 1

    # 9. Calculate sample mean, mu, and sample stddev, sigma of the differences
    # of collision times
    diff_t = [t[i]-t[i-1] for i in range(1,len(t))]
    v = float(len(diff_t)) # float() so it works in Python 2.6,2.7
    #print("\tv =", v)
    mu = sum(diff_t)/v
    sigma = sum([ti*ti for ti in diff_t]) / v - (mu * mu)
    sigma = math.sqrt(sigma)

    # 10. Compute the lower-bound of the confidence interval for the mean
    # based on a normal distribution with confidence level alpha=0.95
    mu_bar = mu -(1.96 * sigma / math.sqrt(v))

    # 11. Define a one-parameter family of prob. distributions parameterized by
    # p ...
    # 12. Using a binary search, solve for the parameter p s.t. Ep equals mu_bar
    valid, p = solve_for_p(mu_bar, n)
    if not valid:
        return False, 0.0, 0.0

    # 13. The min-entropy is the negative logarithm of the parameter p:
    # min-entropy = -log2(p)
    min_entropy = -math.log(p, 2.0)

    return True, p, min_entropy
