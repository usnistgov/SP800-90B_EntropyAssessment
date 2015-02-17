import math

# The Partial Collection Test
# Section 9.3.4 of Draft NIST SP 800-90B (August 2012), 
#
# T. A. Hall
# CSD/ITL/NIST
# 25 Feb 2013
# Updated 30 May 2013


# Compute expected value of Ti
def Ep_pc(p, n):
    # Distribution given Step 6
    q = (1.0 - p)/(n - 1)
    # Formulation of expected value given in Step 7
    return 1.0 - pow((1.0 - p), n) + (n - 1) * (1.0 - pow(1.0 - q, n))

# Binary search for value of 'p' s.t. E{Ti} = mu_bar
def solve_for_p(mu_bar, n, tolerance=.00001):
    p = 0.5
    adj = 0.5
    Ep = Ep_pc(p, n)
    Ep_maxvalid = Ep_pc(1.0/float(n), n)
    if mu_bar > Ep_maxvalid:
        print("\tpartial collection test not run:")
        print("\tmu_bar = %g, max valid value for this test and model = %g" % (mu_bar, Ep_maxvalid))
        return False, 0.0
    while abs(mu_bar - Ep) > tolerance:
        adj /= 2.0
        if mu_bar < Ep:
            p += adj
        else:
            p -= adj
        Ep = Ep_pc(p, n)

    return True, p

# Section 9.3.4.3 - Partial Collection Test Details
def partial_collection_test(dataset, n):
    # 1. Consider dataset as v non-overlapping data subsets of length n, where
    # n is the size of the output space
    # 2. Count the number of distinct values seen in each data subset (Ti)
    # 3. Repeat Step 2 until end of dataset is reached. If a minimum of 500
    # events have not been observed, the noise source outputs will be mapped
    # down (NOTE: this check and mapping are not performed inside this
    # function)
    v = len(dataset) // n
    counts = [len(set(dataset[i:i+n])) for i in range(0, v*n, n)]

    # 4. Calculate the sample mean, mu, and the sample std dev, sigma
    mu = sum(counts) / float(v)
    
    sigma = sum([Ti * Ti for Ti in counts]) / float(v) - (mu * mu)
    sigma = math.sqrt(sigma)

    # 5. Lower bound of confidence interval for the mean
    mu_bar = mu - (1.96 * sigma)/math.sqrt(v)

    # 6. Define one-parameter family of prob distributions...
    # 7. Solve for p s.t. E{Ti} equals mu_bar
    valid, p = solve_for_p(mu_bar, n)
    if not valid:
        return False, 0.0, 0.0

    # 8. The min-entropy is negative log base 2 of p
    min_entropy = -math.log(p, 2.0)

    return True, p, min_entropy
